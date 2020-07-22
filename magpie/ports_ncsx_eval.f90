!-------------------------------------------------------------------------------
! ncsx_ports_eval.f90
!
! Subroutines to evaluate whether user-provided points are inside NCSX ports.
! In the normal usage case, these will be called many times during a program
! run.
!-------------------------------------------------------------------------------

module ports_ncsx_eval

use magpie_globals, only: dp

implicit none

contains

!-------------------------------------------------------------------------------
! in_ncsx_port(x, y, z, in_port)
!
! Determines whether a user-supplied point is inside an NCSX port
!
! Input parameters:
!     real(dp)    :: x, y, z  -> the x, y, and z coordinate of the point to be 
!                            assessed
!
! Output parameters:
!     logical :: in_port  -> true if (x, y, z) is in any of the NCSX ports,
!                            except those ports explicitly declared to be 
!                            excluded in the global variables.
!
! Required modules:
!     ports_ncsx_globals
!     ncsx_ports_load
!
! Subroutines:
!     in_circular_ports()
!
! Functions:
!     in_port_4()
!     in_port_12()
!     in_nb_port()
!     in_dome()
!-------------------------------------------------------------------------------
subroutine in_ncsx_port(x, y, z, in_port)

    use magpie_globals,     only: incl_port_ncsx_nb, incl_port_ncsx_04,      &
                                  incl_port_ncsx_12, incl_port_ncsx_dome
    use ports_ncsx_globals, only: ncsx_ports_on, ncsx_ports_initialized, &
                                  nfp_ncsx
    use ports_ncsx_load,    only: load_ncsx_ports

    implicit none

    real(dp),    intent(IN)  :: x, y, z
    logical, intent(OUT) :: in_port
    real(dp)                 :: x_sym, y_sym, z_sym
    logical              :: in_cir = .false.,   in_p04 = .false.,   &
                            in_p12 = .false.,   in_pnb = .false.,   &
                            in_dom = .false.

    ! Load the port data into memory if necessary
    if (.not. ncsx_ports_initialized) then
        call load_ncsx_ports()
        ncsx_ports_initialized = .true.
    end if

    ! Collapse the point to its stellarator-symmetric equivalent
    call stellarator_symmetric_equiv(x, y, z, nfp_ncsx, x_sym, y_sym, z_sym)

    ! Check whether point is inside each port or port set if applicable
    call in_circular_ports(x_sym, y_sym, z_sym, in_cir)
    if (incl_port_ncsx_04)   in_p04 = in_port_4(x_sym, y_sym, z_sym)
    if (incl_port_ncsx_12)   in_p12 = in_port_12(x_sym, y_sym, z_sym)
    if (incl_port_ncsx_nb)   in_pnb = in_nb_port(x_sym, y_sym, z_sym)
    if (incl_port_ncsx_dome) in_dom = in_dome(x_sym, y_sym, z_sym)

    in_port = in_cir .or. in_p04 .or. in_p12 .or. in_pnb .or. in_dom

end subroutine in_ncsx_port

!-------------------------------------------------------------------------------
! stellarator_symmetric_equiv(x, y, z, nfp, x_sym, y_sym, z_sym)
!
! Returns the stellarator-symmetric equivalent coordinates in the first 
! half-period for an input set of coordinates.
!
! Input parameters:
!     real(dp) :: x, y, z     -> Cartesian x, y, and z coordinates 
!     integer :: nfp      -> Number of field periods
!
! Output parameters:
!     real(dp) :: x_sym, y_sym, z_sym  -> Equivalent point in first half-period
!-------------------------------------------------------------------------------
subroutine stellarator_symmetric_equiv(x, y, z, nfp, x_sym, y_sym, z_sym)

    use ports_ncsx_globals, only: pi

    implicit none

    real(dp), intent(IN)    :: x, y, z
    integer, intent(IN) :: nfp
    real(dp), intent(OUT)   :: x_sym, y_sym, z_sym
    real(dp)                :: r, phi, phi_interval, phi_halfInt, period

    ! Calculate toroidal angle and ensure it is in the interval [0, 2*pi)
    phi = atan2(y, x)
    if (phi < 0) phi = phi + 2.0 * pi

    ! Collapse to a single field period
    phi_interval = 2.0 * pi / nfp
    period = floor(phi / phi_interval)
    phi = phi - phi_interval * period

    ! Reflect point if in 2nd half-period
    phi_halfInt  = phi_interval / 2.0
    if (phi > phi_halfInt) then
        phi = phi_interval - phi
        z_sym = -z
    else
        z_sym = z
    end if

    ! Calculate equivalent x and y coordinates
    r = sqrt(x*x + y*y)
    x_sym = r * cos(phi)
    y_sym = r * sin(phi)

    return

end subroutine stellarator_symmetric_equiv

!-------------------------------------------------------------------------------
! in_circular_ports(x, y, z, in_cir)
!
! Wrapper that calls in_circular_port() for each of the included circular
! ports
!-------------------------------------------------------------------------------
subroutine in_circular_ports(x, y, z, in_cir)

    use ports_ncsx_globals, only: incl_circular_ports, nCirPrt, id_cir, &
                              xo_cir, yo_cir, zo_cir, ax_cir, ay_cir, az_cir, &
                              ir_cir, thick_cir, l0_cir, l1_cir

    implicit none

    real(dp),    intent(IN)  :: x, y, z
    logical, intent(OUT) :: in_cir
    logical              :: in_cir_i
    integer              :: i

    in_cir = .false.

    do i = 1, nCirPrt
        if (incl_circular_ports(i)) then
            in_cir_i = in_circular_port(x, y, z, &
                                        xo_cir(i), yo_cir(i), zo_cir(i), &
                                        ax_cir(i), ay_cir(i), az_cir(i), &
                                        ir_cir(i), thick_cir(i), &
                                        l0_cir(i), l1_cir(i))
            in_cir = in_cir .or. in_cir_i
        end if
    end do

end subroutine in_circular_ports

!-------------------------------------------------------------------------------
! in_circular_port(x, y, z, xo, yo, zo, ax, ay, az, ir, thick, l0, l1)
!
! Determines whether a given point in 3D space is inside a circular port
! of given parameters.
!
! Input parmaters:
!     real(dp) :: x, y, z    -> Cartesian x, y, and z coordinates of the point to
!                           be tested.
!     real(dp) :: xo, yo, zo -> x, y, and z coordinates of the port origin; i.e., a
!                           reference point through which the port axis passes
!     real(dp) :: ax, ay, az -> x, y, and z components of a unit vector in the
!                           direction of the axis
!     real(dp) :: ir         -> Inner radius of the port
!     real(dp) :: thick      -> Thickness of the port wall
!     real(dp) :: l0, l1     -> Distances along the axis from the origin to the
!                           beginning and end of the exclusion volume, 
!                           respectively
!
! Return value: logical indicating whether x, y, and z are in the exclusion
!               volume of the port
!-------------------------------------------------------------------------------
logical function in_circular_port(x, y, z, xo, yo, zo, ax, ay, az, & 
                                  ir, thick, l0, l1)

    implicit none

    ! Input parameters (see documentation string for more information)
    real(dp), intent(IN) :: x, y, z, xo, yo, zo, ax, ay, az, ir, thick, l0, l1

    real(dp) :: l_proj
    real(dp) :: proj_x, proj_y, proj_z
    real(dp) :: rvec_x, rvec_y, rvec_z, r2, or, or2

    ! Initialize return value
    in_circular_port = .false.

    ! Project test points onto the port axis
    l_proj = (x - xo) * ax + (y - yo) * ay + (z - zo) * az
    if (l_proj < l0 .or. l_proj > l1) return  ! Return if outside axial bounds
    proj_x = l_proj * ax + xo
    proj_y = l_proj * ay + yo
    proj_z = l_proj * az + zo

    ! Radial distance from the axis
    rvec_x = x - proj_x
    rvec_y = y - proj_y
    rvec_z = z - proj_z
    r2 = rvec_x*rvec_x + rvec_y*rvec_y + rvec_z*rvec_z

    ! Compare distance to port outer radius
    or = ir + thick
    or2 = or*or
    in_circular_port = (r2 <= or2)

    return

end function in_circular_port

!-------------------------------------------------------------------------------
! in_dome(x, y, z)
!
! Determines whether a test point is in the dome (which is modeled here as a
! circular port). Considers the domes in both the first and second half-modules,
! as the dome overlaps the triangular cross-section (the boundary between
! half-modules).
!-------------------------------------------------------------------------------
logical function in_dome(x, y, z)

    use ports_ncsx_globals, &
        only: xo_dom,  yo_dom,  zo_dom,  ax_dom,  ay_dom,  az_dom,  &
              xo_dom2, yo_dom2, zo_dom2, ax_dom2, ay_dom2, az_dom2, &
              ir_dom, thick_dom, l0_dom, l1_dom

    implicit none

    real(dp), intent(IN) :: x, y, z

    in_dome = (in_circular_port(x, y, z, xo_dom, yo_dom, zo_dom,          &
                               ax_dom, ay_dom, az_dom, ir_dom,            &
                               thick_dom, l0_dom, l1_dom)                 &
                  .or.                                                    &
               in_circular_port(x, y, z, xo_dom2, yo_dom2, zo_dom2,       &
                                ax_dom2, ay_dom2, az_dom2, ir_dom,        &
                                thick_dom, l0_dom, l1_dom)                 )

    return

end function in_dome

!-------------------------------------------------------------------------------
! in_nb_port(x, y, z)
! 
! Determines whether a test point with coordinates x, y, and z is inside the
! volume of the NB port.
!-------------------------------------------------------------------------------
logical function in_nb_port(x, y, z)

    use ports_ncsx_globals, only: yo1_pnb, zo1_pnb, yo2_pnb, zo2_pnb, &
                              yo3_pnb, zo3_pnb, yo6_pnb, zo6_pnb, &
                              or1sq_pnb, or2sq_pnb, or3sq_pnb, or6sq_pnb, &
                              l0_pnb, l1_pnb, zMax_pnb, yMax_pnb, &
                              yTop_pnb, zTop_pnb, yBot_pnb, zBot_pnb, &
                              yUp_pnb,  zUp_pnb,  yLo_pnb,  zLo_pnb, &
                              slopeUp_pnb, slopeLo_pnb

    implicit none

    real(dp), intent(IN) :: x, y, z      ! Coordinates of test point

    ! Initialize return value
    in_nb_port = .false.

    ! Return false if not within a bounding box
    if (x < l0_pnb    .or.       x > l1_pnb    .or. &
        y > yMax_pnb  .or.  abs(z) > zMax_pnb        )  return

    !---------------------------------------------------------------------------
    ! If in the bounding box, evaluate depending on the z coordinate:
    !---------------------------------------------------------------------------

    ! Between edge of circle 6 and the xz plane:
    if (z < zBot_pnb) then
        in_nb_port = (y - yo6_pnb <= sqrt(or6sq_pnb - (z-zo6_pnb)**2 )) 

    ! Between the lower segment and the xz plane:
    else if (z < zLo_pnb) then
        in_nb_port = (y - yLo_pnb <= slopeLo_pnb * (z - zLo_pnb)) 

    ! Between the edge of circle 1 and the xz plane:
    else if (z < zUp_pnb) then
        in_nb_port = (y - yo1_pnb <= sqrt(or1sq_pnb - (z-zo1_pnb)**2 )) 

    ! Between the upper segment and the xz plane:
    else if (z < zTop_pnb) then
        in_nb_port = (y - yUp_pnb <= slopeUp_pnb * (z - zUp_pnb))

    ! Between the edges of circle 2 (upper) and circle 3 (lower):
    else
        in_nb_port = (y - yo2_pnb <=  sqrt(or2sq_pnb - (z-zo2_pnb)**2) .and. &
                      y - yo3_pnb >= -sqrt(or3sq_pnb - (z-zo3_pnb)**2)        )
    end if

    return

end function in_nb_port

!-------------------------------------------------------------------------------
! in_port_4(x, y, z)
!
! Determines whether an input point is inside the volume of NCSX port 4
!-------------------------------------------------------------------------------
logical function in_port_4(x, y, z)

    use ports_ncsx_globals, only: zMax_p04, z_circ_p04, l0_p04, l1_p04, phi_p04, &
                              theta_a_p04, theta_b_p04, &
                              w1a_p04, w1b_p04, w2a_p04, w2b_p04, &
                              slope_a_p04, slope_b_p04, &
                              l_narrow_p04, orc_p04, orc2_p04, &
                              thick_ai_p04, thick_bi_p04, &
                              thick_ao_p04, thick_bo_p04, &
                              rhat_x_p04, rhat_y_p04, phat_x_p04, phat_y_p04

    implicit none

    real(dp), intent(IN) :: x, y, z
    real(dp)             :: abs_z, r, w, wa, wai, wao, wb, wbi, wbo

    ! Initialize output value
    in_port_4 = .false.

    ! Radial coordinate of the point's orthog. projection onto the port axis
    r = x * rhat_x_p04 + y * rhat_y_p04

    ! Ensure input point is within z bounds; otherwise terminate
    abs_z = abs(z)
    if (abs_z > zMax_p04 .or. r < l0_p04 .or. r > l1_p04) return

    ! Distance from the port axis, parallel to the xy plane
    w = x * phat_x_p04 + y * phat_y_p04

    ! Constraints on w in the -phi (a) direction, dependent on r, 
    ! for "inner" and "outer" sections relative to the narrow point
    wai = w1a_p04 + (l_narrow_p04 - r) * tan(theta_a_p04) + thick_ai_p04
    wao = w1a_p04 + (r - l_narrow_p04) * slope_a_p04      + thick_ao_p04
    wa = max(wai, wao)   

    ! Constraint on w in the +phi (b) direction, dependent on r,
    ! for "inner" and "outer" sections relative to the narrow point
    wbi = w1b_p04 - (l_narrow_p04 - r) * tan(theta_b_p04) + thick_bi_p04
    wbo = w1b_p04 + (r - l_narrow_p04) * slope_b_p04      + thick_bo_p04
    wb = max(wbi, wbo)  

    ! Evaluate depending on whether between straight walls or rounded corners
    if (abs_z <= z_circ_p04) then
        in_port_4 = (w >= -wa .and. w <= wb)
    else
        in_port_4 = (w >= -wa + orc_p04 &
                          - sqrt(orc2_p04 - (abs_z-z_circ_p04)**2) .and. &
                     w <=  wb - orc_p04 &
                          + sqrt(orc2_p04 - (abs_z-z_circ_p04)**2)        )
    end if

    return

end function in_port_4

!-------------------------------------------------------------------------------
! in_port_12(x, y, z)
! 
! Determines whether the input point is within the exclusion volume for port 12
!-------------------------------------------------------------------------------
logical function in_port_12(x, y, z)

    use ports_ncsx_globals, only: xMin_p12, xMax_p12, yMax_p12, &
                              xo1_p12, xo2_p12, slope_p12, &
                              xInb_p12, yInb_p12, xOut_p12, &
                              or1_p12, or1sq_p12, or2_p12, or2sq_p12, &
                              l0_p12, l1_p12

    implicit none

    real(dp), intent(IN) :: x, y, z
    real(dp)             :: abs_z

    ! Initialize output value
    in_port_12 = .false.

    ! Determine whether the point is within a bounding box
    abs_z = abs(z)
    if (x < xMin_p12 .or. x > xMax_p12 .or. y > yMax_p12 .or. &
        abs_z < l0_p12 .or. abs_z > l1_p12  ) return

    !---------------------------------------------------------------------------
    ! Evaluate based on where the point falls along the x-axis
    !---------------------------------------------------------------------------

    ! Between the inboard circle and the x-axis
    if (x < xInb_p12) then
        in_port_12 = (y**2 <= or1sq_p12 - (x - xo1_p12)**2)

    ! Between the line segment and the x-axis
    else if (x <= xOut_p12) then
        in_port_12 = (y <= yInb_p12 + slope_p12 * (x - xInb_p12))

    ! Between the outboard circle and the x-axis
    else
        in_port_12 = (y**2 <= or2sq_p12 - (x - xo2_p12)**2)

    end if

    return

end function in_port_12

end module ports_ncsx_eval


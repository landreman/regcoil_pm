!-------------------------------------------------------------------------------
! qhex_array.f90
! 
! Module with subroutines used to construct a grid of magnets with quadrilateral
! hexahedron geometry around a toroidal surface.
!
! Author: K. C. Hammond
! Contact: khammond@pppl.gov
! Updated: 2020-02-13
!-------------------------------------------------------------------------------
module qhex_array

use magpie_globals, only: dp

implicit none

contains

!-------------------------------------------------------------------------------
! build_qhex_array(ref_surf, lim_surf, bg, incl_upper, adj_interval, qhArr)
!
! Constructs an array of quadrilaterally-faced hexahedrons according to a
! pre-defined background grid.
!
! Input parameters:
!     type(surface) :: ref_surf    -> The reference toward which qhexs "point"
!     type(surface) :: lim_surf    -> Upper limiting surface for qhexs (vessel)
!     type(base_grid_struct) :: bg -> Structure with data for the base grid
!     logical :: incl_upper        -> Whether to add upper qhexs (if possible)
!     real(dp) :: adj_interval     -> Interval for retracting faces when 
!                                     removing overlaps
!     
! Output parameters:
!     type(qhex_config) :: qhArr  -> Structure with info on the constructed
!                                     configuration
!-------------------------------------------------------------------------------
subroutine build_qhex_array(ref_surf, lim_surf, bg, incl_upper, &
                            max_ht_lo, max_ht_up, adj_interval, qhArr)

    use magpie_globals,  only: len_line
    use surface_calc,    only: surface !, surf_r, surf_z, surf_unorm
    use surface_solve,   only: surf_dist_3d
    use base_grid,       only: base_grid_struct
    use qhex_properties, only: quad_hexahedron, qhex_array, qhex_corners, &
                               qhex_volume_ctr, check_qhex
    use qhex_adjust,     only: reset_qhex_height

    implicit none

    type(surface), intent(IN) :: ref_surf, lim_surf
    type(base_grid_struct), intent(IN) :: bg
    logical, intent(IN) :: incl_upper
    real(dp), intent(IN) :: max_ht_lo, max_ht_up, adj_interval
    type(qhex_array), intent(OUT) :: qhArr
    integer :: i, j, nTheta, nPhi
    real(dp) :: r, rnx, rny, rnz, l0=0.0, l, chi2
    character(len=len_line) :: str_nErr

    nTheta = bg%nTheta - 1
    nPhi = bg%nPhi - 1

    ! Set the global properties for the qhex configuration
    qhArr%incl_upper = incl_upper
    qhArr%nErr = 0
    qhArr%nErrLo = 0
    qhArr%nErrUp = 0
    qhArr%nTheta = nTheta
    qhArr%nPhi = nPhi
    qhArr%nLower = nTheta*nPhi
    allocate(qhArr%pairs(nTheta, nPhi))
    if (incl_upper) then
        qhArr%nTotal = 2*qhArr%nLower
        qhArr%nUpper = qhArr%nLower
    else
        qhArr%nTotal = qhArr%nLower
        qhArr%nUpper = 0
    end if

    write(*,fmt='(A)') '    Constructing array of quadrilateral-faced ' &
                       // 'hexahedra'

    do i = 1, nTheta
        do j = 1, nPhi

            qhArr%pairs(i,j)%incl_upper = incl_upper

            ! Define lateral bounding planes 
            call bounding_planes_lateral(i, j, bg, qhArr%pairs(i,j)%lo)

            ! Determine origin point for qhex normal on the reference surface
            call qhex_axis(ref_surf, bg%thetaSurf(i,j), bg%phiSurf(i,j), &
                           qhArr%pairs(i,j)%lo)

            ! Intersection of the qhex normal vector with the vessel surface
            call surf_dist_3d(lim_surf,                                      &
                     qhArr%pairs(i,j)%lo%rx, qhArr%pairs(i,j)%lo%ry,         &
                     qhArr%pairs(i,j)%lo%rz, -qhArr%pairs(i,j)%lo%snx,       &
                     -qhArr%pairs(i,j)%lo%sny, -qhArr%pairs(i,j)%lo%snz,     &
                     l0, qhArr%pairs(i,j)%lo%ref_theta,                      &
                     qhArr%pairs(i,j)%lo%ref_phi,                            &
                     l, qhArr%pairs(i,j)%lo%ves_theta,                       &
                     qhArr%pairs(i,j)%lo%ves_phi,                            &
                     qhArr%pairs(i,j)%lo%vx, qhArr%pairs(i,j)%lo%vy,         &
                     qhArr%pairs(i,j)%lo%vz, chi2                             )

            ! Determine the base bounding plane
            call base_plane(i, j, bg, qhArr%pairs(i,j)%lo)

            ! Use reference surface point as initial top plane reference point
            qhArr%pairs(i,j)%lo%otx = qhArr%pairs(i,j)%lo%rx
            qhArr%pairs(i,j)%lo%oty = qhArr%pairs(i,j)%lo%ry
            qhArr%pairs(i,j)%lo%otz = qhArr%pairs(i,j)%lo%rz

            ! Determine the corners of the qhex using initial bounding planes
            call qhex_corners(qhArr%pairs(i,j)%lo)

            ! Reset the height according to the lateral edges
            call reset_qhex_height(lim_surf, qhArr%pairs(i,j)%lo, max_ht_lo)

            ! Calculate the volume and centroid of the qhex
            call qhex_volume_ctr(qhArr%pairs(i,j)%lo)

            ! Check for erroneous properties
            call check_qhex(qhArr%pairs(i,j)%lo)
            if (qhArr%pairs(i,j)%lo%err) qhArr%nErrLo = qharr%nErrLo + 1

            ! Add an upper qhex if applicable
            if (incl_upper) then
                call add_upper_qhex(lim_surf, qhArr%pairs(i,j)%lo, &
                                    max_ht_up, qhArr%pairs(i,j)%up  )
                if (qhArr%pairs(i,j)%up%err) qhArr%nErrUp = qhArr%nErrUp + 1
            end if

        end do
    end do

    qhArr%nErr = qhArr%nErrLo + qhArr%nErrUp


    ! Check for instances of overlap and shrink qhexes to remove
    write(*,fmt='(A)') '    Adjusting array to remove magnet-magnet overlaps'
    call remove_overlaps(qhArr, adj_interval, max_ht_up, lim_surf)

    ! Summary report
    str_nErr = ''
    if (incl_upper) write(str_nErr, fmt='(A, I0, A, I0, A)') &
        '(', qhArr%nErrLo, ' lower, ', qhArr%nErrUp, ' upper)'
    write(*,fmt='(A, I0, A, I0, A)') &
        '    Grid constructed with ', qhArr%nTotal, ' hexahedra with ', &
        qhArr%nErr, ' erroneous ' // trim(str_nErr)

end subroutine build_qhex_array

!-------------------------------------------------------------------------------
! split_qhex(qhIn, n, gap, qhOut)
!
! Slices a qhex into a set of qhexes of equal height.
!
! Input parameters:
!     type(quad_hexahedron) :: qhIn  
!                      -> the qhex to be "sliced" (note: the input structure is
!                         not modified!)
!     integer :: n     -> Number of slices
!     real(dp) :: gap  -> Gap spacing to include between slices. Note that the
!                          slices will be erroneous if gap*(n-1) > qhIn%height.
!
! Output parameters:
!     type(quad_hexahedron), dimension(:) :: qhOut
!                      -> Array (must be length n) of slices
!-------------------------------------------------------------------------------
subroutine split_qhex(qhIn, n, gap, qhOut)

    use qhex_properties, only: quad_hexahedron, qhex_corners, qhex_volume_ctr, &
                               check_qhex

    implicit none

    type(quad_hexahedron), intent(IN) :: qhIn
    integer, intent(IN) :: n
    real(dp), intent(IN) :: gap
    type(quad_hexahedron), dimension(:), intent(OUT) :: qhOut
    real(dp) :: sub_height, h_bot, h_top
    integer :: i

    sub_height = qhIn%height / n

    do i = 1, n

        ! Initialize each output qhex as the input qhex
        qhOut(i) = qhIn

        ! Modify top and base planes
        h_bot = sub_height*(i-1) + 0.5*gap
        h_top = sub_height*i     - 0.5*gap
        qhOut(i)%obx = qhIn%obx + qhIn%snx * h_bot
        qhOut(i)%oby = qhIn%oby + qhIn%sny * h_bot
        qhOut(i)%obz = qhIn%obz + qhIn%snz * h_bot
        qhOut(i)%otx = qhIn%obx + qhIn%snx * h_top
        qhOut(i)%oty = qhIn%oby + qhIn%sny * h_top
        qhOut(i)%otz = qhIn%obz + qhIn%snz * h_top

        ! Update corner points + height, volume, and center accordingly
        call qhex_corners(qhOut(i))
        call qhex_volume_ctr(qhOut(i))
        call check_qhex(qhOut(i))

        ! If the parent qhex was erroneous, all children will remain as such
        if (qhIn%err) qhOut(i)%err = .true.

    end do

end subroutine split_qhex

!-------------------------------------------------------------------------------
! bounding_planes_lateral()
! 
! Determines parameters of the bounding planes for the toroidal and poloidal 
! sides of each trapezoidal box through calls to the plane_boxcar subroutine.
!
! Most of the work in this wrapper subroutine concerns selecting subsets of
! points in the vert_{x,y,z} and ves_{x,y,z} arrays to be used in the boxcar
! averaging, handling cases in which the input trapezoid is near the edges
! of the grid of planes and poloidal vertices.
!-------------------------------------------------------------------------------
subroutine bounding_planes_lateral(i_pol, j_tor, bg, qhex)

    use magpie_globals,  only: gap_qhx, pi, nfp, qhex_nBoxPol, qhex_nBoxTor, &
                                  qhex_poloidally_closed, &
                                  qhex_toroidally_closed
    use base_grid,       only: base_grid_struct
    use qhex_properties, only: quad_hexahedron

    implicit none

    integer, intent(IN) :: i_pol, j_tor
    type(base_grid_struct) :: bg
    type(quad_hexahedron), intent(INOUT) :: qhex
    integer :: n_arr_pol_max, n_arr_pol, n_arr_tor_max, n_arr_tor
    integer :: ind_min_pol, ind_max_pol, ind_min_tor, ind_max_tor
    integer :: ind_ref_pol, ind_ref_tor
    integer :: m, n, tor_index, disp, disp_min
    integer :: pol_index, pol_ind_offs 
    logical :: reflect
    real(dp) :: phiref, phi, zero = 0.0
    real(dp), allocatable :: arr_base_pb_x(:), arr_base_pb_y(:), arr_base_pb_z(:)
    real(dp), allocatable :: arr_base_pf_x(:), arr_base_pf_y(:), arr_base_pf_z(:)
    real(dp), allocatable :: arr_base_tb_x(:), arr_base_tb_y(:), arr_base_tb_z(:)
    real(dp), allocatable :: arr_base_tf_x(:), arr_base_tf_y(:), arr_base_tf_z(:)
    real(dp), allocatable :: arr_top_pb_x(:),  arr_top_pb_y(:),  arr_top_pb_z(:)
    real(dp), allocatable :: arr_top_pf_x(:),  arr_top_pf_y(:),  arr_top_pf_z(:)
    real(dp), allocatable :: arr_top_tb_x(:),  arr_top_tb_y(:),  arr_top_tb_z(:)
    real(dp), allocatable :: arr_top_tf_x(:),  arr_top_tf_y(:),  arr_top_tf_z(:)
    real(dp) :: ang_dx, ang_dy, ang_dz

    !---------------------------------------------------------------------------
    ! Part 0: Allocate arrays to store vertex points for plane boxcar averaging.
    !         Use the largest possible size for each array.
    !---------------------------------------------------------------------------

    n_arr_pol_max = 2*qhex_nBoxPol + 2
    n_arr_tor_max = 2*qhex_nBoxTor + 2

    allocate(arr_base_tb_x(n_arr_pol_max), arr_base_tb_y(n_arr_pol_max), &
             arr_base_tb_z(n_arr_pol_max), arr_base_tf_x(n_arr_pol_max), &
             arr_base_tf_y(n_arr_pol_max), arr_base_tf_z(n_arr_pol_max), &
             arr_top_tb_x(n_arr_pol_max), arr_top_tb_y(n_arr_pol_max),   &
             arr_top_tb_z(n_arr_pol_max), arr_top_tf_x(n_arr_pol_max),   &
             arr_top_tf_y(n_arr_pol_max), arr_top_tf_z(n_arr_pol_max),   &
             arr_base_pb_x(n_arr_tor_max), arr_base_pb_y(n_arr_tor_max), &
             arr_base_pb_z(n_arr_tor_max), arr_base_pf_x(n_arr_tor_max), &
             arr_base_pf_y(n_arr_tor_max), arr_base_pf_z(n_arr_tor_max), &
             arr_top_pb_x(n_arr_tor_max), arr_top_pb_y(n_arr_tor_max),   &
             arr_top_pb_z(n_arr_tor_max), arr_top_pf_x(n_arr_tor_max),   &
             arr_top_pf_y(n_arr_tor_max), arr_top_pf_z(n_arr_tor_max)     )

    !---------------------------------------------------------------------------
    ! Part 1: For the input trapezoid with a location on the grid defined
    !         by i_plane and i_pol, determine the number of points to 
    !         average over in the poloidal dimension (for the toroidal 
    !         back/front planes) and in the toroidal dimension (for the 
    !         poloidal back/front planes)
    !
    !         In cases where periodicity is assumed, the number of points will
    !         always be the same. Otherwise, the number of points is reduced
    !         for trapezoids near the edge of the grid (depending on the
    !         n_boxcar variables).
    !---------------------------------------------------------------------------
    if (qhex_poloidally_closed) then
        n_arr_pol = n_arr_pol_max
        ind_ref_pol = qhex_nBoxPol + 1
    else
        ind_min_pol = max(1, i_pol-qhex_nBoxPol)
        ind_max_pol = min(bg%nTheta, i_pol+1+qhex_nBoxPol)
        n_arr_pol = ind_max_pol - ind_min_pol + 1
        ind_ref_pol = i_pol - ind_min_pol + 1
    end if

    if (qhex_toroidally_closed) then
        n_arr_tor = n_arr_tor_max
        ind_ref_tor = qhex_nBoxTor + 1
    else
        ind_min_tor = max(1, j_tor-qhex_nBoxTor)
        ind_max_tor = min(bg%nPhi+1, j_tor+1+qhex_nBoxTor)
        n_arr_tor = ind_max_tor - ind_min_tor + 1
        ind_ref_tor = j_tor - ind_min_tor + 1
    end if

    !---------------------------------------------------------------------------
    ! Part 2: Populate the arrays of points to be used for boxcar averaging
    !         for the toroidal back/front and poloidal back/front planes.
    !
    !         In cases where periodicity is assumed, boxcar points that lie
    !         "over the edge of the grid" are taken from stellarator-symmetric
    !         locations in the grid. Otherwise, the arrays of points are 
    !         truncated at the edge of the grid.
    !---------------------------------------------------------------------------
    do m = 1, n_arr_pol

        ! Determine indices of the points to be taken from the vert/top arrays
        if (qhex_poloidally_closed) then

            disp = (m-1) - qhex_nBoxPol
            if (i_pol + disp <= 0) then
                pol_index = (bg%nTheta-1) + i_pol + disp
            else if (i_pol + disp > bg%nTheta) then
                pol_index = i_pol + disp - (bg%nTheta-1)
            else
                pol_index = i_pol + disp
            end if

            if (pol_index > bg%nTheta .or. pol_index <= 0) then
                stop 'bounding_planes_lateral: qhex_nBoxPol must be < bg%nTheta'
            end if

        else

            pol_index = ind_min_pol + m - 1

        end if

        ! Populate array parameters for the toroidal front and back planes
        arr_base_tb_x(m) = bg%x(pol_index, j_tor)
        arr_base_tb_y(m) = bg%y(pol_index, j_tor)
        arr_base_tb_z(m) = bg%z(pol_index, j_tor)
        arr_base_tf_x(m) = bg%x(pol_index, j_tor+1)
        arr_base_tf_y(m) = bg%y(pol_index, j_tor+1)
        arr_base_tf_z(m) = bg%z(pol_index, j_tor+1)
        arr_top_tb_x(m) = bg%xSurf(pol_index, j_tor)
        arr_top_tb_y(m) = bg%ySurf(pol_index, j_tor)
        arr_top_tb_z(m) = bg%zSurf(pol_index, j_tor)
        arr_top_tf_x(m) = bg%xSurf(pol_index, j_tor+1)
        arr_top_tf_y(m) = bg%ySurf(pol_index, j_tor+1)
        arr_top_tf_z(m) = bg%zSurf(pol_index, j_tor+1)

    end do

    do n = 1, n_arr_tor

        ! Determine indices of the points to be taken from the vert/top arrays
        if (qhex_toroidally_closed) then

            disp = (n-1) - qhex_nBoxTor
            if (j_tor + disp <= 0) then
                tor_index = -disp - (j_tor-2)
                pol_index = bg%nTheta - (i_pol-1)
                pol_ind_offs  = -1
                reflect = .true.
                phiref = 0.0
            else if (j_tor + disp > bg%nPhi) then
                tor_index = bg%nPhi - (j_tor + disp - bg%nPhi)
                pol_index = bg%nTheta - (i_pol-1)
                pol_ind_offs = -1
                reflect = .true.
                phiref = pi/nfp
            else
                tor_index = j_tor + disp
                pol_index = i_pol
                pol_ind_offs = 1
                reflect = .false.
                phiref = 0.0
            end if

            if (tor_index > (bg%nPhi+1) .or. tor_index <= 0) then
                stop 'bounding_planes_lateral: qhex_nBoxTor must be < bg%nPhi+1'
            end if

        else

            tor_index = ind_min_tor + n - 1
            pol_index = i_pol
            pol_ind_offs = 1
            reflect = .false.

        end if

        ! Populate array parameters for the toroidal front and back planes
        call stell_reflect(reflect, phiref, pol_index, tor_index, &
                           bg%x, bg%y, bg%z, &
                           arr_base_pb_x(n), arr_base_pb_y(n), arr_base_pb_z(n))
        call stell_reflect(reflect, phiref, pol_index+pol_ind_offs, tor_index, &
                           bg%x, bg%y, bg%z, &
                           arr_base_pf_x(n), arr_base_pf_y(n), arr_base_pf_z(n))
        call stell_reflect(reflect, phiref, pol_index, tor_index, &
                           bg%xSurf, bg%ySurf, bg%zSurf, &
                           arr_top_pb_x(n), arr_top_pb_y(n), arr_top_pb_z(n))
        call stell_reflect(reflect, phiref, pol_index+pol_ind_offs, tor_index, &
                           bg%xSurf, bg%ySurf, bg%zSurf, &
                           arr_top_pf_x(n), arr_top_pf_y(n), arr_top_pf_z(n))

    end do

    !---------------------------------------------------------------------------
    ! Part 3: Calculate the plane parameters using the plane_boxcar subroutine,
    !         given the input points determined in the previous steps.
    !---------------------------------------------------------------------------

    ! Toroidal back plane
    call plane_boxcar(n_arr_pol,                                   &
                      arr_base_tb_x, arr_base_tb_y, arr_base_tb_z, &
                      arr_top_tb_x, arr_top_tb_y, arr_top_tb_z,    &
                      qhex%ntbx, qhex%ntby, qhex%ntbz,             &
                      qhex%otbx, qhex%otby, qhex%otbz               )
    qhex%ntbx = -qhex%ntbx
    qhex%ntby = -qhex%ntby
    qhex%ntbz = -qhex%ntbz
    qhex%otbx = &
        0.5 * (bg%x(i_pol, j_tor) + bg%x(i_pol+1, j_tor)) & 
        - 0.5 * gap_qhx * qhex%ntbx
    qhex%otby = &
        0.5 * (bg%y(i_pol, j_tor) + bg%y(i_pol+1, j_tor)) &
        - 0.5 * gap_qhx * qhex%ntby
    qhex%otbz = &
        0.5 * (bg%z(i_pol, j_tor) + bg%z(i_pol+1, j_tor)) &
        - 0.5 * gap_qhx * qhex%ntbz

    ! Toroidal front plane
    call plane_boxcar(n_arr_pol,                                   & 
                      arr_base_tf_x, arr_base_tf_y, arr_base_tf_z, &
                      arr_top_tf_x, arr_top_tf_y, arr_top_tf_z,    &
                      qhex%ntfx, qhex%ntfy, qhex%ntfz,             &
                      qhex%otfx, qhex%otfy, qhex%otfz               )
    qhex%otfx = & 
        0.5 * (bg%x(i_pol, j_tor+1) + bg%x(i_pol+1, j_tor+1)) &
        - 0.5 * gap_qhx * qhex%ntfx
    qhex%otfy = &
        0.5 * (bg%y(i_pol, j_tor+1) + bg%y(i_pol+1, j_tor+1)) &
        - 0.5 * gap_qhx * qhex%ntfy
    qhex%otfz = &
        0.5 * (bg%z(i_pol, j_tor+1) + bg%z(i_pol+1, j_tor+1)) &
        - 0.5 * gap_qhx * qhex%ntfz

    ! Poloidal back plane
    call plane_boxcar(n_arr_tor,                                   &
                      arr_base_pb_x, arr_base_pb_y, arr_base_pb_z, &
                      arr_top_pb_x, arr_top_pb_y, arr_top_pb_z,    &
                      qhex%npbx, qhex%npby, qhex%npbz,             &
                      qhex%opbx, qhex%opby, qhex%opbz               )
    qhex%opbx = qhex%opbx - 0.5 * gap_qhx * qhex%npbx
    qhex%opby = qhex%opby - 0.5 * gap_qhx * qhex%npby
    qhex%opbz = qhex%opbz - 0.5 * gap_qhx * qhex%npbz
        
    ! Poloidal front plane
    call plane_boxcar(n_arr_tor,                                   & 
                      arr_base_pf_x, arr_base_pf_y, arr_base_pf_z, &
                      arr_top_pf_x, arr_top_pf_y, arr_top_pf_z,    &
                      qhex%npfx, qhex%npfy, qhex%npfz,             &
                      qhex%opfx, qhex%opfy, qhex%opfz               )
    qhex%npfx = -qhex%npfx
    qhex%npfy = -qhex%npfy
    qhex%npfz = -qhex%npfz
    qhex%opfx = qhex%opfx - 0.5 * gap_qhx * qhex%npfx
    qhex%opfy = qhex%opfy - 0.5 * gap_qhx * qhex%npfy
    qhex%opfz = qhex%opfz - 0.5 * gap_qhx * qhex%npfz

    deallocate(arr_base_tb_x, arr_base_tb_y, arr_base_tb_z, &
               arr_base_tf_x, arr_base_tf_y, arr_base_tf_z, &
               arr_top_tb_x, arr_top_tb_y, arr_top_tb_z,    &
               arr_top_tf_x, arr_top_tf_y, arr_top_tf_z,    &
               arr_base_pb_x, arr_base_pb_y, arr_base_pb_z, &
               arr_base_pf_x, arr_base_pf_y, arr_base_pf_z, &
               arr_top_pb_x, arr_top_pb_y, arr_top_pb_z,    &
               arr_top_pf_x, arr_top_pf_y, arr_top_pf_z      )

end subroutine bounding_planes_lateral

!-------------------------------------------------------------------------------
! stell_reflect(reflect, i, j, x, y, z, xr, yr, zr)
!
! Wrapper function that selects elements from input coordinate arrays, 
! performing a stellarator-symmetric reflection if desired
!-------------------------------------------------------------------------------
subroutine stell_reflect(reflect, phiref, i, j, x, y, z, xout, yout, zout)

    implicit none

    logical, intent(IN) :: reflect
    integer, intent(IN) :: i, j
    real(dp), intent(IN) :: phiref, x(:,:), y(:,:), z(:,:)
    real(dp), intent(OUT) :: xout, yout, zout
    real(dp) :: phi, r
    
    if (reflect) then

        r = sqrt(x(i,j)**2 + y(i,j)**2)
        phi = atan2(y(i,j), x(i,j))
        xout = r * cos(phiref + (phiref-phi))
        yout = r * sin(phiref + (phiref-phi))
        zout = -z(i,j)

    else

        xout = x(i,j)
        yout = y(i,j)
        zout = z(i,j)
 
    end if

end subroutine stell_reflect 
        
        
!-------------------------------------------------------------------------------
! plane_boxcar(n, base_x, base_y, base_z, top_x, top_y, top_z, 
!              nx, ny, nz, ox, oy, oz)
!
! Determines the normal vector and reference point for a trapezoid bounding
! plane by averaging preselected arrays of input points
!-------------------------------------------------------------------------------
subroutine plane_boxcar(n, base_x, base_y, base_z, top_x, top_y, top_z, &
                        nx, ny, nz, ox, oy, oz)

    use geometry, only: cross_prod, unit_vector

    implicit none

    integer, intent(IN) :: n
    real(dp), allocatable, intent(IN) :: base_x(:), base_y(:), base_z(:)
    real(dp), allocatable, intent(IN) :: top_x(:), top_y(:), top_z(:)
    real(dp), intent(OUT) :: nx, ny, nz, ox, oy, oz
    integer :: i
    real(dp) :: radial_dx, radial_dy, radial_dz, norm_r
    real(dp) :: angular_dx, angular_dy, angular_dz, norm_a, adx, ady, adz
    real(dp) :: perp_x, perp_y, perp_z, norm_perp
    real(dp), allocatable :: rdx(:), rdy(:), rdz(:)

    allocate(rdx(n), rdy(n), rdz(n))
    radial_dx  = 0.0; radial_dy  = 0.0; radial_dz  = 0.0
    angular_dx = 0.0; angular_dy = 0.0; angular_dz = 0.0
    ox         = 0.0; oy         = 0.0; oz         = 0.0

    do i = 1, n

        ! Construct the reference point as the average of all input points
        ox = ox + base_x(i) + top_x(i)
        oy = oy + base_y(i) + top_y(i)
        oz = oz + base_z(i) + top_z(i)

        ! Construct the radial vector as the avg. of top-to-base unit vectors
        rdx(i) = base_x(i) - top_x(i)
        rdy(i) = base_y(i) - top_y(i)
        rdz(i) = base_z(i) - top_z(i)

        norm_r = sqrt(rdx(i)**2 + rdy(i)**2 + rdz(i)**2)

        radial_dx = radial_dx + (rdx(i)/norm_r)
        radial_dy = radial_dy + (rdy(i)/norm_r)
        radial_dz = radial_dz + (rdz(i)/norm_r)

        if (i > 1) then

            ! Construct the angular vector as the avg. of angular unit vectors
            adx = (top_x(i) + 0.5*rdx(i)) - (top_x(i-1) + 0.5*rdx(i-1))
            ady = (top_y(i) + 0.5*rdy(i)) - (top_y(i-1) + 0.5*rdy(i-1))
            adz = (top_z(i) + 0.5*rdz(i)) - (top_z(i-1) + 0.5*rdz(i-1))

            norm_a = sqrt(adx**2 + ady**2 + adz**2)

            angular_dx = angular_dx + adx/norm_a
            angular_dy = angular_dy + ady/norm_a
            angular_dz = angular_dz + adz/norm_a

        end if

    end do

    ! Division for the average calculations
    ox = ox/(2*n)
    oy = oy/(2*n)
    oz = oz/(2*n)
    radial_dx = radial_dx/n
    radial_dy = radial_dy/n
    radial_dz = radial_dz/n
    angular_dx = angular_dx/(n-1)
    angular_dy = angular_dy/(n-1)
    angular_dz = angular_dz/(n-1)

    ! Plane-normal vectors
    call cross_prod(angular_dx, angular_dy, angular_dz, &
                    radial_dx,  radial_dy,  radial_dz, perp_x, perp_y, perp_z)
    call unit_vector(perp_x, perp_y, perp_z, nx, ny, nz)

    deallocate(rdx, rdy, rdz)

end subroutine plane_boxcar

!-------------------------------------------------------------------------------
! base_plane(i, j, bg, qhex)
!
! Determines the base plane parameters for a (lower-layer) qhex according to
! base grid and reference surface parameters
!
! Input parameters:
!     integer :: i, j        -> poloidal and toroidal index of qhex in grid
!     type(base_grid_struct) -> base grid that determines the base plane height
!
! Updated parameters:
!     type(quad_hexahedron)  -> qhex whose base is to be determined. Must 
!                               already have the rx, ry, rz, snx, sny, and snz
!                               fields determined.
!-------------------------------------------------------------------------------
subroutine base_plane(i, j, bg, qhex)

    use geometry,        only: plane_elev
    use qhex_properties, only: quad_hexahedron
    use base_grid,       only: base_grid_struct

    implicit none

    integer, intent(IN) :: i, j
    type(base_grid_struct), intent(IN) :: bg
    type(quad_hexahedron), intent(INOUT) :: qhex
    real(dp) :: h1, h2, h3, h4, h_base
  
    ! Reference point for base plane (from lowest vertex point)
    h1 = plane_elev(bg%x(i,j),  bg%y(i,j), bg%z(i,j),  &
                    qhex%rx,    qhex%ry,   qhex%rz,    &
                    -qhex%snx, -qhex%sny,  -qhex%snz    )
    h2 = plane_elev(bg%x(i,j+1), bg%y(i,j+1), bg%z(i,j+1), &
                    qhex%rx,     qhex%ry,     qhex%rz,     &
                    -qhex%snx,   -qhex%sny,   -qhex%snz     )
    h3 = plane_elev(bg%x(i+1,j), bg%y(i+1,j), bg%z(i+1,j), &
                    qhex%rx,     qhex%ry,     qhex%rz,     &
                    -qhex%snx,   -qhex%sny,   -qhex%snz     )
    h4 = plane_elev(bg%x(i+1,j+1), bg%y(i+1,j+1), bg%z(i+1,j+1), &
                    qhex%rx,       qhex%ry,       qhex%rz,       &
                    -qhex%snx,     -qhex%sny,     -qhex%snz       )
    h_base = min(h1, h2, h3, h4)
    qhex%obx = qhex%rx - qhex%snx * h_base
    qhex%oby = qhex%ry - qhex%sny * h_base
    qhex%obz = qhex%rz - qhex%snz * h_base

end subroutine base_plane

!-------------------------------------------------------------------------------
! qhex_axis(surf, theta0, phi0, qh)
!
! Defines a point on the reference surface from which the axis of a qhex should
! originate. 
!
! Input parameters:
!     type(surface) :: surf -> reference surface
!     real(dp) :: theta0, phi0 -> initial guesses for theta and phi for the 
!                                 reference point
!
! Updating parameters:
!     type(quad_hexahedron) :: qh 
!         -> qhex for which the reference point and axis are to be determined.
!            The lateral face parameters (origin points and normal vectors) must
!            be determined before calling. This subroutine assigns values to
!            the following fields:
!                ref_theta, ref_phi, rx, ry, rz, snx, sny, snz
!-------------------------------------------------------------------------------
subroutine qhex_axis(surf, theta0, phi0, qh)

    use geometry,        only: two_plane_intersect
    use surface_calc,    only: surface, surf_r, surf_z, surf_unorm
    use surface_solve,   only: surf_dist_3d
    use base_grid,       only: base_grid_struct
    use qhex_properties, only: quad_hexahedron

    implicit none

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta0, phi0
    type(quad_hexahedron), intent(INOUT) :: qh
    real(dp) :: ox1, oy1, oz1, ox2, oy2, oz2, ox3, oy3, oz3, ox4, oy4, oz4
    real(dp) :: ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3, ax4, ay4, az4
    real(dp) :: ix1, iy1, iz1, ix2, iy2, iz2, ix3, iy3, iz3, ix4, iy4, iz4
    real(dp) :: theta1, theta2, theta3, theta4, phi1, phi2, phi3, phi4
    real(dp) :: l0=0.0, sl, sx, sy, sz, chi2, r, rnx, rny, rnz

    ! Edge at corner 1: toroidal back plane and poloidal back plane
    call two_plane_intersect( &
             qh%otbx, qh%otby, qh%otbz, qh%ntbx, qh%ntby, qh%ntbz, &
             qh%opbx, qh%opby, qh%opbz, qh%npbx, qh%npby, qh%npbz, &
             ox1, oy1, oz1, ax1, ay1, az1)

    ! Edge at corner 2: toroidal front plane and poloidal back plane
    call two_plane_intersect( &
             qh%otfx, qh%otfy, qh%otfz, qh%ntfx, qh%ntfy, qh%ntfz, &
             qh%opbx, qh%opby, qh%opbz, qh%npbx, qh%npby, qh%npbz, &
             ox2, oy2, oz2, ax2, ay2, az2)

    ! Edge at corner 3: toroidal front plane and poloidal front plane
    call two_plane_intersect( &
             qh%otfx, qh%otfy, qh%otfz, qh%ntfx, qh%ntfy, qh%ntfz, &
             qh%opfx, qh%opfy, qh%opfz, qh%npfx, qh%npfy, qh%npfz, &
             ox3, oy3, oz3, ax3, ay3, az3)

    ! Edge at corner 4: toroidal back plane and poloidal front plane
    call two_plane_intersect( &
             qh%otbx, qh%otby, qh%otbz, qh%ntbx, qh%ntby, qh%ntbz, &
             qh%opfx, qh%opfy, qh%opfz, qh%npfx, qh%npfy, qh%npfz, &
             ox4, oy4, oz4, ax4, ay4, az4)

    ! Vessel intersection points of each edge
    call surf_dist_3d(surf, ox1, oy1, oz1, ax1, ay1, az1, l0, theta0, phi0, &
                      sl, theta1, phi1, sx, sy, sz, chi2)
    call surf_dist_3d(surf, ox2, oy2, oz2, ax2, ay2, az2, l0, theta0, phi0, &
                      sl, theta2, phi2, sx, sy, sz, chi2)
    call surf_dist_3d(surf, ox3, oy3, oz3, ax3, ay3, az3, l0, theta0, phi0, &
                      sl, theta3, phi3, sx, sy, sz, chi2)
    call surf_dist_3d(surf, ox4, oy4, oz4, ax4, ay4, az4, l0, theta0, phi0, &
                      sl, theta4, phi4, sx, sy, sz, chi2)

    ! Average theta and phi among the intersection points
    qh%ref_theta = 0.25*(theta1 + theta2 + theta3 + theta4)
    qh%ref_phi   = 0.25*(phi1   + phi2   + phi3   + phi4  )

    ! Coordinates of the average intersection point
    r = surf_r(surf, qh%ref_theta, qh%ref_phi)
    qh%rx = r*cos(qh%ref_phi)
    qh%ry = r*sin(qh%ref_phi)
    qh%rz = surf_z(surf, qh%ref_theta, qh%ref_phi)

    ! Unit normal vector
    call surf_unorm(surf, qh%ref_theta, qh%ref_phi, rnx, rny, rnz)
    qh%snx = -rnx
    qh%sny = -rny
    qh%snz = -rnz

end subroutine qhex_axis

!-------------------------------------------------------------------------------
! add_upper_qhex(qhex, uqhex)
!
! If space permits, adds a second qhex that abuts the top surface of 
! an existing qhex. Primarily useful for adding magnet mass in concave regions
! of the vessel.
!
! Input parameters:
!     type(surface)         :: surf  -> Vessel surface or other reference
!                                       surface that limits the qhex size
!     type(quad_hexahedron) :: qhex  -> The qhex on top of which a new one
!                                       is to be added
!
! Output parameters:
!     type(quad_hexahedron) :: uqhex -> The new qhex.
!-------------------------------------------------------------------------------
subroutine add_upper_qhex(surf, qhex, max_height, uqhex)

    use magpie_globals,  only: gap_rad, gap_lim
    use geometry,        only: cross_prod, unit_vector
    use surface_calc,    only: surface
    use surface_solve,   only: surf_dist_3d
    use qhex_properties, only: quad_hexahedron, qhex_corners, qhex_volume_ctr, &
                               check_qhex
    use qhex_adjust,     only: reset_qhex_height

    implicit none

    type(surface), intent(IN) :: surf
    type(quad_hexahedron), intent(IN)  :: qhex
    real(dp), intent(IN) :: max_height
    type(quad_hexahedron), intent(OUT) :: uqhex
    real(dp) :: xvec_tb, yvec_tb, zvec_tb, xvec_tf, yvec_tf, zvec_tf
    real(dp) :: xvec_pb, yvec_pb, zvec_pb, xvec_pf, yvec_pf, zvec_pf
    real(dp) :: targ_dist, targ_x, targ_y, targ_z
    real(dp), dimension(4) :: cvd, cx, cy, cz, cvix, cviy, cviz
    integer :: i, i_opt
    real(dp) :: zero = 0.0, theta0, phi0, theta, phi, chi2 
    real(dp) :: cxmean, cymean, czmean

    ! Properties carried over from the input qhexezoid
    uqhex%vx = qhex%vx;    uqhex%vy = qhex%vy;    uqhex%vz = qhex%vz
    uqhex%rx = qhex%rx;    uqhex%ry = qhex%ry;    uqhex%rz = qhex%rz
    uqhex%snx = qhex%snx;  uqhex%sny = qhex%sny;  uqhex%snz = qhex%snz
    uqhex%obx = qhex%otx + gap_rad*uqhex%snx
    uqhex%oby = qhex%oty + gap_rad*uqhex%sny
    uqhex%obz = qhex%otz + gap_rad*uqhex%snz
    uqhex%otx = uqhex%vx - gap_lim*uqhex%snx
    uqhex%oty = uqhex%vy - gap_lim*uqhex%sny
    uqhex%otz = uqhex%vz - gap_lim*uqhex%snz
    uqhex%ves_theta = qhex%ves_theta
    uqhex%ves_phi = qhex%ves_phi
    uqhex%ref_theta = qhex%ref_theta
    uqhex%ref_phi = uqhex%ref_phi

    cx = (/ qhex%ctx1, qhex%ctx2, qhex%ctx3, qhex%ctx4 /)
    cy = (/ qhex%cty1, qhex%cty2, qhex%cty3, qhex%cty4 /)
    cz = (/ qhex%ctz1, qhex%ctz2, qhex%ctz3, qhex%ctz4 /)
    cxmean = 0.0; cymean = 0.0; czmean = 0.0
    do i = 1, 4
        call surf_dist_3d(surf, cx(i), cy(i), cz(i), qhex%snx, qhex%sny, &
                          qhex%snz, zero, qhex%ves_theta, qhex%ves_phi,      &
                          cvd(i), theta, phi, cvix(i), cviy(i), cviz(i), chi2)
        cxmean = cxmean + cx(i)*cvd(i)
        cymean = cymean + cy(i)*cvd(i)
        czmean = czmean + cz(i)*cvd(i)
    end do
    i_opt = maxloc(cvd, 1)
    targ_dist = 2*cvd(i_opt)
    cxmean = cxmean / sum(cvd)
    cymean = cymean / sum(cvd)
    czmean = czmean / sum(cvd)
    targ_x = cxmean + targ_dist * qhex%snx
    targ_y = cymean + targ_dist * qhex%sny
    targ_z = czmean + targ_dist * qhex%snz

    !---------------------------------------------------------------------------
    ! Normal vectors and reference points of the lateral faces
    !---------------------------------------------------------------------------

    ! Toroidal back face
    call cross_prod( &
             qhex%ctx4-qhex%ctx1, qhex%cty4-qhex%cty1, qhex%ctz4-qhex%ctz1, &
             targ_x-qhex%ctx4,  targ_y-qhex%cty4,  targ_z-qhex%ctz4,  &
             xvec_tb, yvec_tb, zvec_tb)
    call unit_vector(xvec_tb,    yvec_tb,    zvec_tb, &
                     uqhex%ntbx, uqhex%ntby, uqhex%ntbz)
    uqhex%otbx = qhex%ctx4
    uqhex%otby = qhex%cty4
    uqhex%otbz = qhex%ctz4

    ! Toroidal front face
    call cross_prod( &
             qhex%ctx2-qhex%ctx3, qhex%cty2-qhex%cty3, qhex%ctz2-qhex%ctz3, &
             targ_x-qhex%ctx2,  targ_y-qhex%cty2,  targ_z-qhex%ctz2,  &
             xvec_tf, yvec_tf, zvec_tf)
    call unit_vector(xvec_tf,    yvec_tf,    zvec_tf, &
                     uqhex%ntfx, uqhex%ntfy, uqhex%ntfz)
    uqhex%otfx = qhex%ctx2
    uqhex%otfy = qhex%cty2
    uqhex%otfz = qhex%ctz2

    ! Poloidal back face
    call cross_prod( &
             qhex%ctx1-qhex%ctx2, qhex%cty1-qhex%cty2, qhex%ctz1-qhex%ctz2, &
             targ_x-qhex%ctx1,  targ_y-qhex%cty1,  targ_z-qhex%ctz1,  &
             xvec_pb, yvec_pb, zvec_pb)
    call unit_vector(xvec_pb,    yvec_pb,    zvec_pb, &
                     uqhex%npbx, uqhex%npby, uqhex%npbz)
    uqhex%opbx = qhex%ctx1
    uqhex%opby = qhex%cty1
    uqhex%opbz = qhex%ctz1

    ! Poloidal front face
    call cross_prod( &
             qhex%ctx3-qhex%ctx4, qhex%cty3-qhex%cty4, qhex%ctz3-qhex%ctz4, &
             targ_x-qhex%ctx3,  targ_y-qhex%cty3,  targ_z-qhex%ctz3,  &
             xvec_pf, yvec_pf, zvec_pf)
    call unit_vector(xvec_pf,    yvec_pf,    zvec_pf, &
                     uqhex%npfx, uqhex%npfy, uqhex%npfz)
    uqhex%opfx = qhex%ctx3
    uqhex%opfy = qhex%cty3
    uqhex%opfz = qhex%ctz3

    call qhex_corners(uqhex)
    call reset_qhex_height(surf, uqhex, max_height)
    call qhex_volume_ctr(uqhex)
    call check_qhex(uqhex)

end subroutine add_upper_qhex

!-------------------------------------------------------------------------------
! remove_overlaps(qhArr)
!
! Checks each adjacent pair of qhex-pair structures in a toroidal grid for 
! overlaps and shrinks them until they no longer overlap or one becomes 
! erroneous.
!
! Presently does not check for overlaps between the edge qhex pairs on adjacent
! half-periods.
!
! Updating parameter:
!     type(qhex_array) :: qhArr  -> the qhex array to be cleared of overlaps.
!
! Input parameter:
!     real(dp) :: adj_interval   -> the adjustment interval for retracting qhex
!                                   faces when attempting to remove overlaps
!     real(dp) :: max_ht_up      -> Maximum allowable height for upper qhex if
!                                   applicable
!     type(surface) :: surf      -> Limiting surface for upper qhex
!-------------------------------------------------------------------------------
subroutine remove_overlaps(qhArr, adj_interval, max_ht_up, surf)

    use magpie_globals,  only: nfp, pi, stell_symm, tor_symm
    use surface_calc,    only: surface
    use qhex_properties, only: qhex_pair, qhex_array
    use qhex_adjust,     only: stell_transf_mode_len, stell_qhPair_transform

    implicit none

    type(qhex_array), intent(INOUT) :: qhArr
    real(dp), intent(IN) :: adj_interval, max_ht_up
    type(surface), intent(IN) :: surf
    integer :: i, j, nGridPoints, i_pol, i_tor, j_pol, j_tor
    character(len=stell_transf_mode_len) :: mode 
    type(qhex_pair) :: transf_pair_phi0, transf_pair_phi1
    integer :: ind0, ind1, inv_sign
    real(dp) :: phi0, phi1

    nGridPoints = qhArr%nTheta * qhArr%nPhi

    ! Perform the operation for all unique pairs
    do i = 1, nGridPoints
        do j = i+1, nGridPoints
            
            ! Toroidal and poloidal indices for the i- and j- grid points
            i_pol = mod(i-1, qhArr%nTheta) + 1
            i_tor = (i-1)/qhArr%nTheta + 1
            j_pol = mod(j-1, qhArr%nTheta) + 1
            j_tor = (j-1)/qhArr%nTheta + 1

            ! Check for overlaps between each unique pair of qhex-pairs
            call adjust_pair_via_lower(qhArr%pairs(i_pol, i_tor), &
                                       qhArr%pairs(j_pol, j_tor), &
                                       adj_interval, max_ht_up, surf)
        end do
    end do

    ! If symmetry is assumed, transform toroidal edge pairs and fix overlaps
    if (stell_symm .or. tor_symm) then

        ! Indices for the toroidal start and end of the magnet array
        ind0 = 1
        ind1 = qhArr%nPhi

        ! Set input values that depend on the type of assumed symmetry
        if (stell_symm) then
            mode = 'reflect'
            phi0 = 0          ! Reflection plane, toroidal start of magnet array
            phi1 = pi/nfp     ! Reflection plane, toroidal end of magnet array
            inv_sign = 1      ! Reflection plane remains the same for inversion
        else if (tor_symm) then
            mode = 'translate'
            phi0 = 2*pi/nfp   ! Toroidal interval, toroidal start of magnet arr
            phi1 = phi0       ! Toroidal interval, toroidal end of magnet arr
            inv_sign = -1     ! Toroidal interval must be negated for inversion
        end if
            
        do i = 1, qhArr%nTheta

            ! Create transformed version of an original pair
            call stell_qhPair_transform( &
                     mode, phi0, qhArr%pairs(i,ind0), transf_pair_phi0)
            call stell_qhPair_transform( &
                     mode, phi1, qhArr%pairs(i,ind1), transf_pair_phi1)

            ! Adjust the transformed pair for overlaps with all original pairs
            do j = 1, qhArr%nTheta
                call adjust_pair_via_lower(qhArr%pairs(j,ind0), &
                         transf_pair_phi0, adj_interval, max_ht_up, surf)
                call adjust_pair_via_lower(qhArr%pairs(j,ind1), &
                         transf_pair_phi1, adj_interval, max_ht_up, surf)
            end do

            ! Replace original pair with inverse-transform of transformed pair
            call stell_qhPair_transform( &
                     mode, inv_sign*phi0, transf_pair_phi0, qhArr%pairs(i,ind0))
            call stell_qhPair_transform( &
                     mode, inv_sign*phi1, transf_pair_phi1, qhArr%pairs(i,ind1))
        end do

    end if 

    ! Update the tally of erroneous magnets
    qhArr%nErrLo = 0
    qhArr%nErrUp = 0
    do i = 1, qhArr%nTheta
        do j = 1, qhArr%nPhi

            if (qhArr%pairs(i,j)%lo%err) &
                qhArr%nErrLo = qhArr%nErrLo + 1
            if (qhArr%incl_upper .and. qhArr%pairs(i,j)%up%err) &
                qhArr%nErrUp = qhArr%nErrUp + 1

        end do
    end do
    qhArr%nErr = qhArr%nErrLo + qhArr%nErrUp

end subroutine remove_overlaps

!-------------------------------------------------------------------------------
! adjust_pair_via_lower(pair1, pair2, interval, max_ht_up, surf)
!
! Removes overlaps between qhex pairs by retracting lateral faces of lower
! quex and then replacing the upper qhex with a call to add_upper_qhex.
!
! Input parameters:
!     real(dp) :: interval   -> Distance by which lateral faces of the lower
!                               qhex are to be retracted each time an overlap
!                               is found
!     real(dp) :: max_ht_up  -> Maximum height for the upper qhexes
!     type(surface) :: surf  -> Limiting surface for upper qhex
! 
! Updating parameters
!     type(qhex_pair) :: pair1, pair2 -> qhex pairs for adjusting
!-------------------------------------------------------------------------------
subroutine adjust_pair_via_lower(pair1, pair2, interval, max_ht_up, surf)

    use surface_calc,    only: surface
    use qhex_properties, only: qhex_pair, nHexFaces
    use qhex_overlap,    only: qhex_overlapping
    use qhex_adjust,     only: adjust_lateral_faces, check_and_adjust

    implicit none

    type(qhex_pair), intent(INOUT) :: pair1, pair2
    real(dp), intent(IN) :: interval, max_ht_up
    type(surface), intent(IN) :: surf
    logical, dimension(nHexFaces) :: ovl1, ovl2
    logical :: lowers_adjusted

    lowers_adjusted = .false.

    ! Resolve overlaps between the lower qhex of each pair
    do while(qhex_overlapping(pair1%lo, pair2%lo, ovl1, ovl2))
        if (pair1%lo%err .or. pair2%lo%err) exit
        call adjust_lateral_faces(pair1%lo, pair2%lo, ovl1, ovl2, interval)
        lowers_adjusted = .true.
    end do

    ! Resolve overlaps between pair1-upper and pair2-lower if applicable
    if (pair1%incl_upper) then
        if (lowers_adjusted) &
            call add_upper_qhex(surf, pair1%lo, max_ht_up, pair1%up)
        do while(qhex_overlapping(pair1%up, pair2%lo, ovl1, ovl2))
            if (pair1%up%err .or. pair2%lo%err) exit
            call adjust_lateral_faces(pair1%lo, pair2%lo, ovl1, ovl2, interval)
            call add_upper_qhex(surf, pair1%lo, max_ht_up, pair1%up)
        end do
    end if

    ! Resolve overlaps between pair1-lower and pair2-upper if applicable
    if (pair2%incl_upper) then
        if (lowers_adjusted) &
            call add_upper_qhex(surf, pair2%lo, max_ht_up, pair2%up)
        do while(qhex_overlapping(pair1%lo, pair2%up, ovl1, ovl2))
            if (pair1%lo%err .or. pair2%up%err) exit
            call adjust_lateral_faces(pair1%lo, pair2%lo, ovl1, ovl2, interval)
            call add_upper_qhex(surf, pair2%lo, max_ht_up, pair2%up)
        end do
    end if

    ! Resolve overlaps between the upper quexes of each pair if applicable
    if (pair1%incl_upper .and. pair2%incl_upper) then
        do while(qhex_overlapping(pair1%up, pair2%up, ovl1, ovl2))
            if (pair1%up%err .or. pair2%up%err) exit
            call adjust_lateral_faces(pair1%lo, pair2%lo, ovl1, ovl2, interval)
            call add_upper_qhex(surf, pair1%lo, max_ht_up, pair1%up)
            call add_upper_qhex(surf, pair2%lo, max_ht_up, pair2%up)
        end do
    end if

end subroutine adjust_pair_via_lower

end module qhex_array


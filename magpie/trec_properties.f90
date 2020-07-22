module trec_properties

use magpie_globals, only: dp

implicit none

type trec
    real(dp) :: xib, yib, rib        ! x, y, r coords, inner back  corners
    real(dp) :: xif, yif, rif        ! x, y, r coords, inner front corners
    real(dp) :: xob, yob, rob        ! x, y, r coords, outer back  corners
    real(dp) :: xof, yof, rof        ! x, y, r coords, outer front corners
    real(dp) :: zl,  zu              ! z coords, lower and upper corners
    real(dp) :: cmx, cmy, cmz, cmr   ! x, y, z coords, centroid
    real(dp) :: axx, axy, axz        ! x, y, z components of polarization axis
    real(dp) :: rhat_x, rhat_y       ! x, y components, radial unit vector
                                     ! (normal to inboard/outboard faces)
    real(dp) :: phat_x, phat_y       ! x, y components, toroidal unit vector
                                     ! (normal to front/back faces)
    real(dp) :: lr, lp, lz           ! dimensions (meters) in r, "phi", and z
    real(dp) :: vol                  ! volume (meters^3)
end type trec

! Indices labeling the toroidal angles for the vertices and the centroid
integer, parameter :: iIB = 1, iOB = 2, iCM = 3, iOF = 4, iIF = 5

integer, parameter :: nTrecFaces = 6

contains

!-------------------------------------------------------------------------------
! new_trec(cmr, cmp, cmz, lr, lp, lz)
!
! Constructs a trapezoidally-enclosed rectangular prism (trec) structure.
!
! Input parameters
!     real(dp) :: cmr, cmp, cmz    -> r, phi, z coordinates of the centroid
!     real(dp) :: lr, lp, lz       -> lengths of the radial, "toroidal", and 
!                                     vertical dimensions
!-------------------------------------------------------------------------------
function new_trec(cmr, cmp, cmz, lr, lp, lz) result(tr)

    implicit none

    real(dp), intent(IN) :: cmr, cmp, cmz, lr, lp, lz
    real(dp) :: rhat_x, rhat_y, phat_x, phat_y
    type(trec) :: tr

    ! Centroid coordinates
    tr%cmr = cmr
    tr%cmx = cmr*cos(cmp)
    tr%cmy = cmr*sin(cmp)
    tr%cmz = cmz

    ! Dimensions
    tr%lr = lr
    tr%lp = lp
    tr%lz = lz

    tr%vol = lr*lp*lz

    ! Unit vectors in the radial and toroidal directions (at the centroid)
    rhat_x =  cos(cmp)
    rhat_y =  sin(cmp)
    phat_x = -sin(cmp)
    phat_y =  cos(cmp)

    ! Inner back corner (back = side with lower toroidal angle)
    tr%xib = tr%cmx - 0.5 * lr * rhat_x - 0.5 * lp * phat_x
    tr%yib = tr%cmy - 0.5 * lr * rhat_y - 0.5 * lp * phat_y
    tr%rib = sqrt(tr%xib**2 + tr%yib**2)

    ! Inner front corner (front = side with greater toroidal angle)
    tr%xif = tr%cmx - 0.5 * lr * rhat_x + 0.5 * lp * phat_x
    tr%yif = tr%cmy - 0.5 * lr * rhat_y + 0.5 * lp * phat_y
    tr%rif = sqrt(tr%xif**2 + tr%yif**2)

    ! Outer back corner (back = side with lower toroidal angle)
    tr%xob = tr%cmx + 0.5 * lr * rhat_x - 0.5 * lp * phat_x
    tr%yob = tr%cmy + 0.5 * lr * rhat_y - 0.5 * lp * phat_y
    tr%rob = sqrt(tr%xob**2 + tr%yob**2)

    ! Outer front corner (front = side with greater toroidal angle)
    tr%xof = tr%cmx + 0.5 * lr * rhat_x + 0.5 * lp * phat_x
    tr%yof = tr%cmy + 0.5 * lr * rhat_y + 0.5 * lp * phat_y
    tr%rof = sqrt(tr%xof**2 + tr%yof**2)

    ! Z coordinates of lower and upper faces
    tr%zl = cmz - 0.5 * lz
    tr%zu = cmz + 0.5 * lz

    ! Initial axis (radial direction)
    tr%axx = rhat_x
    tr%axy = rhat_y
    tr%axz = 0

    ! Radial and toroidal unit vectors
    tr%rhat_x = rhat_x
    tr%rhat_y = rhat_y
    tr%phat_x = phat_x
    tr%phat_y = phat_y

end function new_trec

!-------------------------------------------------------------------------------
! trec_grid_params(dr, dz, dphi, gap, lpMin, rMin, rMax, zMin, zMax, 
!                  rTri, r0, z0, nr, nz)
!
! Determines the radial and vertical grid dimensions and starting points for
! an array of trapezoidally-enclosed rectangular prisms.
!
! If the upper and lower vertical limits zMin and zMax are on opposite sides of
! the z=0 plane, the cell boundaries will be placed such that one boundary is
! guaranteed to coincide with the z=0 plane.
!
! Input parameters:
!     real(dp) :: dr, dz     -> Radial and vertical dimensions of grid cells
!     real(dp) :: dphi       -> Subtended angle (radians) of each trapezoid 
!                               (i.e., toroidal sector)
!     real(dp) :: gap        -> Minimum clearance between prism edges and the
!                               sector boundary
!     real(dp) :: lpMin      -> Minimum allowable toroidal dimension of prisms
!     real(dp) :: rMin, rMax -> Min, max allowable radial grid point
!     real(dp) :: zMin, zMax -> Min, max allowable vertical grid point
!
! Output parameters:
!     real(dp) :: rTri       -> Minimum radial coordinate imposed by gap, dphi
!                               (may supersede rMin)
!     real(dp) :: r0, z0     -> Radial, vertical coordinates of first cells
!     integer :: nr, nz      -> Radial, vertical dimensions of grid
!-------------------------------------------------------------------------------
subroutine trec_grid_params(dr, dz, dphi, gap, lpMin, rMin, rMax, zMin, zMax, &
                            rTri, r0, z0, nr, nz)

    implicit none

    real(dp), intent(IN)  :: dr, dz, dphi, gap, lpMin, rMin, rMax, zMin, zMax
    real(dp), intent(OUT) :: rTri, r0, z0
    integer, intent(OUT)  :: nr, nz
    real(dp) :: r_for_lpMin
    integer :: nzLower, nzUpper

    ! Radial coordinate at which inner walls of trapezoids meet
    rTri = gap / sin(0.5*dphi)

    ! Minimum radial coordinate that permits prisms with minimum toroidal length
    r_for_lpMin = rTri + 0.5*lpMin / cos(0.5*dphi)

    ! Minimum/starting r coordinate of prisms in grid
    r0 = max(rMin, r_for_lpMin)
    nr = floor((rMax - r0)/dr)
 
    ! Ensure that the z=0 plane coincides with a grid cell boundary
    if (zMin < 0 .and. zMax > 0) then
        nzLower = floor(-zMin/dz)
        nzUpper = floor(zMax/dz)
        nz = nzLower + nzUpper
        z0 = -real(nzLower,dp)*dz
    else
        nz = floor((zMax-zMin)/dz)
        z0 = zMin        
    end if

end subroutine trec_grid_params

!-------------------------------------------------------------------------------
! trec_in_polygon_5plane(tr, np, rp, zp)
!
! Wrapper for in_polygon() designed to determine whether a rectangular prism
! lies within a toroidal surface by checking 12 reference points: the eight
! vertices as well as four points where the edges intersect the plane of the
! centroid whose normal vector is the toroidal unit vector.
!
! "5plane" refers to the fact that the 12 query points have one of 5 distinct
! toroidal angles; hence, 5 corresponding poloidal cross-sections of the 
! toroidal surface must be supplied to check the query points.
!
! Input parameters:
!     type(trec) :: tr  -> the prism whose reference points are to be checked
!     integer :: np     -> the number of vertices in each cross-section
!                             of the toroidal surface to be provided as input
!     real(dp), dimension(np,5) :: rp, zp
!                       -> r, z coordinates of the five cross-sections of the
!                          toroidal surface. Each column corresponds to a 
!                          different toroidal angle. The order of the angles
!                          is as determined by the index parameters supplied
!                          in this module:
!                              iIB: radially inner vertices, toroidal back side
!                              iOB: radially outer vertices, toroidal front side
!                              iCM: plane of the centroid of the prism
!                              iOF: radially outer vertices, toroidal front side
!                              iIF: radially inner vertices, toroidal back side 
!
! Returns:
!     logical :: tip    -> True if any of the reference points are found to
!                          be inside the cross-section for their respective
!                          toroidal angle.
!-------------------------------------------------------------------------------
function trec_in_polygon_5plane(tr, np, rp, zp) result(tip)

    use geometry, only: in_polygon

    implicit none

    type(trec), intent(IN) :: tr
    integer, intent(IN) :: np
    real(dp), dimension(np,5), intent(IN) :: rp, zp
    logical :: tip
    integer, parameter :: nq = 12  ! total number of query points
    real(dp), dimension(nq) :: rq, zq
    integer(dp), dimension(nq) :: ind_phiq
    real(dp) :: ric, roc
    integer :: i

    ! Radial coordinates of edges in plane of the centroid (with toroidal normal
    ! vector)
    ric = tr%cmr - 0.5*tr%lr
    roc = tr%cmr + 0.5*tr%lr

    ! Arrays of query points and indices for the corresponding phi angles
    rq = (/ tr%rib, tr%rib, tr%rob, tr%rob,   ric,    ric,    roc, &
            roc,   tr%rif, tr%rif, tr%rof, tr%rof /)
    zq = (/ tr%zl,  tr%zu,  tr%zl,  tr%zu,  tr%zl,  tr%zu,  tr%zl,  &
            tr%zu,  tr%zl,  tr%zu,  tr%zl,  tr%zu  /)
    ind_phiq = (/ iIB, iIB, iOB, iOB, iCM, iCM, iCM, iCM, iIF, iIF, iOF, iOF /)

    do i = 1, nq
        if (in_polygon(np, rp(:,ind_phiq(i)), zp(:,ind_phiq(i)), &
                       rq(i), zq(i))                              ) then
            tip = .true.
            return
        end if
    end do

    tip = .false.

end function trec_in_polygon_5plane

!-------------------------------------------------------------------------------
! trec_from_polygon_5plane(tr, np, rp, zp)
! 
! Estimates the minimum poloidal distance of a rectangular prism from a
! toroidal surface using 12 query points at a total of 5 unique toroidal angles 
! (hence the name "5plane").
!
! Input parameters:
!     type(trec) :: tr  -> the prism whose reference points are to be checked
!     integer :: np     -> the number of vertices in each cross-section
!                             of the toroidal surface to be provided as input
!     real(dp), dimension(np,5) :: rp, zp
!                       -> r, z coordinates of the five cross-sections of the
!                          toroidal surface. See documentation for 
!                          trec_in_polygon_5plane for more details.
!
! Returns:
!     real(dp) :: minDist -> Minimum distance determined for any of the 
!                            query points to a cross-section of the toroidal
!                            surface in its respective poloidal plane.
!-------------------------------------------------------------------------------
function trec_from_polygon_5plane(tr, np, rp, zp) result(minDist)

    implicit none

    type(trec), intent(IN) :: tr
    integer, intent(IN) :: np
    real(dp), dimension(np,5), intent(IN) :: rp, zp
    integer, parameter :: nq = 12
    real(dp) :: ric, roc, minDist, minDistSq, distSq_i, distSq_j
    real(dp), dimension(nq) :: rq, zq
    integer,  dimension(nq) :: ind_phiq
    integer :: i, j

    ! Radial coordinates of edges in plane of the centroid (with toroidal normal
    ! vector)
    ric = tr%cmr - 0.5*tr%lr
    roc = tr%cmr + 0.5*tr%lr

    ! Arrays of query points and indices for the corresponding phi angles
    rq = (/ tr%rib, tr%rib, tr%rob, tr%rob,   ric,    ric,    roc, &
            roc,   tr%rif, tr%rif, tr%rof, tr%rof /)
    zq = (/ tr%zl,  tr%zu,  tr%zl,  tr%zu,  tr%zl,  tr%zu,  tr%zl,  &
            tr%zu,  tr%zl,  tr%zu,  tr%zl,  tr%zu  /)
    ind_phiq = (/ iIB, iIB, iOB, iOB, iCM, iCM, iCM, iCM, iIF, iIF, iOF, iOF /)

    do i = 1, nq
        do j = 1, np
            distSq_j = (rq(i) - rp(j,ind_phiq(i)))**2 + &
                       (zq(i) - zp(j,ind_phiq(i)))**2 
            if (j == 1 .or. distSq_j < distSq_i) distSq_i = distSq_j
        end do
        if (i == 1 .or. distSq_i < minDistSq) minDistSq = distSq_i
    end do

    minDist = sqrt(minDistSq)

end function trec_from_polygon_5plane

!-------------------------------------------------------------------------------
! estimate_trec_theta(tr, np, thetap, rp, zp)
!
! Estimates the poloidal angle (theta) value for a rectangular prism by 
! setting it equal to that of the closest point on a poloidal cross-section
! of a toroidal surface at the corresponding toroidal angle. 
!
! The main intended use case is to come up with an initial guess of theta 
! for the perpendicular intersection solver (surf_perp_intersect_3d).
!
! Input parameters:
!     type(trec) :: tr    -> the prism whose theta value is to be estimated
!     integer :: np       -> number of vertices supplied in a cross-section
!                            of the toroidal surface
!     real(dp) :: thetap  -> theta values of each vertex in the cross-section
!     real(dp) :: rp, zp  -> r, z coords of each vertex in the cross-section
!
! Returns:
!     real(dp) :: theta   -> Theta value of the point in the cross-section 
!                            closest to the centroid of the prism
!-------------------------------------------------------------------------------
function estimate_trec_theta(tr, np, thetap, rp, zp) result(theta)

    implicit none

    type(trec), intent(IN) :: tr
    integer, intent(IN) :: np
    real(dp), dimension(np) :: thetap, rp, zp
    real(dp) :: theta
    integer :: i, minDistInd
    real(dp) :: distSq_i, minDistSq

    
    do i = 1, np
        distSq_i = (rp(i) - tr%cmr)**2 + (zp(i) - tr%cmz)**2
        if (i == 1 .or. distSq_i < minDistSq) then
            minDistSq = distSq_i
            minDistInd = i
        end if
    end do

    theta = thetap(minDistInd)

end function estimate_trec_theta

!-------------------------------------------------------------------------------
! trec_set_perpendicular(tr, surf, l0, theta0, phi0)
!
! Sets the axis of a trec structure to be locally perpendicular to a toroidal
! surface. 
!
! This is essentially a wrapper for surf_perp_intersect_3d.
!
! Input parameters
!     type(surface) :: surf     -> The toroidal surface
!     real(dp) :: l0            -> Initial guess for the distance of centroid 
!                                  from where a line along its axis would 
!                                  intersect the toroidal surface
!     real(dp) :: theta0, phi0  -> Initial guesses for the poloidal, toroidal
!                                  angles where a line along the axis would
!                                  intersect the toroidal surface
!
! Updating parameters
!     type(trec) :: tr          -> The rectangular prism whose axis is to be
!                                  set perpendicular
!-------------------------------------------------------------------------------
subroutine trec_set_perpendicular(tr, surf, l0, theta0, phi0)

    use surface_calc,  only: surface
    use surface_solve, only: surf_perp_intersect_3d

    implicit none

    type(trec), intent(INOUT) :: tr
    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: l0, theta0, phi0
    real(dp) :: l, theta, phi, sx, sy, sz, ux, uy, uz, chi2

    call surf_perp_intersect_3d(surf, tr%cmx, tr%cmy, tr%cmz, &
                                l0, theta0, phi0, l, theta, phi, &
                                sx, sy, sz, ux, uy, uz, chi2)
    tr%axx = -ux
    tr%axy = -uy
    tr%axz = -uz

end subroutine trec_set_perpendicular

!-------------------------------------------------------------------------------
! trec_obj_overlap(tr, maxCheckGap, ovl)
!
! Checks whether the faces of a rectangular prism intersect with any applicable 
! object (port, coil, etc.) as specified by the user.
!
! NOTE: for correct functionality, calling procedures must call the
! intersections_to_check subroutine once prior to calls to trec_obj_overlap.
!
! Input parameters
!     type(trec) :: tr        -> trec whose faces are to be checked
!     real(dp) :: maxCheckGap -> Max spacing between adjacent query points on 
!                                each face for overlap checking
!
! Output parameters:
!     logical :: ovl          -> True if any faces are found to intersect an
!                                applicable object
!-------------------------------------------------------------------------------
subroutine trec_obj_overlap(tr, maxCheckGap, ovl)

    use intersections, only: check_intersections

    implicit none

    type(trec), intent(IN) :: tr
    real(dp), intent(IN) :: maxCheckGap
    logical, intent(OUT) :: ovl
    real(dp), dimension(nTrecFaces) :: ax, ay, az, bx, by, bz, la, lb
    real(dp), dimension(nTrecFaces) :: ox, oy, oz
    real(dp) :: delta_a, delta_b, qx, qy, qz
    integer :: i, j, k, na, nb

    ! Unit vector for the first ("a") dimension along each of the six faces
    ax = (/ tr%phat_x, tr%phat_x, tr%rhat_x, tr%rhat_x, tr%rhat_x, tr%rhat_x /)
    ay = (/ tr%phat_y, tr%phat_y, tr%rhat_y, tr%rhat_y, tr%rhat_y, tr%rhat_y /)
    az = (/     0._dp,     0._dp,     0._dp,     0._dp,     0._dp,     0._dp /)

    ! Unit vector second ("b") dimension along each face
    bx = (/     0._dp,     0._dp,     0._dp,     0._dp, tr%phat_x, tr%phat_x /)
    by = (/     0._dp,     0._dp,     0._dp,     0._dp, tr%phat_y, tr%phat_y /)
    bz = (/     1._dp,     1._dp,     1._dp,     1._dp,     0._dp,     0._dp /)

    ! Dimensions of the "a" and "b" dimensions of each face
    la = (/     tr%lp,     tr%lp,     tr%lr,     tr%lr,     tr%lr,     tr%lr /)
    lb = (/     tr%lz,     tr%lz,     tr%lz,     tr%lz,     tr%lp,     tr%lp /)

    ! Origin point for each face
    ox = (/    tr%xib,    tr%xob,    tr%xib,    tr%xif,    tr%xib,    tr%xib /)
    oy = (/    tr%yib,    tr%yob,    tr%yib,    tr%yif,    tr%yib,    tr%yib /)
    oz = (/     tr%zl,     tr%zl,     tr%zl,     tr%zl,     tr%zl,     tr%zu /)

    ! Check for overlap on a grid of query points on each face
    do i = 1, nTrecFaces
        na = ceiling(la(i)/maxCheckGap) + 1
        nb = ceiling(lb(i)/maxCheckGap) + 1
        delta_a = la(i) / (real(na, dp) - 1.0)
        delta_b = lb(i) / (real(nb, dp) - 1.0)
        do j = 1, na
            do k = 1, nb
                qx = ox(i) + (j-1)*ax(i)*delta_a + (k-1)*bx(i)*delta_b
                qy = oy(i) + (j-1)*ay(i)*delta_a + (k-1)*by(i)*delta_b
                qz = oz(i) + (j-1)*az(i)*delta_a + (k-1)*bz(i)*delta_b
                call check_intersections(qx, qy, qz, ovl)
                if (ovl) return
            end do
        end do
    end do

    ovl = .false.

end subroutine trec_obj_overlap

end module trec_properties


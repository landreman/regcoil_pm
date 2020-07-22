module cbrick_properties

use magpie_globals, only: dp

implicit none

type cbrick
    ! r = radial coordinate
    ! z = vertical coordinate
    ! i = inboard face
    ! o = outboard face
    ! b = back face (toroidally)
    ! f = front face
    ! l = lower face
    ! u = upper face
    ! ax[x,y,z] = radial unit vector
    ! cm[x,y,z] = center of mass (uniform density presumed)
    real(dp) :: ri, ro
    real(dp) :: zl, zu
    real(dp) :: phib, phif, phi
    integer :: indb, indf
    real(dp) :: vol
    real(dp) :: cmx, cmy, cmz, cmr
    real(dp) :: axx, axy, axz
end type cbrick

contains

!-------------------------------------------------------------------------------
! new_cbrick(ri, ro, zl, zu, phib, phif, indb, indf)
!
! Constructs a cbrick structure. The axis vector is initialized as the radial
! unit vector at the center of mass.
!
! Input parameters
!     real(dp) :: ri, ro       -> inner, outer radius
!     real(dp) :: zl, zu       -> lower, upper z coordinate of horizontal faces
!     real(dp) :: phib, phif   -> back, front toroidal angle
!-------------------------------------------------------------------------------
function new_cbrick(ri, ro, zl, zu, phib, phif) result(cb)

    implicit none

    real(dp), intent(IN) :: ro, ri, zl, zu, phib, phif
    type(cbrick) :: cb

    ! Directly input quantities
    cb%ri = ri
    cb%ro = ro
    cb%zl = zl
    cb%zu = zu
    cb%phib = phib
    cb%phif = phif

    ! Derived quantities
    cb%phi = 0.5*(phib + phif)
    cb%vol = 0.5*(phif - phib)*(ro**2 - ri**2)*(zu - zl)
    cb%cmx = (1./(3.*cb%vol))*(ro**3 - ri**3)*(sin(phif) - sin(phib))*(zu - zl)
    cb%cmy = (1./(3.*cb%vol))*(ro**3 - ri**3)*(cos(phib) - cos(phif))*(zu - zl)
    cb%cmz = (1./(4.*cb%vol))*(ro**2 - ri**2)*(phif - phib)*(zu**2 - zl**2)
    cb%cmr = sqrt(cb%cmx**2 + cb%cmy**2)
    cb%axx = cos(cb%phi)
    cb%axy = sin(cb%phi)
    cb%axz = 0

end function new_cbrick

!-------------------------------------------------------------------------------
! corner_in_polygon(cb, np, rp, zp)
!
! Determines whether the corners of an annular sector lie within a polygon
! in the r-z plane.
!
! Input parameters
!     type(cbrick) :: cb    -> curved brick whose corners are to be checked
!     integer :: np         -> number of vertices on the polygon
!     real(dp) :: rp, zp    -> r, z coordinates of the polygon's vertices
!
! Returns
!     logical :: cip  -> true if any of the corners lie within the polygon
!-------------------------------------------------------------------------------
function corner_in_polygon(cb, np, rp, zp) result(cip)

    use geometry, only: in_polygon

    implicit none

    type(cbrick), intent(IN) :: cb
    integer, intent(IN) :: np
    logical :: cip
    real(dp), dimension(np), intent(IN) :: rp, zp
    real(dp), dimension(4) :: rq, zq
    integer :: i

    rq = (/ cb%ri, cb%ro, cb%ro, cb%ri /)
    zq = (/ cb%zl, cb%zl, cb%zu, cb%zu /)

    do i = 1, 4
        if (in_polygon(np, rp, zp, rq(i), zq(i))) then
            cip = .true.
            return
        end if
    end do

    cip = .false.

end function corner_in_polygon

!-------------------------------------------------------------------------------
! corners_from_surf_xsect(cb, surf, phi, theta0) result(lvec)
!
! Determines the distances of the corners of a curved brick from a poloidal
! cross-section of a toroidal surface along a line that intersects the
! cross-section at a right angle in the poloidal plane.
!
! Input parameters:
!     type(cbrick) :: cb    -> structure for the curved brick
!     type(surface) :: surf -> structure for the toroidal surface
!     real(dp) :: phi       -> tor. angle of the pol. plane in which to evaluate
!     real(dp) :: theta0    -> initial guess, poloidal angle at which the 
!                              connecting line intersects the cross-section
!
! Returns:
!     real(dp), dimension(4) :: lvec  -> Distances of each of the four corners
!-------------------------------------------------------------------------------
function corners_from_surf_xsect(cb, surf, phi, theta0) result(lvec)

    use surface_calc,  only: surface, surf_r, surf_z
    use surface_solve, only: surf_perp_intersect_2d

    type(cbrick), intent(IN) :: cb
    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: phi, theta0
    real(dp) :: theta, rInt, zInt, chi2
    real(dp), dimension(4) :: lvec
    real(dp) :: l0

    l0 = sqrt((surf_r(surf, theta0, phi)-cb%cmr)**2 + &
              (surf_z(surf, theta0, phi)-cb%cmz)**2)

    call surf_perp_intersect_2d(surf, phi, cb%ri, cb%zl, l0, theta0, &
                                lvec(1), theta, rInt, zInt, chi2)
    call surf_perp_intersect_2d(surf, phi, cb%ro, cb%zl, l0, theta0, &
                                lvec(2), theta, rInt, zInt, chi2)
    call surf_perp_intersect_2d(surf, phi, cb%ro, cb%zu, l0, theta0, &
                                lvec(3), theta, rInt, zInt, chi2)
    call surf_perp_intersect_2d(surf, phi, cb%ri, cb%zu, l0, theta0, &
                                lvec(4), theta, rInt, zInt, chi2)

end function corners_from_surf_xsect

!-------------------------------------------------------------------------------
! corners_from_polygon(cb, np, rp, zp) result(lvec)
!
! Estimates the distances of the corners of a curved brick from a polygon in
! the r-z plane by finding the distance between each corner and the closest
! vertex of the polygon to the respective corner.
!
! Input parameters
!     type(cbrick) :: cb    -> curved brick whose corners are to be checked
!     integer :: np         -> number of vertices on the polygon
!     real(dp) :: rp, zp    -> r, z coordinates of the polygon's vertices
!
! Returns:
!     real(dp), dimension(4) :: lvec  -> Distances of each of the four corners
!-------------------------------------------------------------------------------
function corners_from_polygon(cb, np, rp, zp) result(lvec)

    use surface_solve, only: surf_perp_intersect_2d

    type(cbrick), intent(IN) :: cb
    integer, intent(IN) :: np
    real(dp), dimension(np), intent(IN) :: rp, zp
    real(dp), dimension(np) :: dists
    real(dp), dimension(4) :: lvec, rvec, zvec
    integer :: i, j

    rvec = (/ cb%ri, cb%ro, cb%ro, cb%ri /)
    zvec = (/ cb%zl, cb%zl, cb%zu, cb%zu /)

    do i = 1, 4
        do j = 1, np
            dists(j) = sqrt( (rvec(i)-rp(j))**2 + (zvec(i)-zp(j))**2 )
            if (j == 1 .or. lvec(i) > dists(j)) lvec(i) = dists(j)
        end do
    end do

end function corners_from_polygon

function estimate_cbrick_theta(cb, np, thetap, rp, zp) result(theta)

    implicit none

    type(cbrick), intent(IN) :: cb
    integer, intent(IN) :: np
    real(dp), dimension(np), intent(IN) :: thetap, rp, zp
    real(dp), dimension(np) :: dists
    real(dp) :: theta
    integer :: i, indMin

    do i = 1, np
        dists(i) = sqrt( (cb%cmr-rp(i))**2 + (cb%cmz-zp(i))**2 )
        if (i == 1 .or. dists(indMin) > dists(i)) indMin = i
    end do

    theta = thetap(indMin)

end function estimate_cbrick_theta

subroutine set_perpendicular(cb, surf, l0, theta0, phi0)

    use surface_calc,  only: surface
    use surface_solve, only: surf_perp_intersect_3d

    implicit none

    type(cbrick), intent(INOUT) :: cb
    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: l0, theta0, phi0
    real(dp) :: l, theta, phi, sx, sy, sz, ux, uy, uz, chi2

    call surf_perp_intersect_3d(surf, cb%cmx, cb%cmy, cb%cmz, &
                                l0, theta0, phi0, l, theta, phi, &
                                sx, sy, sz, ux, uy, uz, chi2)
    cb%axx = -ux
    cb%axy = -uy
    cb%axz = -uz

end subroutine set_perpendicular

!-------------------------------------------------------------------------------
! cbrick_grid_params(dr, dz, rMin, rMax, zMin, zMax, r0, z0, nr, nz)
!
! Determines the radial and vertical grid dimensions and starting points
! for an array of curved bricks. 
!
! If the upper and lower vertical limits zMin and zMax are on opposite sides
! of the z=0 plane, the cell boundaries will be placed such that one boundary
! is guaranteed to coincide with the z=0 plane.
!
! Input parameters:
!     real(dp) :: dr, dz      -> Radial, vertical dimension of each grid cell
!     real(dp) :: rMin, rMax  -> Min, max allowable radial grid extent
!     real(dp) :: zMin, zMax  -> Min, max allowable vertical grid extent
!
! Output parameters:
!     real(dp), :: r0, z0     -> Radial, vertical coordinates of first cells
!     integer :: nr, nz       -> Radial, vertical dimensions of grid
!-------------------------------------------------------------------------------
subroutine cbrick_grid_params(dr, dz, rMin, rMax, zMin, zMax, r0, z0, nr, nz)

    implicit none

    real(dp), intent(IN) :: dr, dz, rMin, rMax, zMin, zMax
    real(dp), intent(OUT) :: r0, z0
    integer, intent(OUT) :: nr, nz
    integer :: nzLower, nzUpper

    nr = floor((rMax-rMin)/dr)
    r0 = rMin

    ! If grid crosses the z=0 plane, ensure that the z=0 plane coincides with
    ! a vertical boundary between grid cells
    if (zMin < 0 .and. zMax > 0) then
        nzLower = floor(-zMin/dz)
        nzUpper = floor(zMax/dz)
        nz = nzLower + nzUpper
        z0 = -real(nzLower,dp)*dz
    else
        nz = floor((zMax-zMin)/dz)
        z0 = zMin
    endif

end subroutine cbrick_grid_params

!-------------------------------------------------------------------------------
! cbrick_obj_overlap(cb, maxCheckGap, ovl)
!
! Checks whether the faces of a cbrick intersect with any applicable object
! (port, coil, etc.) as specified by the user.
!
! NOTE: for correct functionality, calling procedures must call the
! intersections_to_check subroutine once prior to calls to cbrick_obj_overlap.
!
! Input parameters
!     type(cbrick) :: cb      -> cbrick whose faces are to be checked
!     real(dp) :: maxCheckGap -> Max spacing between adjacent query points on 
!                                each face for overlap checking
!
! Output parameters:
!     logical :: ovl          -> True if any faces are found to intersect an
!                                applicable object
!-------------------------------------------------------------------------------
subroutine cbrick_obj_overlap(cb, maxCheckGap, ovl)

    use intersections,  only: check_intersections

    implicit none

    type(cbrick), intent(IN) :: cb
    real(dp), intent(IN) :: maxCheckGap
    logical, intent(OUT) :: ovl
    integer :: nr, nz, nphi, i, j
    real(dp) :: dr, dz, dphi, rq, zq, phiq

    ovl = .false.
    
    ! Determine quantities and intervals for query points based on maxCheckGap
    nr = ceiling((cb%ro - cb%ri)/maxCheckGap) + 1
    nz = ceiling((cb%zu - cb%zl)/maxCheckGap) + 1
    nphi = ceiling(cb%ro*(cb%phif - cb%phib)/maxCheckGap) + 1
    dr = (cb%ro - cb%ri)/(real(nr,dp) - 1.0)
    dz = (cb%zu - cb%zl)/(real(nz,dp) - 1.0)
    dphi = (cb%phif - cb%phib)/(real(nphi,dp) - 1.0)

    ! Front and back faces: constant phi
    do i = 1, nr
        do j = 1, nz

            rq   = cb%ri + (i-1)*dr
            zq   = cb%zl + (j-1)*dz

            ! Back face
            call check_intersections(rq*cos(cb%phib), rq*sin(cb%phib), zq, ovl)
            if (ovl) return

            ! Front face
            call check_intersections(rq*cos(cb%phif), rq*sin(cb%phif), zq, ovl)
            if (ovl) return

        end do
    end do

    ! Lower and upper faces: constant z
    do i = 1, nr
        do j = 1, nphi

            rq   = cb%ri   + (i-1)*dr
            phiq = cb%phib + (j-1)*dphi

            ! Lower face
            call check_intersections(rq*cos(phiq), rq*sin(phiq), cb%zl, ovl)
            if (ovl) return

            ! Upper face
            call check_intersections(rq*cos(phiq), rq*sin(phiq), cb%zu, ovl)
            if (ovl) return

        end do
    end do

    ! Inboard and outboard faces: constant r
    do i = 1, nphi
        do j = 1, nz

            phiq = cb%phib + (i-1)*dphi
            zq   = cb%zl  + (j-1)*dz

            ! Inboard face
            call check_intersections(cb%ri*cos(phiq), cb%ri*sin(phiq), zq, ovl)
            if (ovl) return

            ! Outboard face
            call check_intersections(cb%ro*cos(phiq), cb%ro*sin(phiq), zq, ovl)
            if (ovl) return

        end do
    end do

end subroutine cbrick_obj_overlap

end module cbrick_properties


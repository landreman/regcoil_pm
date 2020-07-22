!-------------------------------------------------------------------------------
! base_grid.f90
!
! Subroutines for generating the toroidal grid used as a starting point to
! define and constrain the set of quadrilaterally-faced hexahedrons. The
! grid forms a torus in 3D space which encloses the plasma and vacuum vessel.
! The hexahedral magnets will be constrained bo lie in the space between the
! base grid and the vacuum vessel or other defined toroidal winding surface.
!
! When the grid points are projected onto the two-dimensional space of toroidal
! and poloidal angles, the quadrilaterals formed by the grid points (which do
! not necessarily have uniform spacing) each correspond to a single hexahedron
! or upper-/lower-hexahedron pair. Thus, an nxm grid corresponds to an array
! or (n-1)x(m-1) hexahedra or hexahedron pairs.
!
! Author: K. C. Hammond
! Contact: khammond@pppl.gov
! Updated: 2020-01-06
!-------------------------------------------------------------------------------
module base_grid

use magpie_globals, only: dp

implicit none

type base_grid_struct
    integer :: nTheta                   ! Number of points in poloidal dimension
    integer :: nPhi                     ! Number of points in toroidal dimension
    real(dp), allocatable :: theta(:,:) ! Poloidal angle of the vertex points
    real(dp), allocatable :: sep(:,:)   ! Separation distance from surface
    real(dp), allocatable :: phi(:,:)   ! Toroidal angles of each vertex
    real(dp), allocatable :: r(:,:)     ! r coordinate of each vertex
    real(dp), allocatable :: x(:,:)     ! x coordinate of each vertex
    real(dp), allocatable :: y(:,:)     ! y coordinate of each vertex
    real(dp), allocatable :: z(:,:)     ! z coordinate of each vertex
    real(dp), allocatable :: xSurf(:,:) ! x coord @ perp. projection on surface
    real(dp), allocatable :: ySurf(:,:) ! y coord @ perp. projection on surface
    real(dp), allocatable :: zSurf(:,:) ! z coord @ perp. projection on surface
    real(dp), allocatable :: thetaSurf(:,:) ! theta @ perp. projection on surf
    real(dp), allocatable :: phiSurf(:,:)   ! phi @ perp. projection on surf
    real(dp), allocatable :: nx(:,:)    ! x comp., unit surf. norm. @ projection
    real(dp), allocatable :: ny(:,:)    ! y comp., unit surf. norm. @ projection
    real(dp), allocatable :: nz(:,:)    ! z comp., unit surf. norm. @ projection
end type base_grid_struct

contains

!-------------------------------------------------------------------------------
! build_base_grid(surf, nTheta, nPhi, unif_poloidal_spacing, bg)
!
! Populates the arrays of coordinates for the base grid used for generating
! the magnet array. Also finds the locations and direction vectors of the 
! perpendicular intersections of each grid point with the input surface.
!
! Input parameters:
!     type(surface) :: surf    -> Structure containing the surface data
!     integer :: nTheta        -> Number of points in the poloidal array
!     real(dp) :: phi          -> Toroidal angles of the vertices
!     logical :: unif_poloidal_spacing  
!                              -> True if the values of theta should be
!                                 updated to have uniform spacing per unit
!                                 (input) poloidal angle. See docstring for
!                                 vertices_in_poloidal_plane for more info.
!
! Output parameters:
!     type(base_grid_struct) :: bg -> structure with grid data arrays
!-------------------------------------------------------------------------------
subroutine build_base_grid(surf, nTheta, nPhi, unif_poloidal_spacing, bg)

    use magpie_globals, only: qhex_grid_theta0, qhex_grid_theta1, &
                              qhex_grid_phi0, qhex_grid_phi1, &
                              radial_extent, qhex_cust_vert_theta, &
                              qhex_cust_vert_phi, qhex_cust_vert_sep, &
                              qhex_vert_theta, qhex_vert_phi, qhex_vert_sep
    use surface_calc,   only: surface
    use surface_solve,  only: surf_perp_intersect_3d

    implicit none

    type(surface), intent(IN) :: surf
    integer, intent(IN) :: nTheta, nPhi
    logical, intent(IN) :: unif_poloidal_spacing
    type(base_grid_struct), intent(OUT) :: bg
    real(dp), dimension(nTheta) :: rSurf, zSurf
    real(dp) :: l0 = 0.0, l, chi2
    integer :: i, j

    write(*,fmt='(A)') '    Building base grid'

    ! Allocate the arrays of the base_grid structure
    bg%nTheta = nTheta
    bg%nPhi   = nPhi
    allocate( bg%theta(nTheta,nPhi), bg%phi(nTheta,nPhi),             &
              bg%sep(nTheta,nPhi), bg%r(nTheta,nPhi), bg%x(nTheta,nPhi),   &
              bg%y(nTheta,nPhi), bg%z(nTheta,nPhi), bg%xSurf(nTheta,nPhi), &
              bg%ySurf(nTheta,nPhi), bg%zSurf(nTheta,nPhi),                &
              bg%thetaSurf(nTheta,nPhi), bg%phiSurf(nTheta,nPhi),          &
              bg%nx(nTheta,nPhi), bg%ny(nTheta,nPhi), bg%nz(nTheta,nPhi)    )

    ! Calculate automatic values of sufr_theta, phi, and sep
    call generate_base_grid_angles(nTheta, nPhi, &
                                   qhex_grid_theta0, qhex_grid_theta1, &
                                   qhex_grid_phi0,   qhex_grid_phi1,   &
                                   bg%theta, bg%phi)
    bg%sep = radial_extent

    ! Overwrite with custom input values if required
    if (qhex_cust_vert_theta) bg%theta = qhex_vert_theta(1:nTheta,1:nPhi)
    if (qhex_cust_vert_phi)   bg%phi   = qhex_vert_phi(1:nTheta,1:nPhi)
    if (qhex_cust_vert_sep)   bg%sep   = qhex_vert_sep(1:nTheta,1:nPhi)

    ! Find the r and z coordinates of each grid point, adjusting if necessary
    do i = 1, nPhi
        call vertices_in_poloidal_plane( &
                 surf, nTheta, bg%phi(:,i), unif_poloidal_spacing, &
                 bg%theta(:,i), bg%sep(:,i), bg%r(:,i), bg%z(:,i), &
                 rSurf, zSurf)
    end do

    ! Transform polar coordinates to Cartesian x and y coordinates
    bg%x = bg%r*cos(bg%phi)
    bg%y = bg%r*sin(bg%phi)

    ! Find the points of perpendicular intersection with the surface
    do i = 1, nTheta
        do j = 1, nPhi
            call surf_perp_intersect_3d(surf, bg%x(i,j), bg%y(i,j), bg%z(i,j), &
                     l0, bg%theta(i,j), bg%phi(i,j),                           &
                     l, bg%thetaSurf(i,j), bg%phiSurf(i,j),                    &
                     bg%xSurf(i,j), bg%ySurf(i,j), bg%zSurf(i,j),              &
                     bg%nx(i,j), bg%ny(i,j), bg%nz(i,j), chi2             )
        end do
    end do

end subroutine build_base_grid

!-------------------------------------------------------------------------------
! vertices_in_poloidal_plane(surf, nTheta, phi, unif_spacing, 
!                            theta, sep, r, z, rSurf, zSurf)
!
! Determines the coordinates of the base grid vertices in a poloidal plane
! according to input poloidal angles (theta) and separation distances (sep)
! from an input surface (e.g., representing the vacuum vessel).
!
! This subroutine has the option to modify the surface poloidal angles in order
! to attain uniform spacing between successive poloidal vertices per unit
! input poloidal angle. Thus, if the input poloidal angles are all uniformly
! spaced, the updates poloidal angles will not necessarily be uniformly spaced
! but rather chosen such that the associated vertices will have equal spacing
! in the poloidal plane.
!
! Note: the set of points in the "poloidal plane" need not have the same 
! value of toroidal angle phi, although this is the most common usage case.
!
! Input parameters:
!     type(surface) :: surf    -> Structure containing the surface data
!     integer :: nTheta        -> Number of points in the poloidal array
!     real(dp) :: phi          -> Toroidal angles of the vertices
!     logical :: unif_spacing  -> True if the values of theta should be
!                                 updated to have uniform spacing per unit
!                                 (input) poloidal angle. Sep will also be
!                                 updated through interpolation in case it
!                                 varies with poloidal angle.
!
! Updated parameters:
!     real(dp) :: theta   -> Surface theta values of the vertex points
!     real(dp) :: sep          -> Separation distance of vertex from vessel
!
! Output parameters:
!     real(dp) :: r, z         -> r and z coordinates of the output vertex pts  
!     real(dp) :: rSurf, zSurf -> r and z coordinates of corresponding locations
!                                 on the vacuum vessel
!-------------------------------------------------------------------------------
subroutine vertices_in_poloidal_plane(surf, nTheta, phi, unif_spacing, &
                                      theta, sep, r, z, rSurf, zSurf)

    use magpie_globals, only: maxIter, dist_tol
    use algebra,        only: linear_interpolate
    use geometry,       only: point_segment_lengths
    use surface_calc,   only: surface

    implicit none

    type(surface), intent(IN) :: surf
    integer, intent(IN) :: nTheta
    real(dp), dimension(nTheta), intent(IN) :: phi
    logical, intent(IN) :: unif_spacing
    real(dp), dimension(nTheta), intent(INOUT) :: theta, sep
    real(dp), dimension(nTheta), intent(OUT)   :: r, z, rSurf, zSurf
    real(dp), dimension(nTheta) :: theta_in, theta_prev, sep_in
    real(dp), dimension(nTheta) :: s_ideal, s
    integer :: i, j
    logical :: adjust_vertex_spacing
    real(dp) :: meanSqResid, drdt, dzdt, alpha, nr, nz
    
    ! Ideal normalized segment length (proportional to user-input angle)
    s_ideal = (theta - theta(1))/(theta(nTheta) - theta(1))

    ! Store initial values of vessel theta and segment length
    theta_in = theta
    sep_in = sep

    ! Iteratively adjust theta until the uniform-spacing condition is
    ! attained (if desired)
    do i = 1, maxIter

        ! Update the r and z coordinates of the vertex points
        do j = 1, nTheta
            call thetaSep_to_rz(surf, theta(j), phi(j), sep(j), &
                                r(j), z(j), rSurf(i), zSurf(i))
        end do

        ! Determine the variance in arclength/angle ratio between segments
        ! according to the user-input vessel theta angles
        call point_segment_lengths(nTheta, s_ideal, r, z, s, meanSqResid)

        if (isnan(meanSqResid)) then
            stop 'ves_rz_in_poloidal_plane: NaN encountered in segment ' // &
                 'length resituals. Check values of qhex_vert_theta and ' // &
                 'nVertices.'
        end if

        ! Terminate if adequately converged or if not enforcing equal ratios
        if (meanSqResid < dist_tol**2 .or. .not. unif_spacing) exit

        ! Reset the surface poloidal angles corresponding to the vertex points
        theta_prev = theta
        call linear_interpolate(nTheta, s,       theta_prev, &
                                nTheta, s_ideal, theta)

        ! Update the radial separation for the vertex points
        call linear_interpolate(nTheta, theta_in, sep_in,    &
                                nTheta, theta,    sep)

    end do

end subroutine vertices_in_poloidal_plane

!-------------------------------------------------------------------------------
! generate_base_grid_angles(ntheta, nphi, theta0, theta1, phi0, phi1, 
!                           theta, phi)
!
! Populates two-dimensional arrays with regularly-spaced values for theta and
! phi angles, with theta varying along the first dimension (rows) and phi 
! varying along the second dimension (columns).
!
! This essentially implements the functionality of the "linspace" and 
! "meshgrids" functions in Matlab, NumPy, etc.
!
! Input parameters:
!     integer :: ntheta, nphi    -> Dimensions of the grid in theta and phi
!     real(dp) :: theta0, theta1 -> Start and end points of the grid in theta
!     real(dp) :: phi0, phi1     -> Start and end points of the grid in phi
!
! Output parameters (must be arrays of dimension ntheta, nphi):
!     real(dp) :: theta(:,:)     -> Array of theta values (uniform across rows)
!     real(dp) :: phi(:,:)       -> Array of phi values (uniform down columns)
!-------------------------------------------------------------------------------
subroutine generate_base_grid_angles(ntheta, nphi, theta0, theta1, phi0, phi1, &
                                     theta, phi)

    implicit none

    integer, intent(IN) :: ntheta, nphi
    real(dp), intent(IN) :: theta0, theta1, phi0, phi1
    real(dp), dimension(ntheta, nphi), intent(OUT) :: theta, phi
    integer :: i, j
    real(dp) :: dtheta, dphi

    dtheta = (theta1 - theta0)/(ntheta - 1)
    dphi   = (phi1   - phi0  )/(nphi   - 1)

    do i = 1, ntheta
        do j = 1, nphi
            theta(i,j) = theta0 + (i-1)*dtheta
            phi(i,j)   = phi0   + (j-1)*dphi
        end do
    end do

end subroutine generate_base_grid_angles

!-------------------------------------------------------------------------------
! thetaSep_to_rz(surf, theta, phi, sep, r, z)
!
! Calculates the poloidal r and z coordinates of points separated from the
! vessel at a specified distance along the (2d, poloidal) normal vector 
! originating from a specificied location on the vessel.
!
! Input parameters:
!     type(surface) :: surf    -> Structure with the surface data
!     real(dp) :: theta, phi   -> Poloidal and toroidal angle of the reference
!                                 point on the surface
!     real(dp) :: sep          -> Separation of the output point from the 
!                                 reference point along the surface normal
!
! Output parameters:
!     real(dp) :: r, z         -> Poloidal r, z coordinates of the output point
!     real(dp) :: rSurf, zSurf -> Poloidal r, z coordinates evaluated on the
!                                 surface for the input theta and phi
!-------------------------------------------------------------------------------
subroutine thetaSep_to_rz(surf, theta, phi, sep, r, z, rSurf, zSurf)

    use surface_calc, only: surface, surf_r, surf_z, surf_poloidal_unorm

    implicit none

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi, sep
    real(dp), intent(OUT) :: r, z, rSurf, zSurf
    real(dp) :: ur, uz

    ! Coordinates and normal vector of the surface at the desired angle
    rSurf = surf_r(surf, theta, phi)
    zSurf = surf_z(surf, theta, phi)
    call surf_poloidal_unorm(surf, theta, phi, ur, uz)

    ! Location displaced outward from the surface along the normal direction
    r = rSurf + sep*ur
    z = zSurf + sep*uz

end subroutine thetaSep_to_rz

end module base_grid

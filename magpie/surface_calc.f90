!-------------------------------------------------------------------------------
! surface_calc.f90
!
! Module containing functions for computing spatial coordinates and derivatives
! on a toroidal surface parametrized by a poloidal angle theta and a toroidal
! (azimuthal) angle phi. The surface is assumed to be represented as a Fourier 
! series whose radial (r) and vertical (z) coordinates at a given azimuthal 
! angle phi are sums of sinusoidal modes characterized by
!
!   r = sum_{m,n} rc(m,n) * cos(m*theta + nfp*n*phi) +
!                 rs(m,n) * sin(m*theta + nfp*n*phi)
!   z = sum_{m,n} zs(m,n) * sin(m*theta + nfp*n*phi) + 
!                 zc(m,n) * cos(m*theta + nfp*n*phi)
!     
! Here, nfp is the number of field periods; i.e., how many times the geometry
! of the surface repeats itself during one toroidal revolution. Also, if the 
! surface has stellarator symmetry, then rs(m,n) = zc(m,n) = 0 for all m,n.
!
! The rc, rs, zc, and zs coefficients are each stored in one-dimensional arrays
! in a structure of type "surface" defined in this module. Accompanying integer
! arrays provide the respective m and n numbers for each mode. 
!
! Author: K. C. Hammond
! Email:  khammond@pppl.gov
! Last updated: 2020-01-05
!-------------------------------------------------------------------------------
module surface_calc

use magpie_globals, only: dp

implicit none

! Type definition for a structure containing parameters for a toroidal surface
! NOTE: phase for m,n mode is assumed to be m*theta + nfp*n*phi
type surface
    integer :: nModes                   ! Total number of modes
    integer :: nfp                      ! Number of field periods
    logical :: symm = .false.           ! True if stellarator-symmetric
    real(dp), allocatable :: m(:)       ! Poloidal number for each mode
    real(dp), allocatable :: n(:)       ! Toroidal number for each mode
    real(dp), allocatable :: rc(:)      ! Cosine coefficient, r coordinate
    real(dp), allocatable :: zs(:)      ! Sine   coefficient, z coordinate
    real(dp), allocatable :: rs(:)      ! Sine   coefficient, r coordinate
    real(dp), allocatable :: zc(:)      ! Cosine coefficient, z coordinate
end type surface

contains

! R coordinate at a given theta and phi angle
real(dp) function surf_r(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_r = 0
    do i = 1, surf%nModes
        phase = surf%m(i) * theta + surf%n(i) * surf%nfp * phi
        surf_r = surf_r + surf%rc(i)*cos(phase) 
        if (.not. surf%symm) surf_r = surf_r + surf%rs(i)*sin(phase)
    end do

end function surf_r

! Z coordinate at a given theta and phi angle
real(dp) function surf_z(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_z = 0
    do i = 1, surf%nModes
        phase = surf%m(i) * theta + surf%n(i) * surf%nfp * phi
        surf_z = surf_z + surf%zs(i)*sin(phase) 
        if (.not. surf%symm) surf_z = surf_z + surf%zc(i)*cos(phase)
    end do

end function surf_z

! Derivative of the R coordinate wrt theta at a given theta and phi angle
real(dp) function surf_drdt(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_drdt = 0
    do i = 1, surf%nModes
        phase = surf%m(i) * theta + surf%n(i) * surf%nfp * phi
        surf_drdt = surf_drdt - surf%m(i)*surf%rc(i)*sin(phase) 
        if (.not. surf%symm) &
            surf_drdt = surf_drdt + surf%m(i)*surf%rs(i)*cos(phase)
    end do

end function surf_drdt

! Derivative of the Z coordinate wrt theta at a given theta and phi angle
real(dp) function surf_dzdt(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_dzdt = 0
    do i = 1, surf%nModes
        phase = surf%m(i) * theta + surf%n(i) * surf%nfp * phi
        surf_dzdt = surf_dzdt + surf%m(i)*surf%zs(i)*cos(phase) 
        if (.not. surf%symm) &
            surf_dzdt = surf_dzdt - surf%m(i)*surf%zc(i)*sin(phase)
    end do

end function surf_dzdt

! Derivative of the R coordinate wrt phi at a given theta and phi angle
real(dp) function surf_drdp(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_drdp = 0
    do i = 1, surf%nModes
        phase = surf%m(i) * theta + surf%n(i) * surf%nfp * phi
        surf_drdp = surf_drdp - surf%n(i)*surf%rc(i)*sin(phase) 
        if (.not. surf%symm) &
            surf_drdp = surf_drdp + surf%n(i)*surf%rs(i)*cos(phase)
    end do

    surf_drdp = surf_drdp * surf%nfp

end function surf_drdp

! Derivative of the Z coordinate wrt phi at a given theta and phi angle
real(dp) function surf_dzdp(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_dzdp = 0
    do i = 1, surf%nModes
        phase = surf%m(i)*theta + surf%n(i)*surf%nfp*phi
        surf_dzdp = surf_dzdp + surf%n(i)*surf%zs(i)*cos(phase) 
        if (.not. surf%symm) &
            surf_dzdp = surf_dzdp - surf%n(i)*surf%zc(i)*sin(phase)
    end do

    surf_dzdp = surf_dzdp * surf%nfp

end function surf_dzdp

! Second derivative of the R coordinate wrt theta at a given theta and phi angle
real(dp) function surf_d2rdt2(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_d2rdt2 = 0
    do i = 1, surf%nModes
        phase = surf%m(i) * theta + surf%n(i) * surf%nfp * phi
        surf_d2rdt2 = surf_d2rdt2 - surf%m(i)**2*surf%rc(i)*cos(phase) 
        if (.not. surf%symm) &
            surf_d2rdt2 = surf_d2rdt2 - surf%m(i)**2*surf%rs(i)*sin(phase)
    end do

end function surf_d2rdt2

! Second derivative of the Z coordinate wrt theta at a given theta and phi angle
real(dp) function surf_d2zdt2(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_d2zdt2 = 0
    do i = 1, surf%nModes
        phase = surf%m(i) * theta + surf%n(i) * surf%nfp * phi
        surf_d2zdt2 = surf_d2zdt2 - surf%m(i)**2*surf%zs(i)*sin(phase) 
        if (.not. surf%symm) &
            surf_d2zdt2 = surf_d2zdt2 - surf%m(i)*surf%zc(i)*cos(phase)
    end do

end function surf_d2zdt2

! Second derivative of the R coordinate wrt phi at a given theta and phi angle
real(dp) function surf_d2rdp2(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_d2rdp2 = 0
    do i = 1, surf%nModes
        phase = surf%m(i) * theta + surf%n(i) * surf%nfp * phi
        surf_d2rdp2 = surf_d2rdp2 - surf%n(i)**2*surf%rc(i)*cos(phase) 
        if (.not. surf%symm) &
            surf_d2rdp2 = surf_d2rdp2 - surf%n(i)**2*surf%rs(i)*sin(phase)
    end do

    surf_d2rdp2 = surf_d2rdp2 * surf%nfp**2

end function surf_d2rdp2

! Second derivative of the Z coordinate wrt phi at a given theta and phi angle
real(dp) function surf_d2zdp2(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_d2zdp2 = 0
    do i = 1, surf%nModes
        phase = surf%m(i)*theta + surf%n(i)*surf%nfp*phi
        surf_d2zdp2 = surf_d2zdp2 - surf%n(i)**2*surf%zs(i)*sin(phase) 
        if (.not. surf%symm) &
            surf_d2zdp2 = surf_d2zdp2 - surf%n(i)**2*surf%zc(i)*cos(phase)
    end do

    surf_d2zdp2 = surf_d2zdp2 * surf%nfp**2

end function surf_d2zdp2

! Cross-derivative of the R coordinate wrt phi & theta at a given theta and phi 
real(dp) function surf_d2rdtdp(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_d2rdtdp = 0
    do i = 1, surf%nModes
        phase = surf%m(i) * theta + surf%n(i) * surf%nfp * phi
        surf_d2rdtdp = surf_d2rdtdp - surf%m(i)*surf%n(i)*surf%rc(i)*cos(phase) 
        if (.not. surf%symm) &
            surf_d2rdtdp = &
                surf_d2rdtdp - surf%m(i)*surf%n(i)*surf%rs(i)*sin(phase)
    end do

    surf_d2rdtdp = surf_d2rdtdp * surf%nfp

end function surf_d2rdtdp

! Cross-derivative of the Z coordinate wrt phi & theta at a given theta and phi 
real(dp) function surf_d2zdtdp(surf, theta, phi)

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    integer :: i
    real(dp) :: phase

    surf_d2zdtdp = 0
    do i = 1, surf%nModes
        phase = surf%m(i)*theta + surf%n(i)*surf%nfp*phi
        surf_d2zdtdp = surf_d2zdtdp - surf%m(i)*surf%n(i)*surf%zs(i)*sin(phase) 
        if (.not. surf%symm) &
            surf_d2zdtdp = &
                surf_d2zdtdp - surf%m(i)*surf%n(i)*surf%zc(i)*cos(phase)
    end do

    surf_d2zdtdp = surf_d2zdtdp * surf%nfp

end function surf_d2zdtdp

! Normal vector at a given location (NOT a unit vector)
subroutine surf_normal(surf, theta, phi, nx, ny, nz)

    use geometry, only: cross_prod

    implicit none

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    real(dp), intent(OUT) :: nx, ny, nz
    real(dp) :: r, drdt, dxdt, dydt, dzdt, drdp, dxdp, dydp, dzdp

    ! Calculate required coordinates and derivatives for r and z
    r = surf_r(surf, theta, phi)
    drdt = surf_drdt(surf, theta, phi)
    dzdt = surf_dzdt(surf, theta, phi)
    drdp = surf_drdp(surf, theta, phi)
    dzdp = surf_dzdp(surf, theta, phi)

    ! Transform r derivatives to x and y derivatives
    dxdt = drdt*cos(phi)
    dydt = drdt*sin(phi)
    dxdp = drdp*cos(phi) - r*sin(phi)
    dydp = drdp*sin(phi) + r*cos(phi)
    
    ! Calculate the normal vector
    call cross_prod(dxdp, dydp, dzdp, dxdt, dydt, dzdt, nx, ny, nz)

end subroutine surf_normal

! Unit normal vector at a given location
subroutine surf_unorm(surf, theta, phi, ux, uy, uz)

    use geometry, only: unit_vector

    implicit none

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    real(dp), intent(OUT) :: ux, uy, uz
    real(dp) :: nx, ny, nz

    call surf_normal(surf, theta, phi, nx, ny, nz)
    call unit_vector(nx, ny, nz, ux, uy, uz)

end subroutine surf_unorm

! Unit normal vector to a poloidal cross-section of a toroidal surface
subroutine surf_poloidal_unorm(surf, theta, phi, ur, uz)

    implicit none

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: theta, phi
    real(dp), intent(OUT) :: ur, uz
    real(dp) :: drdt, dzdt, angle

    drdt = surf_drdt(surf, theta, phi)
    dzdt = surf_dzdt(surf, theta, phi)
    angle = atan2(dzdt, drdt)
    ur = sin(angle)
    uz = -cos(angle)

end subroutine surf_poloidal_unorm

end module surface_calc


!-------------------------------------------------------------------------------
! surface_solve.f90
!
! Contains a few iterative Newton-Raphson solvers for various systems of 
! nonlinear equations useful for addressing geometric questions about toroidal 
! surfaces.
!
! Author:  K. C. Hammond
! Contact: khammod@pppl.gov
! Updated: 2019-12-30
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Module with type definitions for the data structures used by the solver 
! routines to store values that are used for calculating residuals
!-------------------------------------------------------------------------------
module surface_solve_types

use magpie_globals, only: dp
use surface_calc,   only: surface

implicit none

!-------------------------------------------------------------------------------
! Constant variables supplied as inputs (used as needed depending on the solver)
!-------------------------------------------------------------------------------
type input_vars
    real(dp) :: ox, oy, oz, or    ! x, y, z coords, input reference point
    real(dp) :: ax, ay, az, ar    ! x, y, z components, input ref. axis
    real(dp) :: phi               ! toroidal angle (if held fixed)
    type(surface) :: surf         ! surface for evaluation
end type input_vars

!-------------------------------------------------------------------------------
! Variables that update on each iteration (used as needed depending on the 
! solver)
!-------------------------------------------------------------------------------
type updating_vars
    real(dp) :: l, theta, phi            ! length/angle variables to solve for
    real(dp) :: fx, fy, fz, fr           ! residual components
    real(dp) :: sx, sy, sz, sr           ! x, y, z, and r coordinates on surface
    real(dp) :: nx, ny, nz, nr           ! components of surface normal
    real(dp) :: ln                       ! length of surface normal
    real(dp) :: ux, uy, uz               ! components of surface unit normal
    real(dp) :: drdt, drdp, dzdt, dzdp   ! partial derivs wrt angles on surf
    real(dp) :: sl, st, sp               ! l, theta, phi comps of Newton step
    real(dp) :: chi2                     ! sum of squares of residuals
    real(dp) :: grad_l, grad_t, grad_p   ! components of chi2 gradient
    real(dp) :: slope                    ! deriv of chi2 in Newton direction
end type updating_vars

end module surface_solve_types

!-------------------------------------------------------------------------------
! Main module with solvers and associated subroutines
!-------------------------------------------------------------------------------
module surface_solve

use surface_solve_types, only: input_vars, updating_vars
use magpie_globals,      only: dp
use surface_calc,        only: surface

implicit none

interface
    subroutine calcF(inp, upd)
        use surface_solve_types, only: input_vars, updating_vars
        implicit none
        type(input_vars), intent(IN) :: inp
        type(updating_vars), intent(INOUT) :: upd
    end subroutine calcF
end interface

contains

!-------------------------------------------------------------------------------
! surf_dist_3d(surf, ox, oy, oz, ax, ay, az, l0, theta0, phi0, 
!              l, theta, phi, x, y, z)
!
! Computes the distance between a point and a surface along a particular
! line in 3D space, using the Newton-Raphson method.
!
! Input parameters:
!     type(surface) :: surf   -> structure with the surface parameters
!     real(dp) :: ox, oy, oz  -> x, y, and z coordinates, reference point
!     real(dp) :: ax, ay, az  -> x, y, and z components, unit vector defining
!                                the direction from test point to the surface
!     real(dp) :: l0, theta0, 
!                 phi0        -> initial guesses of the distance l, as well as
!                                the theta and phi coordinates of the 
!                                intersection point on the surface
!
! Return parameters:
!     real(dp) :: l           -> distance from the point to the vessel
!     real(dp) :: theta, phi  -> theta and phi coordinates, intersection point
!     real(dp) :: sx, sy, sz  -> x, y, and z coordinates, intersection point
!-------------------------------------------------------------------------------
subroutine surf_dist_3d(surf, ox, oy, oz, ax, ay, az, l0, theta0, phi0, &
                        l, theta, phi, sx, sy, sz, chi2)

    use magpie_globals, only: dist_tol, maxIter
    use algebra,        only: inverse_3x3
    use surface_calc,   only: surf_drdt, surf_dzdt, surf_drdp, surf_dzdp

    implicit none

    type(surface), intent(IN) :: surf
    real(dp), intent(IN) :: ox, oy, oz, ax, ay, az, l0, theta0, phi0
    real(dp), intent(OUT) :: l, theta, phi, sx, sy, sz, chi2
    procedure(calcF), pointer :: calcF_func => null()
    type(input_vars) :: inp
    type(updating_vars) :: upd
    integer :: i
    real(dp) :: sinPhi, cosPhi
    real(dp) ::  jac_xl,  jac_xt,  jac_xp, &
                 jac_yl,  jac_yt,  jac_yp, &
                 jac_zl,  jac_zt,  jac_zp
    real(dp) :: ijac_lx, ijac_ly, ijac_lz, &
                ijac_px, ijac_py, ijac_pz, &
                ijac_tx, ijac_ty, ijac_tz

    ! Initialize the constant (input) variables in structure
    inp%ox = ox; inp%oy = oy; inp%oz = oz
    inp%ax = ax; inp%ay = ay; inp%az = az
    inp%surf = surf

    ! Initialize the updating variables in the relevant structure
    upd%l = l0
    upd%theta = theta0
    upd%phi = phi0

    ! constant components of the Jacobian (first column)
    jac_xl = -ax
    jac_yl = -ay
    jac_zl = -az

    ! Calculate fx, fy, fz (should be zero if a solution is found)
    ! as well as coordinates on the vessel at (theta, phi)
    calcF_func => calcF_surf_dist_3d
    call calcF_func(inp, upd)

    do i = 1, maxIter+1

        ! terminate if solution has sufficiently converged
        upd%chi2 = upd%fx**2 + upd%fy**2 + upd%fz**2

        ! for debugging
        !write(*,fmt='(A, I0, A, F8.4, A, F8.4, A, F8.4, A, ES10.3)') &
        !    '    Iter ', i, ':  l = ', l, &
        !    '; theta = ', upd%theta*180./3.14159, &
        !    '; phi = ', upd%phi*180./3.14159, '; chi2 = ', upd%chi2

        if (upd%chi2 < dist_tol**2 .or. i == maxIter) then
            ! Update the values of the output variables
            l = upd%l
            theta = upd%theta
            phi = upd%phi
            sx = upd%sx
            sy = upd%sy
            sz = upd%sz
            chi2 = upd%chi2
            exit
        end if

        ! Update sinusoidal terms
        sinPhi = sin(upd%phi)
        cosPhi = cos(upd%phi)
        
        ! compute the Jacobian, second column
        upd%drdt = surf_drdt(inp%surf, upd%theta, upd%phi)
        upd%dzdt = surf_dzdt(inp%surf, upd%theta, upd%phi)
        jac_xt = upd%drdt * cosPhi
        jac_yt = upd%drdt * sinPhi
        jac_zt = upd%dzdt

        ! compute the Jacobian, third column
        upd%drdp = surf_drdp(inp%surf, upd%theta, upd%phi)
        upd%dzdp = surf_dzdp(inp%surf, upd%theta, upd%phi)
        jac_xp = upd%drdp * cosPhi - upd%sr * sinPhi
        jac_yp = upd%drdp * sinPhi + upd%sr * cosPhi
        jac_zp = upd%dzdp

        ! invert the Jacobian
        call inverse_3x3(jac_xl, jac_xt, jac_xp, jac_yl, jac_yt, jac_yp, &
                         jac_zl, jac_zt, jac_zp, ijac_lx, ijac_ly, ijac_lz, &
                         ijac_tx, ijac_ty, ijac_tz, ijac_px, ijac_py, ijac_pz)

        ! The Newton step
        upd%sl = -ijac_lx * upd%fx - ijac_ly * upd%fy - ijac_lz * upd%fz
        upd%st = -ijac_tx * upd%fx - ijac_ty * upd%fy - ijac_tz * upd%fz
        upd%sp = -ijac_px * upd%fx - ijac_py * upd%fy - ijac_pz * upd%fz

        ! The gradient of chi2 in the Newton direction
        upd%grad_l = upd%fx*jac_xl + upd%fy*jac_yl + upd%fz*jac_zl
        upd%grad_t = upd%fx*jac_xt + upd%fy*jac_yt + upd%fz*jac_zt
        upd%grad_p = upd%fx*jac_xp + upd%fy*jac_yp + upd%fz*jac_zp
        upd%slope = upd%grad_l*upd%sl + upd%grad_t*upd%st &
                        + upd%grad_p*upd%sp 

        ! Update l, theta, and phi, refining the Newton step if necessary
        call updateStep(calcF_func, inp, upd)


    end do

end subroutine surf_dist_3d

subroutine calcF_surf_dist_3d(inp, upd)

    use surface_calc, only: surface, surf_r, surf_z

    implicit none

    type(input_vars), intent(IN) :: inp
    type(updating_vars), intent(INOUT) :: upd

    ! r, x, y, and z coordinates on the vessel corresponding to theta and phi
    upd%sr = surf_r(inp%surf, upd%theta, upd%phi)
    upd%sx = upd%sr * cos(upd%phi)
    upd%sy = upd%sr * sin(upd%phi)
    upd%sz = surf_z(inp%surf, upd%theta, upd%phi)

    ! Components that should ideally be zero
    upd%fx = upd%sx - (inp%ox + upd%l*inp%ax)
    upd%fy = upd%sy - (inp%oy + upd%l*inp%ay)
    upd%fz = upd%sz - (inp%oz + upd%l*inp%az)

end subroutine calcF_surf_dist_3d

!-------------------------------------------------------------------------------
! surf_dist_2d(surf, phi, or, oz, ar, az, l0, theta0, l, theta, r, z, chi2)
!
! Computes the distance between a point and a toroidal surface along a 
! particular line in a 2D (poloidal, constant-phi) cross-section
!
! Input parameters:
!     type(surface) :: surf   -> structure with the surface parameters
!     real(dp) :: phi         -> toroidal angle (radians) of the cross-section
!     real(dp) :: or, oz      -> r and z coordinates (meters) of the ref. point
!     real(dp) :: ar, az      -> r and z components of a unit vector defining
!                                the direction from the test point to the vessel
!     real(dp) :: l0, theta0  -> initial guesses of the distance l, as well as
!                                the theta coordinate of the intersection point
!                                on the vessel
!
! Return parameters:
!     real(dp) :: l           -> distance from the point to the vessel
!     real(dp) :: theta       -> theta coordinate of the intersection point
!     real(dp) :: sr, sz      -> r and z coordinates of the intersection point
!     real(dp) :: chi2        -> sum of squares of residuals for the solution
!-------------------------------------------------------------------------------
subroutine surf_dist_2d(surf, phi, or, oz, ar, az, l0, theta0, &
                        l, theta, sr, sz, chi2)

    use magpie_globals, only: dist_tol, maxIter
    use algebra,        only: inverse_2x2
    use surface_calc,   only: surface, surf_r, surf_z, surf_drdt, surf_dzdt

    implicit none

    type(surface), intent(IN) :: surf
    real(dp), intent(IN)  :: phi, or, oz, ar, az, l0, theta0
    real(dp), intent(OUT) :: l, theta, sr, sz, chi2
    type(input_vars) :: inp
    type(updating_vars) :: upd
    integer :: i
    !real(dp) :: fr, fz, sl, st
    real(dp) :: drdt, dzdt
    real(dp) ::  jac_rl,  jac_rt,  jac_zl,  jac_zt
    real(dp) :: ijac_lr, ijac_lz, ijac_tr, ijac_tz

    ! Initialize relevant parameters in the input variable structure
    inp%ar = ar
    inp%az = az
    inp%or = or
    inp%oz = oz
    inp%phi = phi
    inp%surf = surf

    ! Initialize values for which to solve
    upd%l     = l0
    upd%theta = theta0

    ! constant components of the Jacobian (first column)
    jac_rl = -ar
    jac_zl = -az

    do i = 1, maxIter+1

        upd%sr = surf_r(inp%surf, upd%theta, inp%phi)
        upd%sz = surf_z(inp%surf, upd%theta, inp%phi)

        ! components that should ideally be zero
        upd%fr = upd%sr - (inp%or + upd%l*inp%ar)
        upd%fz = upd%sz - (inp%oz + upd%l*inp%az)

        ! terminate if solution has sufficiently converged
        upd%chi2 = upd%fr**2 + upd%fz**2
        if (upd%chi2 < dist_tol**2 .or. i == maxIter) then
            l = upd%l
            theta = upd%theta
            sr = upd%sr
            sz = upd%sz
            chi2 = upd%chi2
            exit
        end if
        
        ! compute the Jacobian, second column
        drdt = surf_drdt(inp%surf, upd%theta, inp%phi)
        dzdt = surf_dzdt(inp%surf, upd%theta, inp%phi)
        jac_rt = drdt 
        jac_zt = dzdt

        ! invert the Jacobian
        call inverse_2x2( jac_rl,  jac_rt,  jac_zl,  jac_zt, &
                         ijac_lr, ijac_lz, ijac_tr, ijac_tz )

        ! The Newton step
        upd%sl = -ijac_lr * upd%fr - ijac_lz * upd%fz
        upd%st = -ijac_tr * upd%fr - ijac_tz * upd%fz

        ! Update the parameters
        upd%l     = upd%l     + upd%sl
        upd%theta = upd%theta + upd%st

    end do

end subroutine surf_dist_2d

!-------------------------------------------------------------------------------
! surf_perp_intersect_3d(surf, ox, oy, oz, l0, theta0, phi0, 
!                       l, theta, phi, x, y, z, nx, ny, nz)
!
! Determines the parameters of a line segment perpendicular to a smooth toroidal
! surface that intersects a given query point.
!
! Input parameters:
!     type(surface) :: surf   -> structure with the surface parameters
!     real(dp) :: ox, oy, oz  -> x, y, and z coordinates of the query point
!     real(dp) :: l0, theta0, 
!                phi0        -> initial guesses of the distance l, as well as
!                               the theta and phi coordinates of the 
!                               intersection point on the surface
!
! Return parameters:
!     real(dp) :: l           -> distance from the point to the surface
!     real(dp) :: theta, phi  -> theta and phi coordinates, intersection point
!     real(dp) :: sx, sy, sz  -> x, y, and z coordinates, intersection point
!     real(dp) :: ux, uy, uz  -> x, y, and z components, unit normal vector
!                                at the surface intersection point
!-------------------------------------------------------------------------------
subroutine surf_perp_intersect_3d(surf, ox, oy, oz, l0, theta0, phi0, &
                                  l, theta, phi, sx, sy, sz, ux, uy, uz, chi2)

    use magpie_globals,      only: dist_tol, maxIter
    use algebra,             only: inverse_3x3
    use surface_calc,        only: surface, &
                                   surf_d2rdt2, surf_d2rdp2, surf_d2rdtdp, &
                                   surf_d2zdt2, surf_d2zdp2, surf_d2zdtdp
    use surface_solve_types, only: input_vars, updating_vars

    implicit none

    type(surface), intent(IN) :: surf
    real(dp), intent(IN)  :: ox, oy, oz, l0, theta0, phi0
    real(dp), intent(OUT) :: l, theta, phi, sx, sy, sz, ux, uy, uz, chi2
    procedure(calcF), pointer :: calcF_func => null()
    type(input_vars) :: inp
    type(updating_vars) :: upd
    integer :: i
    real(dp) :: dxdt, dxdp, dydt, dydp
    real(dp) :: d2rdt2, d2rdp2, d2rdtdp, d2zdt2, d2zdp2, d2zdtdp
    real(dp) :: sinPhi, cosPhi
    real(dp) :: dnxdt, dnydt, dnzdt, dnxdp, dnydp, dnzdp, dlndt, dlndp
    real(dp) :: duxdt, duydt, duzdt, duxdp, duydp, duzdp
    real(dp) ::  jac_xl,  jac_xt,  jac_xp, &
                 jac_yl,  jac_yt,  jac_yp, &
                 jac_zl,  jac_zt,  jac_zp
    real(dp) :: ijac_lx, ijac_ly, ijac_lz, &
                ijac_px, ijac_py, ijac_pz, &
                ijac_tx, ijac_ty, ijac_tz

    ! Initialize input variables in the input variable structure
    inp%ox = ox
    inp%oy = oy
    inp%oz = oz
    inp%surf = surf

    ! Initialize the updating solution variables
    upd%l     = l0
    upd%theta = theta0
    upd%phi   = phi0

    ! Calculate the residuals (and other updating values) for the initial guess
    calcF_func => calcF_surf_perp_intersect_3d
    call calcF_func(inp, upd)

    do i = 1, maxIter+1

        ! terminate if solution has sufficiently converged
        upd%chi2 = upd%fx**2 + upd%fy**2 + upd%fz**2
        if (upd%chi2 < dist_tol**2 .or. i == maxIter) then 
            l = upd%l;   theta = upd%theta; phi = upd%phi
            sx = upd%sx;    sy = upd%sy;     sz = upd%sz
            ux = upd%ux;    uy = upd%uy;     uz = upd%uz
            chi2 = upd%chi2
            exit
        end if

        ! for debugging
        !write(*,fmt='(A, I0, A, F8.4, A, F8.4, A, F8.4, A, ES10.3)') &
        !    '    Iter ', i, ':  l = ', upd%l, &
        !    '; theta = ', upd%theta*180./3.14159, &
        !    '; phi = ', upd%phi*180./3.14159, '; chi2 = ', upd%chi2

        ! Sines and cosines of the current toroidal angles
        sinPhi   = sin(upd%phi)
        cosPhi   = cos(upd%phi)

        ! second derivatives and cross-derivatives for the Jacobian
        d2rdt2  = surf_d2rdt2(inp%surf, upd%theta, upd%phi)
        d2zdt2  = surf_d2zdt2(inp%surf, upd%theta, upd%phi)
        d2rdp2  = surf_d2rdp2(inp%surf, upd%theta, upd%phi)
        d2zdp2  = surf_d2zdp2(inp%surf, upd%theta, upd%phi)
        d2rdtdp = surf_d2rdtdp(inp%surf, upd%theta, upd%phi)
        d2zdtdp = surf_d2zdtdp(inp%surf, upd%theta, upd%phi)

        ! derivatives of x and y coordinates on the vessel
        dxdt = upd%drdt * cosPhi
        dydt = upd%drdt * sinPhi
        dxdp = upd%drdp * cosPhi - upd%sr * sinPhi
        dydp = upd%drdp * sinPhi + upd%sr * cosPhi
        
        ! compute the Jacobian, first column (derivs of residuals wrt l)
        jac_xl = upd%ux
        jac_yl = upd%uy
        jac_zl = upd%uz

        ! compute the Jacobian, second column (derivs of residuals wrt theta)
        dlndt = ( (upd%dzdt*upd%drdp - upd%dzdt*upd%drdt)               &
                      *(  d2zdt2*upd%drdp + upd%dzdt*d2rdtdp    &
                        - d2zdtdp*upd%drdt - upd%dzdt*d2rdt2  ) &
                 + upd%sr*upd%drdt * (upd%dzdt**2 + upd%drdt**2)        &
                 + upd%sr**2 * (upd%dzdt*d2zdt2 + upd%drdt*d2rdt2)) &
                                   / upd%ln
        dnxdt = (d2zdt2*upd%drdp + upd%dzdt*d2rdtdp                &
                    - d2zdtdp*upd%drdt - upd%dzdt*d2rdt2) * sinPhi &
                + (d2zdt2*upd%sr + upd%dzdt*upd%drdt) * cosPhi
        dnydt = (-d2zdt2*upd%drdp - upd%dzdt*d2rdtdp               &
                    + d2zdtdp*upd%drdt + upd%dzdt*d2rdt2) * cosPhi &
                + (d2zdt2*upd%sr + upd%dzdt*upd%drdt) * sinPhi
        dnzdt = -d2rdt2*upd%sr - upd%drdt*upd%drdt
        duxdt = (upd%ln*dnxdt - upd%nx*dlndt) / upd%ln**2
        duydt = (upd%ln*dnydt - upd%ny*dlndt) / upd%ln**2
        duzdt = (upd%ln*dnzdt - upd%nz*dlndt) / upd%ln**2
        jac_xt = dxdt + upd%l*duxdt
        jac_yt = dydt + upd%l*duydt
        jac_zt = upd%dzdt + upd%l*duzdt

        ! compute the Jacobian, third column (derivs of residuals wrt phi)
        dlndp = ( (upd%dzdt*upd%drdp - upd%dzdp*upd%drdt)                   &
                      *(  d2zdtdp*upd%drdp + upd%dzdt*d2rdp2        &
                        - d2zdp2*upd%drdt - upd%dzdp*d2rdtdp  )     &
                 + upd%sr*upd%drdp * (upd%dzdt**2 + upd%drdt**2)            &
                 + upd%sr**2 * (upd%dzdt*d2zdtdp + upd%drdt*d2rdtdp)) &
                                / upd%ln
        dnxdp = (d2zdtdp*upd%drdp + upd%dzdt*d2rdp2                 &
                     - upd%dzdt*upd%sr + upd%dzdp*upd%drdt) * sinPhi        &
                 + (d2zdtdp*upd%sr + 2*upd%dzdt*upd%drdp                &
                     - d2zdp2*upd%drdt - upd%dzdp*d2rdtdp) * cosPhi
        dnydp = (d2zdtdp*upd%sr + 2*upd%dzdt*upd%drdp       &
                     - upd%dzdp*upd%drdt) * sinPhi              &
                 + (-d2zdtdp*upd%drdp - upd%dzdt*d2rdp2 &
                     + upd%dzdt*upd%sr + d2zdp2*upd%drdt    &
                     + upd%dzdp*d2rdtdp) * cosPhi
        dnzdp = -d2rdtdp*upd%sr - upd%drdt*upd%drdp
        duxdp = (upd%ln*dnxdp - upd%nx*dlndp) / upd%ln**2
        duydp = (upd%ln*dnydp - upd%ny*dlndp) / upd%ln**2
        duzdp = (upd%ln*dnzdp - upd%nz*dlndp) / upd%ln**2
        jac_xp = dxdp + upd%l*duxdp
        jac_yp = dydp + upd%l*duydp
        jac_zp = upd%dzdp + upd%l*duzdp

        ! invert the Jacobian
        call inverse_3x3(jac_xl, jac_xt, jac_xp, jac_yl, jac_yt, jac_yp, &
                         jac_zl, jac_zt, jac_zp, ijac_lx, ijac_ly, ijac_lz, &
                         ijac_tx, ijac_ty, ijac_tz, ijac_px, ijac_py, ijac_pz)

        ! The Newton step: -inv(Jacobian).residuals
        upd%sl = -ijac_lx * upd%fx - ijac_ly * upd%fy - ijac_lz * upd%fz
        upd%st = -ijac_tx * upd%fx - ijac_ty * upd%fy - ijac_tz * upd%fz
        upd%sp = -ijac_px * upd%fx - ijac_py * upd%fy - ijac_pz * upd%fz

        upd%grad_l = upd%fx*jac_xl + upd%fy*jac_yl + upd%fz*jac_zl
        upd%grad_t = upd%fx*jac_xt + upd%fy*jac_yt + upd%fz*jac_zt
        upd%grad_p = upd%fx*jac_xp + upd%fy*jac_yp + upd%fz*jac_zp
        upd%slope = upd%grad_l*upd%sl + upd%grad_t*upd%st + upd%grad_p*upd%sp 
        call updateStep(calcF_func, inp, upd)

    end do

end subroutine surf_perp_intersect_3d

subroutine calcF_surf_perp_intersect_3d(inp, upd)

    use surface_calc,        only: surface, surf_r, surf_z, &
                                   surf_drdt, surf_drdp, surf_dzdt, surf_dzdp
    use surface_solve_types, only: input_vars, updating_vars

    implicit none

    type(input_vars), intent(IN) :: inp
    type(updating_vars, intent(INOUT) :: upd
    real :: sinPhi, cosPhi

    sinPhi = sin(upd%phi)
    cosPhi = cos(upd%phi)

    upd%sr = surf_r(inp%surf, upd%theta, upd%phi)
    upd%sx = upd%sr * cosPhi
    upd%sy = upd%sr * sinPhi
    upd%sz = surf_z(inp%surf, upd%theta, upd%phi)

    ! Derivatives for computing the normal vector
    upd%drdt = surf_drdt(inp%surf, upd%theta, upd%phi)
    upd%dzdt = surf_dzdt(inp%surf, upd%theta, upd%phi)
    upd%drdp = surf_drdp(inp%surf, upd%theta, upd%phi)
    upd%dzdp = surf_dzdp(inp%surf, upd%theta, upd%phi)
    
    ! normal vectors (not unit vectors)
    upd%nx =   upd%dzdt*(upd%drdp*sinPhi + upd%sr*cosPhi) &
                   - upd%dzdp*upd%drdt*sinPhi
    upd%ny =  -upd%dzdt*(upd%drdp*cosPhi - upd%sr*sinPhi) &
                   + upd%dzdp*upd%drdt*cosPhi
    upd%nz =   upd%drdt*(upd%drdp*cosPhi - upd%sr*sinPhi)*sinPhi &
                   - upd%drdt*(upd%drdp*sinPhi + upd%sr*cosPhi)*cosPhi

    ! unit normal vectors
    upd%ln = sqrt( (upd%dzdt*upd%drdp - upd%dzdp*upd%drdt)**2 &
                   + upd%sr**2*(upd%dzdt**2 + upd%drdt**2) )
    upd%ux = upd%nx / upd%ln
    upd%uy = upd%ny / upd%ln
    upd%uz = upd%nz / upd%ln

    ! components that should ideally be zero
    upd%fx = upd%sx + upd%l*upd%ux - inp%ox
    upd%fy = upd%sy + upd%l*upd%uy - inp%oy
    upd%fz = upd%sz + upd%l*upd%uz - inp%oz

end subroutine calcF_surf_perp_intersect_3d

!-------------------------------------------------------------------------------
! updateStep(calcF_func, inp, upd)
!
! Adjusts the input solution and associated residuals of a 3d nonlinear problem 
! according to the input Newton step calculated by a solver implemented in the 
! calling subroutine. The length of the step is reduced, if necessary, using
! a line-search algorithm.
! 
! Input parameters:
!     procedure(calcF), 
!         pointer :: calcF_func   -> pointer to a subroutine that updates the
!                                    the residuals and possibly other updating
!                                    variables used by the solver
!     type(input_vars) :: inp     -> structure containing variables input to
!                                    the solver that are input to calcF_func
!
! Updated parameters:
!     type(updating_vars) :: upd  -> structure containing the solution values,
!                                    residuals, and other variables that are
!                                    updated by calcF_func
!
! Line-search algorithm adapted from: W. H. Press, S. A. Teukolsky, 
! W. T. Vetterling, and B. P. Flannery, Numerical Recipes in C, 2nd Edition,
! Cambridge University Press, 1992.
!-------------------------------------------------------------------------------
subroutine updateStep(calcF_func, inp, upd)

    use magpie_globals, only: dist_tol, maxIter

    implicit none

    procedure(calcF), pointer, intent(IN) :: calcF_func
    type(input_vars), intent(IN) :: inp
    type(updating_vars), intent(INOUT) :: upd
    real(dp) :: lambda, fsq, fsq_prev, l_prev, theta_prev, phi_prev, &
                newlambda, lambda2, rhs1, rhs2, a, b, disc, fsq2, test
    real(dp) :: alpha = 1.0e-4, tol = 1.0e-7
    integer :: i, j

    ! Adjust the step according to the line search algorithm if necessary
    lambda = 1.0
    j = 0
    fsq_prev = 0.5 * (upd%fx**2 + upd%fy**2 + upd%fz**2)
    l_prev = upd%l
    theta_prev = upd%theta
    phi_prev = upd%phi

    do i = 1, maxIter
        j = j + 1
        upd%l     =     l_prev + lambda * upd%sl
        upd%theta = theta_prev + lambda * upd%st
        upd%phi   =   phi_prev + lambda * upd%sp

        call calcF_func(inp, upd)

        fsq = 0.5 * (upd%fx**2 + upd%fy**2 + upd%fz**2)

        !if (.true.) exit
        if (fsq < fsq_prev + alpha*lambda*upd%slope) exit

        if (lambda*max(abs(upd%sl), abs(upd%st), abs(upd%sp)) < tol) exit

        if (j == 1) then
            newlambda = -upd%slope/(2 * (fsq - fsq_prev - upd%slope))
        else
            rhs1 = fsq - fsq_prev - lambda*upd%slope
            rhs2 = fsq2 - fsq_prev - lambda2*upd%slope
            a = (rhs1/lambda**2 - rhs2/lambda2**2)/(lambda - lambda2)
            b = (-lambda2*rhs1/lambda**2 + lambda*rhs2/lambda2**2) &
                    /(lambda - lambda2)
            if (a == 0.0) then
                newlambda = -upd%slope/(2.0 * b) 
            else
                disc = b**2 - 3.0*a*upd%slope
                if (disc < 0.0) then
                    stop 'vessel_perp_intersect_3d: roundoff problem'
                else
                    newlambda = (-b + sqrt(disc))/(3.0*a)
                end if
            end if
            if (newlambda > 0.5 * lambda) newlambda = 0.5 * lambda
        end if

        lambda2 = lambda
        lambda = max(0.1*lambda, newlambda)
        fsq2 = fsq
        j = j + 1
    end do

end subroutine updateStep

!-------------------------------------------------------------------------------
! surf_perp_intersect_2d(phi, or, oz, l0, theta0, l, theta, r, z)
!
! Computes the poloidal angle at which a line extending from a given point
! (or, oz) in a poloidal plane interscts the cross-section of a toroidal
! surface at a particular toroidal angle. Also computes the distance from the
! input point to the point of intersection.
!
! Input parameters:
!     type(surface) :: surf   -> structure encapsulating the toroidal surface
!     real(dp) :: phi         -> toroidal angle (radians) of the cross-section
!     real(dp) :: or, oz      -> r and z coordinates (meters) of the ref. point
!     real(dp) :: l0, theta0  -> initial guesses of the distance l, as well as
!                               the theta coordinate of the intersection point
!                               on the surfsel
!
! Return parameters:
!     real(dp) :: l           -> distance from the point to the surfsel
!     real(dp) :: theta       -> theta coordinate of the intersection point
!     real(dp) :: r, z        -> r and z coordinates of the intersection point
!-------------------------------------------------------------------------------
subroutine surf_perp_intersect_2d(surf, phi, or, oz, l0, theta0, &
                                  l, theta, r, z, chi2)

    use magpie_globals, only: dist_tol, maxIter
    use algebra,        only: inverse_2x2
    use surface_calc,   only: surface, surf_r, surf_z, surf_drdt, surf_dzdt, &
                              surf_d2rdt2, surf_d2zdt2

    implicit none

    type(surface), intent(IN) :: surf
    real(dp), intent(IN)  :: phi, or, oz, l0, theta0
    real(dp), intent(OUT) :: l, theta, r, z, chi2
    type(input_vars) :: inp
    type(updating_vars) :: upd
    integer :: i
    real(dp) :: drdt, dzdt, d2rdt2, d2zdt2, norm, dnorm_dt, dnr_dt, dnz_dt
    real(dp) ::  jac_rl,  jac_rt,  jac_zl,  jac_zt
    real(dp) :: ijac_lr, ijac_lz, ijac_tr, ijac_tz

    ! Initialize relevant parameters in the input variable structure
    inp%or = or
    inp%oz = oz
    inp%phi = phi
    inp%surf = surf

    ! Initialize values for which to solve
    upd%l     = l0
    upd%theta = theta0

    do i = 1, maxIter

        ! Compute surface coordinates and derivatives
        upd%sr = surf_r(inp%surf, upd%theta, inp%phi)
        upd%sz = surf_z(inp%surf, upd%theta, inp%phi)
        drdt = surf_drdt(inp%surf, upd%theta, inp%phi)
        dzdt = surf_dzdt(inp%surf, upd%theta, inp%phi)
        d2rdt2 = surf_d2rdt2(inp%surf, upd%theta, inp%phi)
        d2zdt2 = surf_d2zdt2(inp%surf, upd%theta, inp%phi)

        ! Unit vector normal to the surface
        norm = sqrt(drdt**2 + dzdt**2)
        upd%nr = dzdt/norm
        upd%nz = -drdt/norm

        ! components that should ideally be zero
        upd%fr = upd%sr + upd%l*upd%nr - inp%or
        upd%fz = upd%sz + upd%l*upd%nz - inp%oz

        ! terminate if solution has sufficiently converged
        upd%chi2 = upd%fr*upd%fr + upd%fz*upd%fz
        if (upd%chi2 < dist_tol*dist_tol .or. i == maxIter) then
            l     = upd%l
            theta = upd%theta
            chi2  = upd%chi2
            r     = upd%sr
            z     = upd%sz
            !write(*,fmt='(A,I2,A,F5.2,A,F7.2,A,F5.2,A,F7.2,A,ES8.1)') &
            !    '    iter ', i, ': phi=', inp%phi*180./3.14159265358, &
            !    ', theta0=', theta0*180./3.14159265358, ', l=', l, &
            !    ', theta=', theta*180./3.14159265358, ', chi2=', chi2
            exit
        end if

        ! Compute the Jacobian, first column
        jac_rl = upd%nr
        jac_zl = upd%nz
        
        ! compute the Jacobian, second column
        dnorm_dt = (drdt*d2rdt2 + dzdt*d2zdt2) / norm
        dnr_dt =  d2zdt2/norm - dzdt*dnorm_dt/norm**2
        dnz_dt = -d2rdt2/norm + drdt*dnorm_dt/norm**2
        jac_rt = drdt + upd%l*dnr_dt
        jac_zt = dzdt + upd%l*dnz_dt

        ! invert the Jacobian
        call inverse_2x2( jac_rl,  jac_rt,  jac_zl,  jac_zt, &
                         ijac_lr, ijac_lz, ijac_tr, ijac_tz )

        ! The Newton step
        upd%sl = -ijac_lr * upd%fr - ijac_lz * upd%fz
        upd%st = -ijac_tr * upd%fr - ijac_tz * upd%fz

        ! Update the parameters
        upd%l     = upd%l     + upd%sl
        upd%theta = upd%theta + upd%st

    end do

end subroutine surf_perp_intersect_2d

end module surface_solve

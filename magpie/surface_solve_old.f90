!-------------------------------------------------------------------------------
! surface_solve.f90
!
! Contains a few Newton-Raphson solvers for various systems of nonlinear 
! equations useful for addressing geometric questions about toroidal surfaces.
!
! Author:  K. C. Hammond
! Contact: khammod@pppl.gov
! Updated: 2019-12-21
!-------------------------------------------------------------------------------
module surface_solve

use magpie_globals, only: dp
use surface_calc,   only: surface

implicit none

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
    use surface_calc,   only: surface, surf_drdt, surf_dzdt, &
                              surf_drdp, surf_dzdp

    implicit none

    type(surface), intent(IN) :: surf
    real(dp), intent(IN)  :: ox, oy, oz, ax, ay, az, l0, theta0, phi0
    real(dp), intent(OUT) :: l, theta, phi, sx, sy, sz, chi2
    INTEGER :: i
    real(dp) :: sr, fx, fy, fz, sl, st, sp, grad_l, grad_t, grad_p, slope
    real(dp) :: sinPhi, cosPhi
    real(dp) :: drdt, drdp, dzdt, dzdp
    real(dp) ::  jac_xl,  jac_xt,  jac_xp, &
                 jac_yl,  jac_yt,  jac_yp, &
                 jac_zl,  jac_zt,  jac_zp
    real(dp) :: ijac_lx, ijac_ly, ijac_lz, &
                ijac_px, ijac_py, ijac_pz, &
                ijac_tx, ijac_ty, ijac_tz

    ! constant components of the Jacobian (first column)
    jac_xl = -ax
    jac_yl = -ay
    jac_zl = -az

    ! Initialize values for which to solve
    l     = l0
    theta = theta0
    phi   = phi0

    ! Calculate fx, fy, fz (should be zero if a solution is found)
    ! as well as coordinates on the vessel at (theta, phi)
    call calcF_surf_dist_3d(surf, l, theta, phi, ox, oy, oz, ax, ay, az, &
                            fx, fy, fz, sr, sx, sy, sz)

    do i = 1, maxIter

        ! terminate if solution has sufficiently converged
        chi2 = fx*fx + fy*fy + fz*fz

        ! for debugging
        write(*,fmt='(A, I0, A, F8.4, A, F8.4, A, F8.4, A, ES10.3)') &
            '    Iter ', i, ':  l = ', l, '; theta = ', theta*180./3.14159, &
            '; phi = ', phi*180./3.14159, '; chi2 = ', chi2

        if (chi2 < dist_tol*dist_tol) exit

        ! Update sinusoidal terms
        sinPhi = sin(phi)
        cosPhi = cos(phi)
        
        ! compute the Jacobian, second column
        drdt = surf_drdt(surf, theta, phi)
        dzdt = surf_dzdt(surf, theta, phi)
        jac_xt = drdt * cosPhi
        jac_yt = drdt * sinPhi
        jac_zt = dzdt

        ! compute the Jacobian, third column
        drdp = surf_drdp(surf, theta, phi)
        dzdp = surf_dzdp(surf, theta, phi)
        jac_xp = drdp * cosPhi - sr * sinPhi
        jac_yp = drdp * sinPhi + sr * cosPhi
        jac_zp = dzdp

        ! invert the Jacobian
        call inverse_3x3(jac_xl, jac_xt, jac_xp, jac_yl, jac_yt, jac_yp, &
                         jac_zl, jac_zt, jac_zp, ijac_lx, ijac_ly, ijac_lz, &
                         ijac_tx, ijac_ty, ijac_tz, ijac_px, ijac_py, ijac_pz)

        ! The Newton step
        sl = -ijac_lx * fx - ijac_ly * fy - ijac_lz * fz
        st = -ijac_tx * fx - ijac_ty * fy - ijac_tz * fz
        sp = -ijac_px * fx - ijac_py * fy - ijac_pz * fz

        ! The gradient of chi2 in the Newton direction
        grad_l = fx*jac_xl + fy*jac_yl + fz*jac_zl
        grad_t = fx*jac_xt + fy*jac_yt + fz*jac_zt
        grad_p = fx*jac_xp + fy*jac_yp + fz*jac_zp
        slope = grad_l*sl + grad_t*st + grad_p*sp 

        ! Update l, theta, and phi, refining the Newton step if necessary
        call updateStep_surf_dist_3d(sl, st, sp, l, theta, phi, slope, &
                 surf, ox, oy, oz, ax, ay, az, fx, fy, fz, sr, sx, sy, sz)

    end do

end subroutine surf_dist_3d

subroutine calcF_surf_dist_3d(surf, l, theta, phi, ox, oy, oz, ax, ay, az, &
                              fx, fy, fz, sr, sx, sy, sz)

    use surface_calc, only: surface, surf_r, surf_z

    implicit none

    type(surface), intent(IN) :: surf
    real(dp), intent(IN)  :: l, theta, phi, ox, oy, oz, ax, ay, az
    real(dp), intent(OUT) :: fx, fy, fz, sr, sx, sy, sz

    ! r, x, y, and z coordinates on the vessel corresponding to theta and phi
    sr = surf_r(surf, theta, phi)
    sx = sr * cos(phi)
    sy = sr * sin(phi)
    sz = surf_z(surf, theta, phi)

    ! Components that should ideally be zero
    fx = sx - (ox + l*ax)
    fy = sy - (oy + l*ay)
    fz = sz - (oz + l*az)

end subroutine calcF_surf_dist_3d

subroutine updateStep_surf_dist_3d(sl, st, sp, l, theta, phi, &
               slope, surf, ox, oy, oz, ax, ay, az, fx, fy, fz, sr, sx, sy, sz)

    use magpie_globals, only: dist_tol, maxIter
    use surface_calc,   only: surface

    implicit none

    real(dp), intent(IN)  :: sl, st, sp, slope, ox, oy, oz, ax, ay, az
    type(surface), intent(IN) :: surf
    real(dp), intent(OUT) :: l, theta, phi, fx, fy, fz, sr, sx, sy, sz
    real(dp) :: lambda, fsq, fsq_prev, l_prev, theta_prev, phi_prev, &
                newlambda, lambda2, rhs1, rhs2, a, b, disc, fsq2, test
    real(dp) :: alpha = 1.0e-4, tol = 1.0e-7
    INTEGER :: i, j

    ! Adjust the step according to the line search algorithm if necessary
    lambda = 1.0
    j = 0
    fsq_prev = 0.5 * (fx**2 + fy**2 + fz**2)
    l_prev = l
    theta_prev = theta
    phi_prev = phi

    do i = 1, maxIter
        j = j + 1
        l     =     l_prev + lambda * sl
        theta = theta_prev + lambda * st
        phi   =   phi_prev + lambda * sp

        call calcF_surf_dist_3d(surf, l, theta, phi, ox, oy, oz, ax, ay, az, &
                                fx, fy, fz, sr, sx, sy, sz)

        fsq = 0.5 * (fx**2 + fy**2 + fz**2)

        !if (.true.) exit
        if (fsq < fsq_prev + alpha*lambda*slope) exit

        if (lambda*max(abs(sl), abs(st), abs(sp)) < tol) exit

        if (j == 1) then
            newlambda = -slope/(2 * (fsq - fsq_prev - slope))
        else
            rhs1 = fsq - fsq_prev - lambda*slope
            rhs2 = fsq2 - fsq_prev - lambda2*slope
            a = (rhs1/lambda**2 - rhs2/lambda2**2)/(lambda - lambda2)
            b = (-lambda2*rhs1/lambda**2 + lambda*rhs2/lambda2**2) &
                    /(lambda - lambda2)
            if (a == 0.0) then
                newlambda = -slope/(2.0 * b) 
            else
                disc = b**2 - 3.0*a*slope
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

end subroutine updateStep_surf_dist_3d

end module surface_solve

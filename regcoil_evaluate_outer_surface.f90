subroutine regcoil_evaluate_outer_surface()

  use regcoil_compute_outer_surface_mod
  use regcoil_variables
  
  use omp_lib
  
  implicit none
  
  real(dp), dimension(:,:), allocatable :: major_R_outer, z_outer
  real(dp) :: angle, sinangle, cosangle
  integer :: j, itheta, izeta
  real(dp) :: x_new, y_new, z_new, factor, factor2
  integer :: tic, toc, countrate, tic1, toc1

  call system_clock(tic,countrate)
  if (verbose) print "(a)"," Evaluating the outer magnetization surface in terms of the standard toroidal angle."

  ! Make sure we have the Fourier representation of the latest d
  d0 = sum(d)/(ntheta_coil*nzeta_coil)
  dmnc=0
  dmns=0
  factor = (2.0d+0) / (ntheta_coil * nzeta_coil)
  do izeta = 1, nzeta_coil
     do itheta = 1, ntheta_coil
        do j = 1, mnmax_magnetization
           angle = xm_magnetization(j) * theta_coil(itheta) - xn_magnetization(j) * zeta_coil(izeta)
           sinangle = sin(angle)
           cosangle = cos(angle)
           factor2 = factor
           ! The next 2 lines ensure inverse Fourier transform(Fourier transform) = identity
           if (mod(ntheta_coil,2) == 0 .and.     xm_magnetization(j)  ==    (ntheta_coil/2)) factor2 = factor2 / 2
           if (mod( nzeta_coil,2) == 0 .and. abs(xn_magnetization(j)) == nfp*(nzeta_coil/2)) factor2 = factor2 / 2
           dmnc(j) = dmnc(j) + d(itheta, izeta) * cosangle * factor2
           dmns(j) = dmns(j) + d(itheta, izeta) * sinangle * factor2
        end do
     end do
  end do

  allocate(major_R_outer(ntheta_coil, nzeta_coil))
  allocate(z_outer(ntheta_coil, nzeta_coil))
  major_R_outer = 0
  z_outer = 0

  call system_clock(tic1)
  !$OMP PARALLEL DEFAULT(NONE), PRIVATE(x_new,y_new,z_new,itheta), SHARED(ntheta_coil,nzeta_coil,theta_coil,zetal_coil,major_R_outer,z_outer)

  !$OMP DO
  do izeta = 1,nzeta_coil
     do itheta = 1,ntheta_coil           
        ! Compute r:
        call regcoil_compute_outer_surface_xyz_of_thetazeta(theta_coil(itheta), &
             zetal_coil(izeta),x_new,y_new,z_new)
        major_R_outer(itheta,izeta) = sqrt(x_new * x_new + y_new * y_new)
        z_outer(itheta,izeta) = z_new
     end do
  end do
  !$OMP END DO

  !$OMP END PARALLEL

  call system_clock(toc1)
  if (verbose) print *,"  Computing offset points:",real(toc1-tic1)/countrate,"sec"

  ! Fourier transform the result.
  ! This is not a rate-limiting step, so for clarity of code, we don't bother with an FFT.
  call system_clock(tic1)
  allocate(rmnc_outer(mnmax_coil))
  allocate(rmns_outer(mnmax_coil))
  allocate(zmnc_outer(mnmax_coil))
  allocate(zmns_outer(mnmax_coil))
  rmnc_outer = 0
  rmns_outer = 0
  zmnc_outer = 0
  zmns_outer = 0
  factor = (2.0d+0) / (ntheta_coil * nzeta_coil)
  do izeta = 1, nzeta_coil
     do itheta = 1, ntheta_coil
        do j = 2, mnmax_coil
           angle = xm_coil(j) * theta_coil(itheta) - xn_coil(j) * zeta_coil(izeta)
           sinangle = sin(angle)
           cosangle = cos(angle)
           factor2 = factor
           ! The next 2 lines ensure inverse Fourier transform(Fourier transform) = identity
           if (mod(ntheta_coil,2) == 0 .and.     xm_coil(j)  ==    (ntheta_coil/2)) factor2 = factor2 / 2
           if (mod( nzeta_coil,2) == 0 .and. abs(xn_coil(j)) == nfp*(nzeta_coil/2)) factor2 = factor2 / 2
           rmnc_outer(j) = rmnc_outer(j) + major_R_outer(itheta, izeta) * cosangle * factor2
           rmns_outer(j) = rmns_outer(j) + major_R_outer(itheta, izeta) * sinangle * factor2
           zmnc_outer(j) = zmnc_outer(j) + z_outer(itheta, izeta) * cosangle * factor2
           zmns_outer(j) = zmns_outer(j) + z_outer(itheta, izeta) * sinangle * factor2
        end do
     end do
  end do
  rmnc_outer(1) = sum(major_R_outer) / (ntheta_coil * nzeta_coil)
  zmnc_outer(1) = sum(z_outer) / (ntheta_coil * nzeta_coil)
  
  if (.not. lasym) then
     rmns_outer = 0
     zmnc_outer = 0
  end if

  deallocate(major_R_outer, z_outer)

  call system_clock(toc)
  if (verbose) print *,"Done. Took ",real(toc-tic)/countrate,"sec"

end subroutine regcoil_evaluate_outer_surface


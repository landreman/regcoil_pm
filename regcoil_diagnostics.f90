subroutine regcoil_diagnostics(ilambda)

  use regcoil_variables
  use stel_constants
  use stel_kinds

  implicit none

  integer, intent(in) :: ilambda
  integer :: tic, toc, countrate
  integer :: itheta, izeta, js, offset
  real(dp), dimension(:,:), allocatable :: temp_3D
  real(dp) :: chi2_M_alt

  call system_clock(tic,countrate)

  ! Unpack the solution vector.
  abs_M(:,:,:,ilambda) = 0
  do js = 1, ns_magnetization
     offset = 0*ns_magnetization*num_basis_functions + (js-1) * num_basis_functions
     M_R_mn(:,js,ilambda) = solution(offset+1:offset+num_basis_functions)
     M_R(:,:,js,ilambda) = reshape(matmul(basis_functions, solution(offset+1:offset+num_basis_functions)), (/ ntheta_coil, nzeta_coil /))

     offset = 1*ns_magnetization*num_basis_functions + (js-1) * num_basis_functions
     M_zeta_mn(:,js,ilambda) = solution(offset+1:offset+num_basis_functions)
     M_zeta(:,:,js,ilambda) = reshape(matmul(basis_functions, solution(offset+1:offset+num_basis_functions)), (/ ntheta_coil, nzeta_coil /))

     offset = 2*ns_magnetization*num_basis_functions + (js-1) * num_basis_functions
     M_Z_mn(:,js,ilambda) = solution(offset+1:offset+num_basis_functions)
     M_Z(:,:,js,ilambda) = reshape(matmul(basis_functions, solution(offset+1:offset+num_basis_functions)), (/ ntheta_coil, nzeta_coil /))

     abs_M(:,:,js,ilambda) = sqrt(M_R(:,:,js,ilambda) * M_R(:,:,js,ilambda) + M_zeta(:,:,js,ilambda) * M_zeta(:,:,js,ilambda) + M_Z(:,:,js,ilambda) * M_Z(:,:,js,ilambda))
  end do

  KDifference_x = d_x - matmul(f_x, solution)
  KDifference_y = d_y - matmul(f_y, solution)
  KDifference_z = d_z - matmul(f_z, solution)
  KDifference_Laplace_Beltrami = d_Laplace_Beltrami - matmul(f_Laplace_Beltrami, solution)
  this_K2_times_N = reshape(KDifference_x*KDifference_x + KDifference_y*KDifference_y + KDifference_z*KDifference_z, (/ ntheta_coil, nzeta_coil /)) &
       / norm_normal_coil
  chi2_K(ilambda) = nfp * dtheta_coil * dzeta_coil * sum(this_K2_times_N)
  K2(:,:,ilambda) = this_K2_times_N / norm_normal_coil

  chi2_M(ilambda) = dot_product(solution, matmul(matrix_regularization, solution)) * nfp
  ! Compute chi2_M a second way, as a sanity test:
  allocate(temp_3D(ntheta_coil, nzeta_coil, ns_magnetization))
  temp_3D = abs_M(:,:,:,ilambda) * abs_M(:,:,:,ilambda) * Jacobian_coil
  chi2_M_alt = 0
  do js = 1, ns_magnetization
     ! I NEED TO COMPLETE THIS NEXT BIT.
     chi2_M_alt = chi2_M_alt + sum(temp_3D(:,:,js))
  end do
  if (verbose) print "(3(a,es22.14))","2 methods of computing chi2_M that should agree: method 1 =",chi2_M(ilambda),", method 2 =",chi2_M_alt,", difference =",chi2_M(ilambda) - chi2_M_alt

  Bnormal_total(:,:,ilambda) = (reshape(matmul(g,solution),(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) &
       + Bnormal_from_plasma_current + Bnormal_from_net_coil_currents
  
  max_Bnormal(ilambda) = maxval(abs(Bnormal_total(:,:,ilambda)))
  max_M(ilambda) = maxval(abs_M(:,:,:,ilambda))
  
  chi2_B(ilambda) = nfp * dtheta_plasma * dzeta_plasma &
       * sum(Bnormal_total(:,:,ilambda) * Bnormal_total(:,:,ilambda) * norm_normal_plasma)
  
  call system_clock(toc)
  if (verbose) print *,"  Diagnostics: ",real(toc-tic)/countrate," sec."
  if (verbose) print "(a,es10.3,a,es10.3)","   chi2_B:",chi2_B(ilambda),",  chi2_M:",chi2_M(ilambda)
  if (verbose) print "(a,es10.3,a,es10.3,a,es10.3)","   max(B_n):",max_Bnormal(ilambda),",  max(M):",max_M(ilambda),",  rms M:",sqrt(chi2_M(ilambda)/area_coil)


end subroutine regcoil_diagnostics

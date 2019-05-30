subroutine regcoil_diagnostics(ilambda)

  use regcoil_variables
  use stel_constants
  use stel_kinds

  implicit none

  integer, intent(in) :: ilambda
  integer :: tic, toc, countrate
  integer :: itheta, izeta, js, offset
  real(dp) :: chi2_M_alt, chi2_B_alt
  real(dp), dimension(:), allocatable :: temp_array

  call system_clock(tic,countrate)

  ! Unpack the solution vector.
  do js = 1, ns_magnetization
     offset = 0*ns_magnetization*num_basis_functions + (js-1) * num_basis_functions
     magnetization_vector_mn(:,js,1,ilambda) = solution(offset+1:offset+num_basis_functions)
     magnetization_vector( :,:,js,1,ilambda) = reshape(matmul(basis_functions, solution(offset+1:offset+num_basis_functions)), (/ ntheta_coil, nzeta_coil /))

     offset = 1*ns_magnetization*num_basis_functions + (js-1) * num_basis_functions
     magnetization_vector_mn(:,js,2,ilambda) = solution(offset+1:offset+num_basis_functions)
     magnetization_vector( :,:,js,2,ilambda) = reshape(matmul(basis_functions, solution(offset+1:offset+num_basis_functions)), (/ ntheta_coil, nzeta_coil /))

     offset = 2*ns_magnetization*num_basis_functions + (js-1) * num_basis_functions
     magnetization_vector_mn(:,js,3,ilambda) = solution(offset+1:offset+num_basis_functions)
     magnetization_vector( :,:,js,3,ilambda) = reshape(matmul(basis_functions, solution(offset+1:offset+num_basis_functions)), (/ ntheta_coil, nzeta_coil /))

     abs_M(:,:,js,ilambda) = sqrt(magnetization_vector(:,:,js,1,ilambda) * magnetization_vector(:,:,js,1,ilambda) &
          + magnetization_vector(:,:,js,2,ilambda) * magnetization_vector(:,:,js,2,ilambda) &
          + magnetization_vector(:,:,js,3,ilambda) * magnetization_vector(:,:,js,3,ilambda))
  end do

  chi2_M(ilambda) = dot_product(solution, matmul(matrix_regularization, solution)) * nfp
  ! Compute chi2_M a second way, as a sanity test. This second method is probably slower, so it could eventually be commented out.
  chi2_M_alt = 0
  do itheta = 1, ntheta_coil
     do izeta = 1, nzeta_coil
        chi2_M_alt = chi2_M_alt + d(itheta,izeta) * dot_product(s_weights * Jacobian_coil(itheta,izeta,:), &
             matmul(interpolate_magnetization_to_integration, magnetization_vector(itheta,izeta,:,1,ilambda)) ** 2 &
             + matmul(interpolate_magnetization_to_integration, magnetization_vector(itheta,izeta,:,2,ilambda)) ** 2 &
             + matmul(interpolate_magnetization_to_integration, magnetization_vector(itheta,izeta,:,3,ilambda)) ** 2)
     end do
  end do
  chi2_M_alt = chi2_M_alt * nfp * dtheta_coil * dzeta_coil
  if (verbose) print "(3(a,es22.14))","   2 methods of computing chi2_M that should agree: method 1 =",chi2_M(ilambda),", method 2 =",chi2_M_alt,", relative difference =",(chi2_M(ilambda) - chi2_M_alt) / (abs(chi2_M(ilambda)) + abs(chi2_M_alt))


  allocate(temp_array(ntheta_plasma * nzeta_plasma))
  temp_array = 0
  do js = 1, ns_magnetization
     temp_array = temp_array + matmul(g(:,:,js,1), magnetization_vector_mn(:,js,1,ilambda)) + matmul(g(:,:,js,2), magnetization_vector_mn(:,js,2,ilambda)) + matmul(g(:,:,js,3), magnetization_vector_mn(:,js,3,ilambda))
  end do
  !Bnormal_total(:,:,ilambda) = (reshape(matmul(g,solution),(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) &
  Bnormal_total(:,:,ilambda) = (reshape(temp_array,(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) &
       + Bnormal_from_TF_and_plasma_current
  deallocate(temp_array)

  max_Bnormal(ilambda) = maxval(abs(Bnormal_total(:,:,ilambda)))
  max_M(ilambda) = maxval(abs_M(:,:,:,ilambda))
  min_M(ilambda) = minval(abs_M(:,:,:,ilambda))
  
  chi2_B(ilambda) = nfp * dtheta_plasma * dzeta_plasma &
       * sum(Bnormal_total(:,:,ilambda) * Bnormal_total(:,:,ilambda) * norm_normal_plasma)
  ! Compute chi2_B a second way, as a sanity test.
  chi2_B_alt = nfp * dot_product(solution, matmul(matrix_B, solution)) - 2 * nfp * dot_product(RHS_B, solution) &
       + nfp * dtheta_plasma * dzeta_plasma * sum(norm_normal_plasma * Bnormal_from_TF_and_plasma_current * Bnormal_from_TF_and_plasma_current)
  if (verbose) print "(3(a,es22.14))","   2 methods of computing chi2_B that should agree: method 1 =",chi2_B(ilambda),", method 2 =",chi2_B_alt,", relative difference =",(chi2_B(ilambda) - chi2_B_alt) / (abs(chi2_B(ilambda)) + abs(chi2_B_alt))
  
  call system_clock(toc)
  if (verbose) print *,"  Diagnostics: ",real(toc-tic)/countrate," sec."
  if (verbose) print "(2(a,es10.3))","   chi2_B:",chi2_B(ilambda),",  chi2_M:",chi2_M(ilambda)
  if (verbose) print "(3(a,es10.3))","   max(B_n):",max_Bnormal(ilambda),",  max(M):",max_M(ilambda),",  min(M):",min_M(ilambda)


end subroutine regcoil_diagnostics

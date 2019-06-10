subroutine regcoil_diagnostics(isaved)

  use regcoil_variables
  use stel_constants
  use stel_kinds

  implicit none

  integer, intent(in) :: isaved
  integer :: tic, toc, countrate
  integer :: itheta, izeta, js, offset
  real(dp) :: chi2_M_alt, chi2_B_alt
  real(dp), dimension(:), allocatable :: temp_array

  call system_clock(tic,countrate)

  ! Unpack the solution vector.
  do js = 1, ns_magnetization
     offset = 0*ns_magnetization*num_basis_functions + (js-1) * num_basis_functions
     magnetization_vector_mn(:,js,1,isaved) = solution(offset+1:offset+num_basis_functions)
     magnetization_vector( :,:,js,1,isaved) = reshape(matmul(basis_functions_R, solution(offset+1:offset+num_basis_functions)), (/ ntheta_coil, nzeta_coil /))

     offset = 1*ns_magnetization*num_basis_functions + (js-1) * num_basis_functions
     magnetization_vector_mn(:,js,2,isaved) = solution(offset+1:offset+num_basis_functions)
     magnetization_vector( :,:,js,2,isaved) = reshape(matmul(basis_functions_zeta_Z, solution(offset+1:offset+num_basis_functions)), (/ ntheta_coil, nzeta_coil /))

     offset = 2*ns_magnetization*num_basis_functions + (js-1) * num_basis_functions
     magnetization_vector_mn(:,js,3,isaved) = solution(offset+1:offset+num_basis_functions)
     magnetization_vector( :,:,js,3,isaved) = reshape(matmul(basis_functions_zeta_Z, solution(offset+1:offset+num_basis_functions)), (/ ntheta_coil, nzeta_coil /))

     abs_M(:,:,js,isaved) = sqrt(magnetization_vector(:,:,js,1,isaved) * magnetization_vector(:,:,js,1,isaved) &
          + magnetization_vector(:,:,js,2,isaved) * magnetization_vector(:,:,js,2,isaved) &
          + magnetization_vector(:,:,js,3,isaved) * magnetization_vector(:,:,js,3,isaved))
  end do

  d_iterations(:,:,isaved) = d

  do itheta = 1, ntheta_coil
     do izeta = 1, nzeta_coil
        s_averaged_abs_M(itheta,izeta,isaved) = dot_product(s_magnetization_weights, abs_M(itheta,izeta,:,isaved))
     end do
  end do

  chi2_M(isaved) = dot_product(solution, matmul(matrix_regularization, solution)) * nfp
  ! Compute chi2_M a second way, as a sanity test. This second method is probably slower, so it could eventually be commented out.
  chi2_M_alt = 0
  do itheta = 1, ntheta_coil
     do izeta = 1, nzeta_coil
        chi2_M_alt = chi2_M_alt + (d(itheta,izeta) ** regularization_d_exponent) * dot_product(s_weights * Jacobian_coil(itheta,izeta,:), &
             matmul(interpolate_magnetization_to_integration, magnetization_vector(itheta,izeta,:,1,isaved)) ** 2 &
             + matmul(interpolate_magnetization_to_integration, magnetization_vector(itheta,izeta,:,2,isaved)) ** 2 &
             + matmul(interpolate_magnetization_to_integration, magnetization_vector(itheta,izeta,:,3,isaved)) ** 2)
     end do
  end do
  chi2_M_alt = chi2_M_alt * nfp * dtheta_coil * dzeta_coil
  if (verbose) print "(3(a,es22.14))","   2 methods of computing chi2_M that should agree: method 1 =",chi2_M(isaved),", method 2 =",chi2_M_alt,", relative difference =",(chi2_M(isaved) - chi2_M_alt) / (abs(chi2_M(isaved)) + abs(chi2_M_alt))


  allocate(temp_array(ntheta_plasma * nzeta_plasma))
  temp_array = 0
  do js = 1, ns_magnetization
     temp_array = temp_array + matmul(g(:,:,js,1), magnetization_vector_mn(:,js,1,isaved)) + matmul(g(:,:,js,2), magnetization_vector_mn(:,js,2,isaved)) + matmul(g(:,:,js,3), magnetization_vector_mn(:,js,3,isaved))
  end do
  !Bnormal_total(:,:,isaved) = (reshape(matmul(g,solution),(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) &
  Bnormal_total(:,:,isaved) = (reshape(temp_array,(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) &
       + Bnormal_from_TF_and_plasma_current
  deallocate(temp_array)

  max_Bnormal(isaved) = maxval(abs(Bnormal_total(:,:,isaved)))
  max_M(isaved) = maxval(abs_M(:,:,:,isaved))
  min_M(isaved) = minval(abs_M(:,:,:,isaved))
  
  chi2_B(isaved) = nfp * dtheta_plasma * dzeta_plasma &
       * sum(Bnormal_total(:,:,isaved) * Bnormal_total(:,:,isaved) * norm_normal_plasma)
  ! Compute chi2_B a second way, as a sanity test.
  chi2_B_alt = nfp * dot_product(solution, matmul(matrix_B, solution)) - 2 * nfp * dot_product(RHS_B, solution) &
       + nfp * dtheta_plasma * dzeta_plasma * sum(norm_normal_plasma * Bnormal_from_TF_and_plasma_current * Bnormal_from_TF_and_plasma_current)
  if (verbose) print "(3(a,es22.14))","   2 methods of computing chi2_B that should agree: method 1 =",chi2_B(isaved),", method 2 =",chi2_B_alt,", relative difference =",(chi2_B(isaved) - chi2_B_alt) / (abs(chi2_B(isaved)) + abs(chi2_B_alt))
  
  call system_clock(toc)
  if (verbose) print *,"  Diagnostics: ",real(toc-tic)/countrate," sec."
  if (verbose) print "(2(a,es10.3))","   chi2_B:",chi2_B(isaved),",  chi2_M:",chi2_M(isaved)
  if (verbose) print "(3(a,es10.3))","   max(B_n):",max_Bnormal(isaved),",  max(M):",max_M(isaved),",  min(M):",min_M(isaved)


end subroutine regcoil_diagnostics

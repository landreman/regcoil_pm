subroutine lse_sensitivity()

  use global_variables
  use omp_lib

  implicit none

  integer :: ilambda, iomega, iflag, tic, toc, countrate, i_basis_function, izeta,itheta, tic_begin
  real(dp) :: sum_exp, dsum_expdomega, max_KN, dsum_expdphi, dRMS_Kdphi, q_tilde_denom, q_tilde_num, max_val, arg, dintKpdomega, dintKpdPhi
  real(dp), dimension(:), allocatable :: solution, RHS
  real(dp), dimension(:,:), allocatable :: dnorm_Kdomega, norm_K
  real(dp), dimension(:,:), allocatable :: dnorm_Kdphi
  real(dp), dimension(:,:), allocatable :: matrix
  real(dp), dimension(:), allocatable :: KDifference_x, KDifference_y, KDifference_z
  real(dp), dimension(:,:), allocatable :: dKDifferencedomega

  ! Variables needed by LAPACK:
  integer :: INFO, LWORK
  real(dp), dimension(:), allocatable :: WORK
  integer, dimension(:), allocatable :: IPIV

  allocate(WORK(1), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(IPIV(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dRMSKdomega(nomega_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(darea_coildomega(nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dtarget_optiondOmega(nomega_coil,nlambda))
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dtarget_optiondPhi(num_basis_functions,nlambda))
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnorm_Kdomega(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(norm_K(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(solution(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnorm_Kdphi(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(matrix(num_basis_functions,num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(RHS(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(q_tilde(num_basis_functions,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dlambdadomega(nomega_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_x(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_y(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_z(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dKDifferencedomega(3,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dLSE_current_density_with_areadOmega(nlambda,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  q_tilde = 0

  print *,"Beginning lse_sensitivity."
  call system_clock(tic_begin, countrate)

  ! Call LAPACK's DSYSV in query mode to determine the optimal size of the work array
  call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, RHS, num_basis_functions, WORK, -1, INFO)
  LWORK = WORK(1)
  print *,"Optimal LWORK:",LWORK
  deallocate(WORK)
  allocate(WORK(LWORK), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  ilambda = Nlambda
  solution = single_valued_current_potential_mn(:,ilambda)

  KDifference_x = d_x - matmul(f_x, solution)
  KDifference_y = d_y - matmul(f_y, solution)
  KDifference_z = d_z - matmul(f_z, solution)

  norm_K = K2(:,:,ilambda)**(0.5)

  sum_exp = sum(norm_normal_coil*dtheta_coil*dzeta_coil*nfp/area_coil*exp(target_option_p*(norm_K-max_K(ilambda))))

  call system_clock(tic, countrate)
  !$OMP PARALLEL
  !$OMP MASTER
  print *,"  Number of OpenMP threads:",omp_get_num_threads()
  !$OMP END MASTER
  !$OMP DO PRIVATE(dKDifferencedomega,dnorm_Kdomega,dsum_expdomega,dintKpdomega)
  do iomega = 1, nomega_coil
    dKDifferencedomega(1,:) = dddomega(1,iomega,1:ntheta_coil*nzeta_coil)-matmul(dfxdomega(iomega,:,:), solution)
    dKDifferencedomega(2,:) = dddomega(2,iomega,1:ntheta_coil*nzeta_coil)-matmul(dfydomega(iomega,:,:), solution)
    dKDifferencedomega(3,:) = dddomega(3,iomega,1:ntheta_coil*nzeta_coil)-matmul(dfzdomega(iomega,:,:), solution)
    dnorm_Kdomega = (1/(norm_K*(norm_normal_coil**2))) &
      * reshape(KDifference_x*dKDifferencedomega(1,:) + KDifference_y*dKDifferencedomega(2,:) &
      + KDifference_z*dKDifferencedomega(3,:),(/ ntheta_coil, nzeta_coil/)) &
      - norm_K*dnorm_normaldomega(iomega,:,:)/norm_normal_coil
    darea_coildomega(iomega) = dtheta_coil*dzeta_coil*nfp*sum(dnorm_normaldomega(iomega,:,:))
    if (target_option == 8) then
      dsum_expdomega = sum(dtheta_coil*dzeta_coil*nfp*exp(target_option_p*(norm_K-max_K(ilambda))) &
        *(dnorm_normaldomega(iomega,:,:)/area_coil &
        - norm_normal_coil*darea_coildomega(iomega)/(area_coil**2) &
        + norm_normal_coil*dnorm_Kdomega*target_option_p/(area_coil)))
      dtarget_optiondOmega(iomega,ilambda) = (1/target_option_p) * &
        dsum_expdomega/sum_exp
    end if
    if (target_option == 4) then
      dintKpdomega = sum(dtheta_coil*dzeta_coil*nfp*(target_option_p*norm_normal_coil*norm_K**(target_option_p-1)*dnorm_kdomega + dnorm_normaldomega(iomega,:,:)*norm_K**(target_option_p)))
      dtarget_optiondOmega(iomega,ilambda) = -(1/target_option_p)*L_p_norm_with_area(ilambda)*darea_coildomega(iomega)/area_coil + L_p_norm_with_area(ilambda)**(1-target_option_p)*dintKpdomega/(area_coil*target_option_p)
    end if
    if (target_option == 9) then
      dtarget_optiondOmega(iomega,ilambda) = dchi2Bdomega_withoutadjoint(iomega,ilambda)
    end if
    dLSE_current_density_with_areadOmega(iomega,ilambda) = -(1/target_option_p)*L_p_norm_with_area(ilambda)*darea_coildomega(iomega)/area_coil + L_p_norm_with_area(ilambda)**(1-target_option_p)*dintKpdomega/(area_coil*target_option_p)
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  call system_clock(toc)
  print *,"first loop over omega: ",real(toc-tic)/countrate," sec."

  call system_clock(tic,countrate)
  !$OMP PARALLEL
  !$OMP MASTER
  print *,"  Number of OpenMP threads:",omp_get_num_threads()
  !$OMP END MASTER
  !$OMP DO PRIVATE(dnorm_Kdphi,dsum_expdphi,dintKpdPhi)
  do i_basis_function=1,num_basis_functions
    dnorm_Kdphi = -1/((norm_normal_coil**2)*norm_k) &
      * reshape(KDifference_x*f_x(:,i_basis_function) + KDifference_y*f_y(:,i_basis_function) &
      + KDifference_z*f_z(:,i_basis_function), (/ ntheta_coil, nzeta_coil /))
    if (target_option == 8) then
      dsum_expdphi = sum(dtheta_coil*dzeta_coil*norm_normal_coil*nfp/area_coil* &
        exp(target_option_p*(norm_K - max_K(ilambda))) &
        *target_option_p*dnorm_Kdphi)
      dtarget_optiondPhi(i_basis_function,ilambda) = (1/target_option_p) * &
        dsum_expdphi/sum_exp
    end if
    if (target_option == 4) then
      dintKpdPhi = sum(dtheta_coil*dzeta_coil*nfp*norm_normal_coil*target_option_p*norm_K**(target_option_p-1)*dnorm_Kdphi)
      dtarget_optiondPhi(i_basis_function,ilambda) = (1/(target_option_p*area_coil)) * L_p_norm_with_area(ilambda)**(1-target_option_p)*dintKpdPhi
    end if
    if (target_option == 9) then
      dtarget_optiondPhi(i_basis_function,ilambda) = dchi2Bdphi(ilambda,i_basis_function)
    end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  call system_clock(toc)
  print *,"loop over basis_functions: ",real(toc-tic)/countrate," sec."

  call system_clock(tic)
  matrix = matrix_B + lambda(ilambda) * matrix_K
  RHS = dtarget_optiondPhi(:,ilambda)

  call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, &
    RHS, num_basis_functions, WORK, LWORK, INFO)
  if (INFO /= 0) then
    print *, "!!!!!! Error in LAPACK DSYSV: INFO = ", INFO
  end if
  q_tilde(:,ilambda) = RHS

  call system_clock(toc)
  print *,"adjoint solve: ",real(toc-tic)/countrate," sec."

  q_tilde_denom = dot_product((matmul(matrix_K,solution) - RHS_K),q_tilde(:,ilambda))

  call system_clock(tic)
  !$OMP PARALLEL
  !$OMP MASTER
  print *,"  Number of OpenMP threads:",omp_get_num_threads()
  !$OMP END MASTER
  !$OMP DO PRIVATE(q_tilde_num)
  do iomega=1,nomega_coil
    q_tilde_num = dot_product(dFdomega(iomega,:), q_tilde(:,ilambda))
    dlambdadomega(iomega,ilambda) = (1/q_tilde_denom)*(dtarget_optiondOmega(iomega,ilambda) - q_tilde_num)
    dchi2domega(iomega,ilambda) = dchi2domega(iomega,ilambda) + dlambdadomega(iomega,ilambda)*chi2_K(ilambda)
    if (sensitivity_option == 3) then
      dchi2Kdomega(iomega,ilambda) = dchi2Kdomega(iomega,ilambda) &
        - dot_product((matmul(matrix_K,solution) - RHS_K),q_K(:,ilambda))/q_tilde_denom &
        * (dtarget_optiondOmega(iomega,ilambda) - q_tilde_num)
      dchi2Bdomega(iomega,ilambda) = dchi2Bdomega(iomega,ilambda) &
        - dot_product((matmul(matrix_K,solution) - RHS_K),q_B(:,ilambda))/q_tilde_denom &
        * (dtarget_optiondOmega(iomega,ilambda) - q_tilde_num)
    end if
    if (sensitivity_option == 4) then
      dchi2Kdomega(iomega,ilambda) = dchi2Kdomega(iomega,ilambda) &
        - dot_product((matmul(matrix_K,solution) - RHS_K),q_K(:,ilambda))/q_tilde_denom &
        * (dtarget_optiondOmega(iomega,ilambda) - q_tilde_num)
      dchi2Bdomega(iomega,ilambda) = dchi2Bdomega(iomega,ilambda) &
        - dot_product((matmul(matrix_K,solution) - RHS_K),-lambda(ilambda)*q_K(:,ilambda))/q_tilde_denom &
        * (dtarget_optiondOmega(iomega,ilambda) - q_tilde_num)
    end if
    if (sensitivity_option == 5) then
      dchi2Bdomega(iomega,ilambda) = dchi2Bdomega(iomega,ilambda) &
        - dot_product((matmul(matrix_K,solution) - RHS_K),q_B(:,ilambda))/q_tilde_denom &
        * (dtarget_optiondOmega(iomega,ilambda) - q_tilde_num)
      dchi2Kdomega(iomega,ilambda) = dchi2Kdomega(iomega,ilambda) &
        - dot_product((matmul(matrix_K,solution) - RHS_K),-q_B(:,ilambda)/lambda(ilambda))/q_tilde_denom &
        * (dtarget_optiondOmega(iomega,ilambda) - q_tilde_num)
    end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  call system_clock(toc)
  print *,"second loop over omega: ",real(toc-tic)/countrate," sec."

  call system_clock(toc)
  print *,"lse_sensitivity: ",real(toc-tic_begin)/countrate," sec."

contains

end subroutine lse_sensitivity




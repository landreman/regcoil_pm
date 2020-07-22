subroutine regcoil_build_matrices()

  use regcoil_variables
  use stel_constants
  use stel_kinds
  use omp_lib
  use regcoil_init_Fourier_modes_mod
  use magpie_globals, only: stell_symm, tor_symm
  
  implicit none

  integer :: l_coil, l_mag, itheta_plasma, izeta_plasma, itheta_coil, izeta_coil, izetal_coil, indexl_mag, izeta_max
  real(dp) :: x, y, z, dx, dy, dz, dr2inv, dr32inv, reg_matrix_factor
  integer :: index_plasma, index_coil, j, imn
  integer :: tic, toc, tic1, toc1, toc2, countrate, iflag
  integer :: offset
  real(dp) :: angle, sinangle, cosangle, factor, constants, normal_plasma_dot_dr
  real(dp), dimension(:,:,:,:), allocatable :: g_over_N_plasma
  real(dp), dimension(:), allocatable :: norm_normal_plasma_inv1D
  integer :: js, ks, ls, row_offset, col_offset, j_RZetaZ, k_RZetaZ, k_RZetaZ_max, block_size
  real(dp), dimension(:,:), allocatable :: regularization_block, temp_matrix, Jacobian_coil_2D, regularization_without_RZetaZ
  real(dp), dimension(:,:), allocatable :: basis_functions_times_d
  real(dp), dimension(:), allocatable :: Bnormal_to_cancel_1D
  real(dp) :: thresh = 1.0d-17 ! Threshold for symmetry test when symmetry_option=2
  integer :: izeta_plasma_reflect, itheta_plasma_reflect, index_plasma_reflect, izeta_coil_reflect, itheta_coil_reflect, index_coil_reflect

  ! Variables needed by BLAS DGEMM:
  character :: TRANSA, TRANSB
  integer :: M, N, K, LDA, LDB, LDC
  real(dp) :: BLAS_ALPHA=1, BLAS_BETA=0
  real(dp), dimension(:,:), allocatable :: tempMatrix


  !--------------------------------------------------------------
  ! Initialize Jacobian of the magnetization region (if continuous)
  !--------------------------------------------------------------

  select case (trim(magnet_type))
  case ('continuous')
     if (allocated(Jacobian_coil)) deallocate(Jacobian_coil)
     allocate(Jacobian_coil(ntheta_coil, nzeta_coil, ns_integration))
     do js = 1, ns_integration
        Jacobian_coil(:,:,js) = d * (-norm_normal_coil + sign_normal * s_integration(js) * d * 2 * norm_normal_coil * mean_curvature_coil + s_integration(js) * s_integration(js) * (d * d * Jacobian_ssquared_term))
     end do
     if (any(Jacobian_coil >= 0)) then
        print *,"Error! Jacobian for the magnetization region is not negative-definite. max=",maxval(Jacobian_coil)
        print *,"maxval(d - max_d_before_singularity): ",maxval(d - max_d_before_singularity)
        do js = 1, ns_integration
           do izeta_coil = 1, nzeta_coil
              do itheta_coil = 1, ntheta_coil
                 if (Jacobian_coil(itheta_coil,izeta_coil,js) >= 0) then
                    print "(3(a,i3),5(a,es12.4))","itheta=",itheta_coil,"  izeta=",izeta_coil,"  js=",js,"  Jacobian_coil=",Jacobian_coil(itheta_coil,izeta_coil,js), &
                         "  d=",d(itheta_coil,izeta_coil),"  max_d=",max_d_before_singularity(itheta_coil,izeta_coil),"  mean_curvature_coil=",mean_curvature_coil(itheta_coil,izeta_coil),"  ssquared_term=",Jacobian_ssquared_term(itheta_coil,izeta_coil)
                 end if
              end do
           end do
        end do
        print *,Jacobian_coil
        stop
     end if
     Jacobian_coil = abs(Jacobian_coil) ! When we do volume integrals later, we want the absolute value of the Jacobian.
  case ('qhex')
  end select

  !--------------------------------------------------------------

  if (allocated(g)) deallocate(g)
  allocate(g(ntheta_plasma*nzeta_plasma, num_basis_functions, ns_magnetization, 3),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 8!'

  allocate(g_over_N_plasma(ntheta_plasma*nzeta_plasma, num_basis_functions, ns_magnetization, 3),stat=iflag)

  if (allocated(inductance)) deallocate(inductance)
  allocate(inductance(ntheta_plasma*nzeta_plasma, ntheta_coil*nzeta_coil, ns_magnetization, 3),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 9!'

  if (allocated(matrix_B)) deallocate(matrix_B)
  allocate(matrix_B(system_size, system_size),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 13!'

  if (allocated(matrix_regularization)) deallocate(matrix_regularization)
  allocate(matrix_regularization(system_size, system_size),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 14!'

  if (allocated(RHS_B)) deallocate(RHS_B)
  allocate(RHS_B(system_size),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 15!'

  if (allocated(RHS_regularization)) deallocate(RHS_regularization)
  allocate(RHS_regularization(system_size),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 16!'

  ! if (allocated(norm_normal_plasma_inv1D)) deallocate(norm_normal_plasma_inv1D)
  allocate(norm_normal_plasma_inv1D(ntheta_plasma*nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 17!'

  if (trim(magnet_type) == 'continuous') &
      allocate(d_times_unit_normal_coil(3,ntheta_coil,nzetal_coil))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now compute g
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call system_clock(tic,countrate)
  if (verbose) print *,"Building inductance matrix."

  inductance = 0

  allocate(cos_zetal(nzetal_coil))
  allocate(sin_zetal(nzetal_coil))
  do izetal_coil = 1, nzetal_coil
     cos_zetal(izetal_coil) = cos(zetal_coil(izetal_coil))
     sin_zetal(izetal_coil) = sin(zetal_coil(izetal_coil))
  end do
  constants = mu0 / (4 * pi)
      
  select case (trim(magnet_type))
  case('continuous')
     do l_coil = 0, (nfp-1)
        do izeta_coil = 1, nzeta_coil
           izetal_coil = izeta_coil + l_coil*nzeta_coil
           do itheta_coil = 1, ntheta_coil
              d_times_unit_normal_coil(:,itheta_coil,izetal_coil) = sign_normal * d(itheta_coil,izeta_coil) * normal_coil(:,itheta_coil,izetal_coil) / norm_normal_coil(itheta_coil,izeta_coil)
           end do
        end do
     end do
   case('qhex')
   end select

  if (symmetry_option==3 .or. symmetry_option==2) then
     izeta_max = nzeta_plasma
  else
     izeta_max = nzeta_plasma / 2 + 1 ! Note integer division
  end if

  !$OMP PARALLEL

  !$OMP MASTER
  if (verbose) print *,"  Number of OpenMP threads:",omp_get_num_threads()
  !$OMP END MASTER

  ! Note: the outermost loop below must be over the plasma variables rather than over the coil variables.
  ! This ensures the multiple threads write to different indices in h() rather than to the same indices in h(),
  ! in which case the h(index+plasma)=h(index_plasma)+... update does not work properly.
  !$OMP DO PRIVATE(index_plasma,index_coil,x,y,z,izetal_coil,dx,dy,dz,dr2inv,dr32inv,factor,normal_plasma_dot_dr)
  do izeta_plasma = 1, izeta_max
     do itheta_plasma = 1, ntheta_plasma
        index_plasma = (izeta_plasma-1)*ntheta_plasma + itheta_plasma
        x = r_plasma(1,itheta_plasma,izeta_plasma)
        y = r_plasma(2,itheta_plasma,izeta_plasma)
        z = r_plasma(3,itheta_plasma,izeta_plasma)

        do ks = 1, ns_magnetization
           do izeta_coil = 1, nzeta_coil
              do itheta_coil = 1, ntheta_coil
                 index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
                 do l_coil = 0, (nzetal_coil/nzeta_coil-1)
                    izetal_coil = izeta_coil + l_coil*nzeta_coil
                    do js = 1, ns_integration

                       select case (trim(magnet_type))
                       case ('continuous')
                          dx = x - (r_coil(1,itheta_coil,izetal_coil) + s_integration(js) * d_times_unit_normal_coil(1,itheta_coil,izetal_coil))
                          dy = y - (r_coil(2,itheta_coil,izetal_coil) + s_integration(js) * d_times_unit_normal_coil(2,itheta_coil,izetal_coil))
                          dz = z - (r_coil(3,itheta_coil,izetal_coil) + s_integration(js) * d_times_unit_normal_coil(3,itheta_coil,izetal_coil))
                          
                          dr2inv = 1/(dx*dx + dy*dy + dz*dz)
                          dr32inv = dr2inv*sqrt(dr2inv)
                          factor = dr32inv * Jacobian_coil(itheta_coil, izeta_coil, js) * s_weights(js) &
                               * interpolate_magnetization_to_integration(js, ks) * constants
                       case ('qhex')
                          dx = x - qhex_arr(izetal_coil)%ox
                          dy = y - qhex_arr(izetal_coil)%oy
                          dz = z - qhex_arr(izetal_coil)%oz
                          
                          dr2inv = 1/(dx*dx + dy*dy + dz*dz)
                          dr32inv = dr2inv*sqrt(dr2inv)
                          factor = dr32inv * qhex_arr(izetal_coil)%vol * constants
                       end select

                       normal_plasma_dot_dr = normal_plasma(1,itheta_plasma,izeta_plasma)*dx &
                            + normal_plasma(2,itheta_plasma,izeta_plasma)*dy &
                            + normal_plasma(3,itheta_plasma,izeta_plasma)*dz

                       ! R component of magnetization
                       ! e_R = e_X * cos(zeta) + e_Y * sin(zeta)
                       inductance(index_plasma,index_coil,ks,1) = inductance(index_plasma,index_coil,ks,1) - &
                            (normal_plasma(1,itheta_plasma,izeta_plasma)*cos_zetal(izetal_coil) &
                            +normal_plasma(2,itheta_plasma,izeta_plasma)*sin_zetal(izetal_coil) &
                            - (3*dr2inv) * normal_plasma_dot_dr * &
                            (cos_zetal(izetal_coil)*dx &
                            +sin_zetal(izetal_coil)*dy )) * factor

                       ! zeta component of magnetization
                       ! e_zeta = e_X * (-sin(zeta)) + e_Y * cos(zeta)
                       inductance(index_plasma,index_coil,ks,2) = inductance(index_plasma,index_coil,ks,2) - &
                            (normal_plasma(1,itheta_plasma,izeta_plasma)*(-sin_zetal(izetal_coil)) &
                            +normal_plasma(2,itheta_plasma,izeta_plasma)*  cos_zetal(izetal_coil) &
                            - (3*dr2inv) * normal_plasma_dot_dr * &
                            (-sin_zetal(izetal_coil)*dx &
                            + cos_zetal(izetal_coil)*dy )) * factor

                       ! Z component of magnetization
                       inductance(index_plasma,index_coil,ks,3) = inductance(index_plasma,index_coil,ks,3) - &
                            (normal_plasma(3,itheta_plasma,izeta_plasma) &
                            - (3*dr2inv) * normal_plasma_dot_dr * &
                            dz) * factor
                 
                    end do
                 end do
              end do
           end do
        end do

     end do
  end do
  !$OMP END DO


  if (symmetry_option==1 .or. symmetry_option==2) then
     izeta_max = nzeta_plasma / 2 + 1 ! Note integer division. We repeat this line here in case symmetry_option=2.
     !$OMP DO PRIVATE(izeta_plasma_reflect,itheta_plasma_reflect,index_plasma,index_plasma_reflect,izeta_coil_reflect,itheta_coil_reflect,index_coil,index_coil_reflect)
     do izeta_plasma = izeta_max+1, nzeta_plasma
        izeta_plasma_reflect = nzeta_plasma + 2 - izeta_plasma
        do itheta_plasma = 1, ntheta_plasma
           if (itheta_plasma == 1) then
              itheta_plasma_reflect = 1
           else
              itheta_plasma_reflect = ntheta_plasma + 2 - itheta_plasma
           end if
           index_plasma = (izeta_plasma-1)*ntheta_plasma + itheta_plasma
           index_plasma_reflect = (izeta_plasma_reflect-1)*ntheta_plasma + itheta_plasma_reflect

           select case (trim(magnet_type))

              case ('continuous')
                 do izeta_coil = 1, nzeta_coil
                    if (izeta_coil ==1) then
                       izeta_coil_reflect = izeta_coil
                    else
                       izeta_coil_reflect = nzeta_coil + 2 - izeta_coil
                    end if
                    do itheta_coil = 1, ntheta_coil
                       if (itheta_coil == 1) then
                          itheta_coil_reflect = 1
                       else
                          itheta_coil_reflect = ntheta_coil + 2 - itheta_coil
                       end if
                       index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
                       index_coil_reflect = (izeta_coil_reflect-1)*ntheta_coil + itheta_coil_reflect
      
                       if (symmetry_option==1) then
                          inductance(index_plasma,index_coil,:,1) =  inductance(index_plasma_reflect,index_coil_reflect,:,1)
                          inductance(index_plasma,index_coil,:,2) = -inductance(index_plasma_reflect,index_coil_reflect,:,2)
                          inductance(index_plasma,index_coil,:,3) = -inductance(index_plasma_reflect,index_coil_reflect,:,3)
                       else
                          ! symmetry_option must be 2. Verify symmetry:
                          do js = 1, ns_magnetization
                             if (abs(inductance(index_plasma,index_coil,js,1) - inductance(index_plasma_reflect,index_coil_reflect,js,1)) > thresh) then
                                print "(5(a,i3),3(a,es10.2))","discrepancy1 js:",js," izp:",izeta_plasma,' itp:',itheta_plasma,' izc:',izeta_coil,' itc:',itheta_coil, &
                                     ' i1:',inductance(index_plasma,index_coil,js,1),' i2:',inductance(index_plasma_reflect,index_coil_reflect,js,1), &
                                     ' diff:',inductance(index_plasma,index_coil,js,1) - inductance(index_plasma_reflect,index_coil_reflect,js,1)
                             end if
                             if (abs(inductance(index_plasma,index_coil,js,2) + inductance(index_plasma_reflect,index_coil_reflect,js,2)) > thresh) then
                                print "(5(a,i3),3(a,es10.2))","discrepancy2 js:",js," izp:",izeta_plasma,' itp:',itheta_plasma,' izc:',izeta_coil,' itc:',itheta_coil, &
                                     ' i1:',inductance(index_plasma,index_coil,js,2),' i2:',inductance(index_plasma_reflect,index_coil_reflect,js,2), &
                                     ' diff:',inductance(index_plasma,index_coil,js,2) + inductance(index_plasma_reflect,index_coil_reflect,js,2)
                             end if
                             if (abs(inductance(index_plasma,index_coil,js,3) + inductance(index_plasma_reflect,index_coil_reflect,js,3)) > thresh) then
                                print "(5(a,i3),3(a,es10.2))","discrepancy3 js:",js," izp:",izeta_plasma,' itp:',itheta_plasma,' izc:',izeta_coil,' itc:',itheta_coil, &
                                     ' i1:',inductance(index_plasma,index_coil,js,3),' i2:',inductance(index_plasma_reflect,index_coil_reflect,js,3), &
                                     ' diff:',inductance(index_plasma,index_coil,js,3) + inductance(index_plasma_reflect,index_coil_reflect,js,3)
                             end if
                          end do
                       end if
                    end do
                 end do

              case ('qhex')
                 do izeta_coil = 1, nzeta_coil
   
                    if (symmetry_option==1) then
                       inductance(index_plasma,izeta_coil,:,1) =  inductance(index_plasma_reflect,izeta_coil,:,1)
                       inductance(index_plasma,izeta_coil,:,2) = -inductance(index_plasma_reflect,izeta_coil,:,2)
                       inductance(index_plasma,izeta_coil,:,3) = -inductance(index_plasma_reflect,izeta_coil,:,3)
                    else
                    ! symmetry_option must be 2. Verify symmetry:
                       if (abs(inductance(index_plasma,index_coil,1,1) - inductance(index_plasma_reflect,index_coil_reflect,1,1)) > thresh) then
                          print "(3(a,i3),3(a,es10.2))","discrepancy1 izp:",izeta_plasma,' itp:',itheta_plasma,' im:',izeta_coil, &
                               ' i1:',inductance(index_plasma,index_coil,1,1),' i2:',inductance(index_plasma_reflect,index_coil_reflect,1,1), &
                               ' diff:',inductance(index_plasma,index_coil,1,1) - inductance(index_plasma_reflect,index_coil_reflect,1,1)
                       end if
                       if (abs(inductance(index_plasma,index_coil,1,2) + inductance(index_plasma_reflect,index_coil_reflect,1,2)) > thresh) then
                          print "(3(a,i3),3(a,es10.2))","discrepancy2 izp:",izeta_plasma,' itp:',itheta_plasma,' im:',izeta_coil, &
                               ' i1:',inductance(index_plasma,index_coil,1,2),' i2:',inductance(index_plasma_reflect,index_coil_reflect,1,2), &
                               ' diff:',inductance(index_plasma,index_coil,1,2) + inductance(index_plasma_reflect,index_coil_reflect,1,2)
                       end if
                       if (abs(inductance(index_plasma,index_coil,1,3) + inductance(index_plasma_reflect,index_coil_reflect,1,3)) > thresh) then
                          print "(3(a,i3),3(a,es10.2))","discrepancy3 izp:",izeta_plasma,' itp:',itheta_plasma,' im:',izeta_coil, &
                               ' i1:',inductance(index_plasma,index_coil,1,3),' i2:',inductance(index_plasma_reflect,index_coil_reflect,1,3), &
                               ' diff:',inductance(index_plasma,index_coil,1,3) + inductance(index_plasma_reflect,index_coil_reflect,1,3)
                       end if
                    end if

              end do
           
           end select

        end do
     end do
     !$OMP END DO
  end if

  !$OMP END PARALLEL

  deallocate(cos_zetal, sin_zetal)
  if (trim(magnet_type) == 'continuous') deallocate(d_times_unit_normal_coil)

  call system_clock(toc)
  if (verbose) print *,"Done. Took",real(toc-tic)/countrate,"sec."
    
  call system_clock(tic)

  ! For some reason, the BLAS matrix-matrix multiplication function DGEMM sometimes causes the
  ! program to crash on Edison unless you are careful to use the Intel MKL instead of Cray LibSci.
  ! If you like, you can use "matmul" instead which is slower but more reliable.

  !*******************************************************
  ! Call BLAS3 subroutine DGEMM for matrix multiplications:
  !*******************************************************

  allocate(Bnormal_to_cancel_1D(ntheta_plasma*nzeta_plasma))
  Bnormal_to_cancel_1D = reshape(Bnormal_from_TF_and_plasma_current, (/ ntheta_plasma*nzeta_plasma /))

  norm_normal_plasma_inv1D = reshape(1/norm_normal_plasma, (/ ntheta_plasma*nzeta_plasma /))
  if (verbose) print *,"Computing g and RHS_B:"

  g = 0
  do js = 1, ns_magnetization
     do j_RZetaZ = 1, 3

        call system_clock(tic1)

        select case (trim(magnet_type))
   
           case ('continuous')

      
              ! Here we carry out g = inductance * basis_functions
              ! A = inductance
              ! B = basis_functions
              ! C = g
              M = ntheta_plasma*nzeta_plasma ! # rows of A
              N = num_basis_functions ! # cols of B
              K = ntheta_coil*nzeta_coil ! Common dimension of A and B
              LDA = M
              LDB = K
              LDC = M
              TRANSA = 'N' ! No transposes
              TRANSB = 'N'
              BLAS_ALPHA=dtheta_coil*dzeta_coil
              BLAS_BETA=0 ! If other elements of g get zeroed out, try 1.
              if (j_RZetaZ==1) then
                 call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,inductance(:,:,js,j_RZetaZ),LDA,basis_functions_R,     LDB,BLAS_BETA,g(:,:,js,j_RZetaZ),LDC)
              else
                 call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,inductance(:,:,js,j_RZetaZ),LDA,basis_functions_zeta_Z,LDB,BLAS_BETA,g(:,:,js,j_RZetaZ),LDC)
              end if
      
           case ('qhex')

              ! No need to multiply inductance by basis functions (or by dtheta_coil*dzeta_coil)
              ! Also note ns_magnetization=1 for qhex configurations
              g(:,:,js,j_RZetaZ) = inductance(:,:,js,j_RZetaZ)

        end select

        call system_clock(toc1)
   
        do j = 1, num_basis_functions
           g_over_N_plasma(:,j,js,j_RZetaZ) = g(:,j,js,j_RZetaZ) * norm_normal_plasma_inv1D
        end do
      
        offset = (j_RZetaZ - 1)*ns_magnetization*num_basis_functions+ (js-1)*num_basis_functions
        RHS_B(offset+1:offset+num_basis_functions) = -dtheta_plasma*dzeta_plasma*matmul(Bnormal_to_cancel_1D, g(:,:,js,j_RZetaZ))
      
        call system_clock(toc2)
        if (verbose) print "(a,i4,a,i2,a,es11.3,a,es11.3,a)","  js=",js," j_RZetaZ=",j_RZetaZ," dgemm took",real(toc1-tic1)/countrate, &
            " sec, rest took",real(toc2-toc1)/countrate," sec."

     end do
  end do

  deallocate(Bnormal_to_cancel_1D)
  deallocate(norm_normal_plasma_inv1D)

  call system_clock(toc)
  if (verbose) print *,"Total time to compute g and RHS_B:",real(toc-tic)/countrate,"sec."


  matrix_B = 0

  if (verbose) print *,"About to form matrix_B."
  call system_clock(tic)

  allocate(temp_matrix(num_basis_functions,num_basis_functions))
  ! The logic of these loops is designed to reach (row=x,col=y) and (row=y,col=x) only once, to minimize matrix multiplications.
  ! Whichever one of the two cases is not reached directly is handled by adding the transpose of the block to the appropriate location.
  do js = 1, ns_magnetization
     do ks = 1, js
        do j_RZetaZ = 1, 3
           k_RZetaZ_max = 3
           if (js==ks) k_RZetaZ_max = j_RZetaZ
           do k_RZetaZ = 1, k_RZetaZ_max
              row_offset = (j_RZetaZ-1)*ns_magnetization*num_basis_functions + (js-1)*num_basis_functions
              col_offset = (k_RZetaZ-1)*ns_magnetization*num_basis_functions + (ks-1)*num_basis_functions
              call system_clock(tic1)

              ! Here we carry out matrix_B = (dtheta*dzeta)*(g ^ T) * g_over_N_plasma
              ! A = g
              ! B = g_over_N_plasma
              ! C = contribution to matrix_B
              M = num_basis_functions ! # rows of A^T
              N = num_basis_functions ! # cols of B
              K = ntheta_plasma*nzeta_plasma ! Common dimension of A^T and B
              LDA = K ! Would be M if not taking the transpose.
              LDB = K
              LDC = M
              TRANSA = 'T' ! DO take a transpose!
              TRANSB = 'N'
              BLAS_ALPHA = dtheta_plasma*dzeta_plasma
              BLAS_BETA=0
              call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,g(:,:,js,j_RZetaZ),LDA,g_over_N_plasma(:,:,ks,k_RZetaZ),LDB,BLAS_BETA,temp_matrix,LDC)

              matrix_B(row_offset+1:row_offset+num_basis_functions, col_offset+1:col_offset+num_basis_functions) = temp_matrix
              if ((js.ne.ks) .or. (j_RZetaZ.ne.k_RZetaZ)) then
                 !matrix_B(col_offset+1:col_offset+num_basis_functions, row_offset+1:row_offset+num_basis_functions) = temp_matrix
                 matrix_B(col_offset+1:col_offset+num_basis_functions, row_offset+1:row_offset+num_basis_functions) = transpose(temp_matrix)
              end if
              call system_clock(toc1)

              if (verbose) print "(a,i4,a,i4,a,i2,a,i2,a,es11.3,a)","  js=",js,", ks=",ks,", j_RZetaZ=",j_RZetaZ,", k_RZetaZ=",k_RZetaZ, &
                   " dgemm took",real(toc1-tic1)/countrate, " sec."
           end do
        end do
     end do
  end do
  deallocate(temp_matrix)

  call system_clock(toc)
  if (verbose) print *,"Total time for matrix_B:",real(toc-tic)/countrate,"sec."

  deallocate(g_over_N_plasma)
    

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! Assemble the matrix for the regularization term
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  call system_clock(tic)
  RHS_regularization = 0
  matrix_regularization = 0

  select case (trim(magnet_type))
  case ('continuous')
   
     allocate(regularization_block(num_basis_functions,num_basis_functions))
     allocate(temp_matrix(ntheta_coil*nzeta_coil,num_basis_functions))
     allocate(Jacobian_coil_2D(ntheta_coil*nzeta_coil,ns_integration))
     allocate(regularization_without_RZetaZ(num_basis_functions*ns_magnetization, num_basis_functions*ns_magnetization))
     do js = 1, ns_integration
        Jacobian_coil_2D(:,js) = reshape(Jacobian_coil(:,:,js), (/ ntheta_coil*nzeta_coil /))
     end do
   
     allocate(basis_functions_times_d(ntheta_coil*nzeta_coil,num_basis_functions))
     do izeta_coil = 1, nzeta_coil
        do itheta_coil = 1, ntheta_coil
           index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
           basis_functions_times_d(index_coil,:) = (d(itheta_coil,izeta_coil) ** regularization_d_exponent) * basis_functions_R(index_coil,:) * ports_weight(itheta_coil,izeta_coil)
        end do
     end do
   
     regularization_without_RZetaZ = 0
     ! Add contributions from each integration point in s:
     do ls = 1, ns_integration
        regularization_block = 0
        ! Form basis_functions^T * d * sqrt(g) * basis_functions, to do the (theta,zeta) integrations.
        do j = 1, num_basis_functions 
           temp_matrix(:,j) = Jacobian_coil_2D(:,ls) * basis_functions_R(:,j)
        end do
        ! A = basis_functions_times_d
        ! B = temp_matrix
        ! C = regularization_block
        M = num_basis_functions ! # rows of A^T
        N = num_basis_functions ! # cols of B
        K = ntheta_coil*nzeta_coil ! Common dimension of A^T and B
        LDA = K ! Would be M if not taking the transpose.
        LDB = K
        LDC = M
        TRANSA = 'T' ! DO take a transpose!
        TRANSB = 'N'
        BLAS_ALPHA = dtheta_coil * dzeta_coil * s_weights(ls)
        BLAS_BETA=0
        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,basis_functions_times_d,LDA,temp_matrix,LDB,BLAS_BETA,regularization_block,LDC)           
        do js = 1, ns_magnetization
           do ks = 1, js
              row_offset = (js-1)*num_basis_functions
              col_offset = (ks-1)*num_basis_functions
              regularization_without_RZetaZ(row_offset+1:row_offset+num_basis_functions, col_offset+1:col_offset+num_basis_functions) &
                   = regularization_without_RZetaZ(row_offset+1:row_offset+num_basis_functions, col_offset+1:col_offset+num_basis_functions) &
                   + regularization_block * interpolate_magnetization_to_integration(ls,js) * interpolate_magnetization_to_integration(ls,ks)
              if (ks .ne. js) then
                 ! For off-diagonal-in-s blocks, populate the symmetric block.
                 ! Should I take a transpose of regularization_block in the next line? I think it doesn't matter, and it's already symmetric.
                 regularization_without_RZetaZ(col_offset+1:col_offset+num_basis_functions, row_offset+1:row_offset+num_basis_functions) &
                      = regularization_without_RZetaZ(col_offset+1:col_offset+num_basis_functions, row_offset+1:row_offset+num_basis_functions) &
                      + regularization_block * interpolate_magnetization_to_integration(ls,js) * interpolate_magnetization_to_integration(ls,ks)
              end if
           end do
        end do
     end do
     
     block_size = ns_magnetization*num_basis_functions
     select case (symmetry_option)
     case (3)
   
        ! Make copies for (R,Phi,Z) components of the magnetization:
        do j_RZetaZ = 0,2
           row_offset = j_RZetaZ*block_size
           matrix_regularization(row_offset+1:row_offset+block_size, row_offset+1:row_offset+block_size) = regularization_without_RZetaZ
        end do
   
     case (1,2)
        row_offset = 0
        matrix_regularization(row_offset+1:row_offset+block_size, row_offset+1:row_offset+block_size) = regularization_without_RZetaZ
   
        ! Now repeat all that for the zeta_Z basis functions:
        do izeta_coil = 1, nzeta_coil
           do itheta_coil = 1, ntheta_coil
              index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
              basis_functions_times_d(index_coil,:) = (d(itheta_coil,izeta_coil) ** regularization_d_exponent) * basis_functions_zeta_Z(index_coil,:) * ports_weight(itheta_coil,izeta_coil)
           end do
        end do
   
        regularization_without_RZetaZ = 0
        ! Add contributions from each integration point in s:
        do ls = 1, ns_integration
           regularization_block = 0
           ! Form basis_functions^T * d * sqrt(g) * basis_functions, to do the (theta,zeta) integrations.
           do j = 1, num_basis_functions 
              temp_matrix(:,j) = Jacobian_coil_2D(:,ls) * basis_functions_zeta_Z(:,j)
           end do
           ! A = basis_functions_times_d
           ! B = temp_matrix
           ! C = regularization_block
           M = num_basis_functions ! # rows of A^T
           N = num_basis_functions ! # cols of B
           K = ntheta_coil*nzeta_coil ! Common dimension of A^T and B
           LDA = K ! Would be M if not taking the transpose.
           LDB = K
           LDC = M
           TRANSA = 'T' ! DO take a transpose!
           TRANSB = 'N'
           BLAS_ALPHA = dtheta_coil * dzeta_coil * s_weights(ls)
           BLAS_BETA=0
           call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,basis_functions_times_d,LDA,temp_matrix,LDB,BLAS_BETA,regularization_block,LDC)           
           do js = 1, ns_magnetization
              do ks = 1, js
                 row_offset = (js-1)*num_basis_functions
                 col_offset = (ks-1)*num_basis_functions
                 regularization_without_RZetaZ(row_offset+1:row_offset+num_basis_functions, col_offset+1:col_offset+num_basis_functions) &
                      = regularization_without_RZetaZ(row_offset+1:row_offset+num_basis_functions, col_offset+1:col_offset+num_basis_functions) &
                      + regularization_block * interpolate_magnetization_to_integration(ls,js) * interpolate_magnetization_to_integration(ls,ks)
                 if (ks .ne. js) then
                    ! For off-diagonal-in-s blocks, populate the symmetric block.
                    ! Should I take a transpose of regularization_block in the next line? I think it doesn't matter, and it's already symmetric.
                    regularization_without_RZetaZ(col_offset+1:col_offset+num_basis_functions, row_offset+1:row_offset+num_basis_functions) &
                         = regularization_without_RZetaZ(col_offset+1:col_offset+num_basis_functions, row_offset+1:row_offset+num_basis_functions) &
                         + regularization_block * interpolate_magnetization_to_integration(ls,js) * interpolate_magnetization_to_integration(ls,ks)
                 end if
              end do
           end do
        end do
   
        ! Make copies for (Phi,Z) components of the magnetization:
        do j_RZetaZ = 1,2
           row_offset = j_RZetaZ*block_size
           matrix_regularization(row_offset+1:row_offset+block_size, row_offset+1:row_offset+block_size) = regularization_without_RZetaZ
        end do
   
     case default
        stop "Should not get here!"
     end select
   
     deallocate(regularization_block, temp_matrix, Jacobian_coil_2D)
     deallocate(regularization_without_RZetaZ, basis_functions_times_d)
   
     call system_clock(toc)
     if (verbose) print *,"Regularization term:",real(toc-tic)/countrate,"sec."

  case ('qhex')

    ! If the solution represents only one half-period's worth of magnets, multiply by 2 to
    ! weight correctly relative to chiSq_B
    if (stell_symm) then
        reg_matrix_factor = 2._dp
    else
        reg_matrix_factor = 1._dp
    end if

    do j_RZetaZ = 1, 3

       offset = (j_RZetaZ-1)*num_basis_functions

       do j = 1, num_basis_functions
          matrix_regularization(offset+j,offset+j) = reg_matrix_factor * qhex_arr(j)%vol * d(1,j) ** regularization_d_exponent
       end do

    end do

  end select



end subroutine regcoil_build_matrices




! Documentation of BLAS3 DGEMM subroutine for matrix-matrix multiplication:

!!$*       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!!$* 
!!$*       .. Scalar Arguments ..
!!$*       DOUBLE PRECISION ALPHA,BETA
!!$*       INTEGER K,LDA,LDB,LDC,M,N
!!$*       CHARACTER TRANSA,TRANSB
!!$*       ..
!!$*       .. Array Arguments ..
!!$*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!!$*       ..
!!$*  
!!$*
!!$*> \par Purpose:
!!$*  =============
!!$*>
!!$*> \verbatim
!!$*>
!!$*> DGEMM  performs one of the matrix-matrix operations
!!$*>
!!$*>    C := alpha*op( A )*op( B ) + beta*C,
!!$*>
!!$*> where  op( X ) is one of
!!$*>
!!$*>    op( X ) = X   or   op( X ) = X**T,
!!$*>
!!$*> alpha and beta are scalars, and A, B and C are matrices, with op( A )
!!$*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!!$*> \endverbatim
!!$*
!!$*  Arguments:
!!$*  ==========
!!$*
!!$*> \param[in] TRANSA
!!$*> \verbatim
!!$*>          TRANSA is CHARACTER*1
!!$*>           On entry, TRANSA specifies the form of op( A ) to be used in
!!$*>           the matrix multiplication as follows:
!!$*>
!!$*>              TRANSA = 'N' or 'n',  op( A ) = A.
!!$*>
!!$*>              TRANSA = 'T' or 't',  op( A ) = A**T.
!!$*>
!!$*>              TRANSA = 'C' or 'c',  op( A ) = A**T.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] TRANSB
!!$*> \verbatim
!!$*>          TRANSB is CHARACTER*1
!!$*>           On entry, TRANSB specifies the form of op( B ) to be used in
!!$*>           the matrix multiplication as follows:
!!$*>
!!$*>              TRANSB = 'N' or 'n',  op( B ) = B.
!!$*>
!!$*>              TRANSB = 'T' or 't',  op( B ) = B**T.
!!$*>
!!$*>              TRANSB = 'C' or 'c',  op( B ) = B**T.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] M
!!$*> \verbatim
!!$*>          M is INTEGER
!!$*>           On entry,  M  specifies  the number  of rows  of the  matrix
!!$*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] N
!!$*> \verbatim
!!$*>          N is INTEGER
!!$*>           On entry,  N  specifies the number  of columns of the matrix
!!$*>           op( B ) and the number of columns of the matrix C. N must be
!!$*>           at least zero.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] K
!!$*> \verbatim
!!$*>          K is INTEGER
!!$*>           On entry,  K  specifies  the number of columns of the matrix
!!$*>           op( A ) and the number of rows of the matrix op( B ). K must
!!$*>           be at least  zero.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] ALPHA
!!$*> \verbatim
!!$*>          ALPHA is DOUBLE PRECISION.
!!$*>           On entry, ALPHA specifies the scalar alpha.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] A
!!$*> \verbatim
!!$*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!!$*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!!$*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!!$*>           part of the array  A  must contain the matrix  A,  otherwise
!!$*>           the leading  k by m  part of the array  A  must contain  the
!!$*>           matrix A.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDA
!!$*> \verbatim
!!$*>          LDA is INTEGER
!!$*>           On entry, LDA specifies the first dimension of A as declared
!!$*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!!$*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!!$*>           least  max( 1, k ).
!!$*> \endverbatim
!!$*>
!!$*> \param[in] B
!!$*> \verbatim
!!$*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!!$*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!!$*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!!$*>           part of the array  B  must contain the matrix  B,  otherwise
!!$*>           the leading  n by k  part of the array  B  must contain  the
!!$*>           matrix B.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDB
!!$*> \verbatim
!!$*>          LDB is INTEGER
!!$*>           On entry, LDB specifies the first dimension of B as declared
!!$*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!!$*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!!$*>           least  max( 1, n ).
!!$*> \endverbatim
!!$*>
!!$*> \param[in] BETA
!!$*> \verbatim
!!$*>          BETA is DOUBLE PRECISION.
!!$*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!!$*>           supplied as zero then C need not be set on input.
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] C
!!$*> \verbatim
!!$*>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!!$*>           Before entry, the leading  m by n  part of the array  C must
!!$*>           contain the matrix  C,  except when  beta  is zero, in which
!!$*>           case C need not be set on entry.
!!$*>           On exit, the array  C  is overwritten by the  m by n  matrix
!!$*>           ( alpha*op( A )*op( B ) + beta*C ).
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDC
!!$*> \verbatim
!!$*>          LDC is INTEGER
!!$*>           On entry, LDC specifies the first dimension of C as declared
!!$*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
!!$*>           max( 1, m ).
!!$*> \endverbatim




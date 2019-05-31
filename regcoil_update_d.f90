subroutine regcoil_update_d(jd,isaved)

  use regcoil_variables
  use stel_constants

  implicit none

  integer, intent(in) :: jd, isaved
  integer :: Anderson_size, Anderson_k, i, j
  real(dp), dimension(:), allocatable :: Anderson_RHS, Anderson_solution
  real(dp), dimension(:,:), allocatable :: Anderson_matrix

  ! Update thickness:
  last_d = d
  d = last_d * s_averaged_abs_M(:,:,isaved) / (target_mu0_M / mu0)
  ! Take a mixture of the new and old depths:
  d = Picard_alpha * d + (1 - Picard_alpha) * last_d

  select case (trim(d_option))
  case (d_option_uniform)
     stop "Should not get here"
  case (d_option_Picard)
  case (d_option_Anderson)
     do j = 1, Anderson_depth
        Anderson_G(:,:,j) = Anderson_G(:,:,j+1)
        Anderson_u_tilde(:,:,j) = Anderson_u_tilde(:,:,j+1)
     end do
     Anderson_u_tilde(:,:,Anderson_depth+1) = d
     Anderson_G(:,:,Anderson_depth+1) = d - last_d

     if (jd > 1) then
        ! Need to write this part.
        ! Sanchez-Vizuet notation        Notation here
        ! --------------------------------------------
        !                n                       jd-1
        !                k                 Anderson_k
        !                m             Anderson_depth
        !                u                          d
        !         \tilde{u}          Anderson_u_tilde


        Anderson_k = min(Anderson_depth, jd-1)
        Anderson_size = Anderson_k + 2

        allocate(Anderson_matrix(Anderson_size,Anderson_size))
        allocate(Anderson_RHS(Anderson_size))
        allocate(Anderson_solution(Anderson_size))

        Anderson_RHS = 0
        Anderson_RHS(Anderson_size) = 1

        Anderson_matrix = 0
        Anderson_matrix(Anderson_size, 1:Anderson_size-1) = 1
        Anderson_matrix(1:Anderson_size-1, Anderson_size) = 1
        do j = 1, Anderson_size-1
           do i = 1, Anderson_size-1
              Anderson_matrix(j,i) = 2 * dtheta_coil * dzeta_coil * sum(Anderson_G(:,:,j-Anderson_size+Anderson_depth+2) * Anderson_G(:,:,i-Anderson_size+Anderson_depth+2))
           end do
        end do

        print *,"jd:",jd,"  Anderson_depth:",Anderson_depth,"  Anderson_size:",Anderson_size
        print *,"Anderson_matrix:"
        do j = 1,Anderson_size
           print *,Anderson_matrix(j,:)
        end do

        ! Compute solution = matrix \ RHS.
        ! Use LAPACK's DSYSV since matrix is symmetric.
        ! Note: RHS will be over-written with the solution.
        call DSYSV('U',Anderson_size, 1, Anderson_matrix, Anderson_size, LAPACK_IPIV, Anderson_RHS, Anderson_size, LAPACK_WORK, LAPACK_LWORK, LAPACK_INFO)
        if (LAPACK_INFO /= 0) then
           print *, "!!!!!! Error in LAPACK DSYSV: INFO = ", LAPACK_INFO
           !stop
        end if
        Anderson_solution = Anderson_RHS
        print *,"Anderson solution:",Anderson_solution

        d = 0
        do j = 1, Anderson_size-1
           d = d + Anderson_solution(j) * Anderson_u_tilde(:,:,j - Anderson_size + Anderson_depth + 2)
        end do

        ! Take a mixture of the Anderson update and the Picard update
        d = Anderson_alpha * d + (1 - Anderson_alpha) * Anderson_u_tilde(:,:,Anderson_depth+1)

        deallocate(Anderson_matrix,Anderson_RHS,Anderson_solution)

     end if
  end select

  if (verbose) print *,"Updated d. ||d_old - d_new|| = ",dtheta_coil*dzeta_coil*sum((d - last_d) **2)

  print *,"new d:"
  do j = 1,ntheta_coil
     print "(*(f5.2))",d(j,:)
  end do

  
end subroutine regcoil_update_d

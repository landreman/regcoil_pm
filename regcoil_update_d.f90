subroutine regcoil_update_d(jd,isaved)

  use regcoil_variables
  use stel_constants

  implicit none

  integer, intent(in) :: jd, isaved
  integer :: Anderson_size, Anderson_k, i, j, itheta, izeta
  real(dp), dimension(:), allocatable :: Anderson_RHS, Anderson_solution
  real(dp), dimension(:,:), allocatable :: Anderson_matrix
  integer :: M, N, NRHS, RANK, Anderson_LWORK, Anderson_LIWORK
  real(dp) :: RCOND, diff
  real(dp) :: d0, factor, factor2, angle, sinangle, cosangle
  real(dp), dimension(:), allocatable :: singular_values, Anderson_WORK
  integer, dimension(:), allocatable :: Anderson_IWORK
  real(dp), dimension(:), allocatable :: dmnc, dmns

  ! Update thickness:
  last_d = d
  d = last_d * abs_M(:,:,1,isaved) / (target_mu0_M / mu0)
!  d = last_d * s_averaged_abs_M(:,:,isaved) / (target_mu0_M / mu0)

  ! Truncate d wherever it becomes too large:
  diff = maxval(d - max_d_before_singularity)
  if (diff>0) print "(a,es12.4,a)"," Note: lowering d (at step 1) since maxval(d - max_d_before_singularity) =",diff," is positive."
  d = min(d,max_d_before_singularity)

  if (min_d > 0) then
     if (verbose) print *,"Applying min_d=",min_d
     do i = 1, ntheta_coil
        do j = 1, nzeta_coil
           d(i,j) = d(i,j) + min_d * exp(-d(i,j) / min_d)
        end do
     end do
  end if

  ! Fourier-filter d. First we transform from real to Fourier space:
!!$  print *,"d before filtering:"
!!$  do j=1,ntheta_coil
!!$     print "(*(f7.4))",d(j,:)
!!$  end do
  allocate(dmnc(mnmax_magnetization))
  allocate(dmns(mnmax_magnetization))
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
  ! Now inverse transform:
  d = d0
  do izeta = 1, nzeta_coil
     do itheta = 1, ntheta_coil
        do j = 1, mnmax_magnetization
           angle = xm_magnetization(j)*theta_coil(itheta)-xn_magnetization(j)*zeta_coil(izeta)
           sinangle = sin(angle)
           cosangle = cos(angle)
           d(itheta,izeta) = d(itheta,izeta) + dmnc(j)*cosangle + dmns(j)*sinangle
        end do
     end do
  end do
!!$  print *,"d after filtering:"
!!$  do j=1,ntheta_coil
!!$     print "(*(f7.4))",d(j,:)
!!$  end do

  ! Take a mixture of the new and old depths:
  d = Picard_alpha * d + (1 - Picard_alpha) * last_d

  select case (trim(d_option))
  case (d_option_uniform)
     stop "Should not get here"
  case (d_option_Picard)
  case (d_option_Anderson_old)
     do j = 1, Anderson_depth
        Anderson_G(:,:,j) = Anderson_G(:,:,j+1)
        Anderson_u_tilde(:,:,j) = Anderson_u_tilde(:,:,j+1)
     end do
     Anderson_u_tilde(:,:,Anderson_depth+1) = d
     Anderson_G(:,:,Anderson_depth+1) = d - last_d

     if (jd > 1) then
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
           print "(*(es10.3,2x))",Anderson_matrix(j,:)
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

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  case (d_option_Anderson)
     do j = 1, Anderson_depth
        Anderson_G(:,:,j) = Anderson_G(:,:,j+1)
        Anderson_u_tilde(:,:,j) = Anderson_u_tilde(:,:,j+1)
     end do
     Anderson_u_tilde(:,:,Anderson_depth+1) = d
     Anderson_G(:,:,Anderson_depth+1) = d - last_d

     if (jd > 1) then
        ! Sanchez-Vizuet notation        Notation here
        ! --------------------------------------------
        !                n                       jd-1
        !                k                 Anderson_k
        !                m             Anderson_depth
        !                u                          d
        !         \tilde{u}          Anderson_u_tilde

        Anderson_k = min(Anderson_depth, jd-1)

        allocate(Anderson_matrix(ntheta_coil * nzeta_coil, Anderson_k))
        allocate(Anderson_RHS(ntheta_coil * nzeta_coil))
        allocate(Anderson_solution(Anderson_k))
        allocate(singular_values(Anderson_k))

        Anderson_RHS = reshape(Anderson_G(:,:,Anderson_depth+1), (/ ntheta_coil * nzeta_coil /))
        do j = 1, Anderson_k
           Anderson_matrix(:,j) = Anderson_RHS - reshape(Anderson_G(:,:,Anderson_depth+1+j-Anderson_k-1), (/ ntheta_coil * nzeta_coil /))
        end do

        print *,"jd:",jd,"  Anderson_depth:",Anderson_depth,"  Anderson_k:",Anderson_k
!!$        print *,"Anderson_matrix:"
!!$        do j = 1,ntheta_coil*nz
!!$           print "(*(es10.3,2x))",Anderson_matrix(j,:)
!!$        end do

        M = ntheta_coil * nzeta_coil
        N = Anderson_k
        if (N > M) stop "N>M in Anderson acceleration. This case should not happen!"
        
        NRHS = 1
        RCOND = 1.0d-13

        ! Determine sizes of work arrays:
        Anderson_LWORK = -1
        allocate(Anderson_WORK(1))
        allocate(Anderson_IWORK(1))
        ! DGELSD is a LAPACK subroutine for solving potentially-rank-deficient linear-least-squares problems.
        call DGELSD( M, N, NRHS, Anderson_matrix, M, Anderson_RHS, M, singular_values, RCOND, RANK, Anderson_WORK, Anderson_LWORK, Anderson_IWORK, LAPACK_INFO)
        if (LAPACK_INFO /= 0) then
           print *, "!!!!!! Error in LAPACK DGELSD determining work array sizes: INFO = ", LAPACK_INFO
           stop
        end if
        Anderson_LWORK = Anderson_WORK(1) * 2  ! Documentation says "For good performance, LWORK should generally be larger" than the minimum.
        Anderson_LIWORK = max(1000, 3 * N * MAX( 0, INT( LOG( N/(100+1.0) ) / log(2.0) ) + 1 )  + 11 * N)
        print *,"Anderson_LWORK:",Anderson_LWORK," Anderson_LIWORK:",Anderson_LIWORK
        deallocate(Anderson_WORK, Anderson_IWORK)
        allocate(Anderson_WORK(Anderson_LWORK))
        allocate(Anderson_IWORK(Anderson_LIWORK))

        ! Note: Anderson_matrix and Anderson_RHS are over-written by LAPACK.
        call DGELSD( M, N, NRHS, Anderson_matrix, M, Anderson_RHS, M, singular_values, RCOND, RANK, Anderson_WORK, Anderson_LWORK, Anderson_IWORK, LAPACK_INFO)
        if (LAPACK_INFO /= 0) then
           print *, "!!!!!! Error in LAPACK DGELSD: INFO = ", LAPACK_INFO
           stop
        end if
        Anderson_solution = Anderson_RHS(1:Anderson_k)
        print *,"Anderson_LIWORK:",Anderson_IWORK(1)
        print *,"Anderson solution:",Anderson_solution,1-sum(Anderson_solution)
        print *,"Numerical rank:",RANK
        print *,"Singular values:",singular_values

        ! Last alpha = 1 - sum(previous alphas)
        d = (1 - sum(Anderson_solution)) * Anderson_u_tilde(:,:,Anderson_depth+1)
        do j = 1, Anderson_k
           d = d + Anderson_solution(j) * Anderson_u_tilde(:,:,j - Anderson_k + Anderson_depth)
        end do

        ! Take a mixture of the Anderson update and the Picard update
        d = Anderson_alpha * d + (1 - Anderson_alpha) * Anderson_u_tilde(:,:,Anderson_depth+1)

        deallocate(Anderson_matrix,Anderson_RHS,Anderson_solution,singular_values,Anderson_WORK,Anderson_IWORK)

     end if
  end select

  ! Once again, truncate d wherever it becomes too large:
  diff = maxval(d - max_d_before_singularity)
  if (diff>0) print "(a,es12.4,a)"," Note: lowering d (at step 2) since maxval(d - max_d_before_singularity) =",diff," is positive."
  d = min(d,max_d_before_singularity)

  ! Force d to be positive:
  diff = minval(d)
  if (diff <= 0) print "(a,es12.4,a)"," Note: increasing d since min(d) =",diff," is negative."
  do izeta = 1, nzeta_coil
     do itheta = 1, ntheta_coil
        d(itheta,izeta) = max(d(itheta,izeta), 1.0d-14)
     end do
  end do

  if (verbose) print *,"Updated d. ||d_old - d_new|| = ",dtheta_coil*dzeta_coil*sum((d - last_d) **2)

!!$  print *,"new d:"
!!$  do j = 1,ntheta_coil
!!$     print "(*(f5.2))",d(j,:)
!!$  end do

  
end subroutine regcoil_update_d


!!$*> \brief <b> DGELSD computes the minimum-norm solution to a linear least squares problem for GE matrices</b>
!!$*
!!$*  =========== DOCUMENTATION ===========
!!$*
!!$* Online html documentation available at 
!!$*            http://www.netlib.org/lapack/explore-html/ 
!!$*
!!$*> \htmlonly
!!$*> Download DGELSD + dependencies 
!!$*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelsd.f"> 
!!$*> [TGZ]</a> 
!!$*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelsd.f"> 
!!$*> [ZIP]</a> 
!!$*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelsd.f"> 
!!$*> [TXT]</a>
!!$*> \endhtmlonly 
!!$*
!!$*  Definition:
!!$*  ===========
!!$*
!!$*       SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
!!$*                          WORK, LWORK, IWORK, INFO )
!!$* 
!!$*       .. Scalar Arguments ..
!!$*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
!!$*       DOUBLE PRECISION   RCOND
!!$*       ..
!!$*       .. Array Arguments ..
!!$*       INTEGER            IWORK( * )
!!$*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
!!$*       ..
!!$*  
!!$*
!!$*> \par Purpose:
!!$*  =============
!!$*>
!!$*> \verbatim
!!$*>
!!$*> DGELSD computes the minimum-norm solution to a real linear least
!!$*> squares problem:
!!$*>     minimize 2-norm(| b - A*x |)
!!$*> using the singular value decomposition (SVD) of A. A is an M-by-N
!!$*> matrix which may be rank-deficient.
!!$*>
!!$*> Several right hand side vectors b and solution vectors x can be
!!$*> handled in a single call; they are stored as the columns of the
!!$*> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!!$*> matrix X.
!!$*>
!!$*> The problem is solved in three steps:
!!$*> (1) Reduce the coefficient matrix A to bidiagonal form with
!!$*>     Householder transformations, reducing the original problem
!!$*>     into a "bidiagonal least squares problem" (BLS)
!!$*> (2) Solve the BLS using a divide and conquer approach.
!!$*> (3) Apply back all the Householder tranformations to solve
!!$*>     the original least squares problem.
!!$*>
!!$*> The effective rank of A is determined by treating as zero those
!!$*> singular values which are less than RCOND times the largest singular
!!$*> value.
!!$*>
!!$*> The divide and conquer algorithm makes very mild assumptions about
!!$*> floating point arithmetic. It will work on machines with a guard
!!$*> digit in add/subtract, or on those binary machines without guard
!!$*> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!!$*> Cray-2. It could conceivably fail on hexadecimal or decimal machines
!!$*> without guard digits, but we know of none.
!!$*> \endverbatim
!!$*
!!$*  Arguments:
!!$*  ==========
!!$*
!!$*> \param[in] M
!!$*> \verbatim
!!$*>          M is INTEGER
!!$*>          The number of rows of A. M >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] N
!!$*> \verbatim
!!$*>          N is INTEGER
!!$*>          The number of columns of A. N >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] NRHS
!!$*> \verbatim
!!$*>          NRHS is INTEGER
!!$*>          The number of right hand sides, i.e., the number of columns
!!$*>          of the matrices B and X. NRHS >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] A
!!$*> \verbatim
!!$*>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!$*>          On entry, the M-by-N matrix A.
!!$*>          On exit, A has been destroyed.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDA
!!$*> \verbatim
!!$*>          LDA is INTEGER
!!$*>          The leading dimension of the array A.  LDA >= max(1,M).
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] B
!!$*> \verbatim
!!$*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!!$*>          On entry, the M-by-NRHS right hand side matrix B.
!!$*>          On exit, B is overwritten by the N-by-NRHS solution
!!$*>          matrix X.  If m >= n and RANK = n, the residual
!!$*>          sum-of-squares for the solution in the i-th column is given
!!$*>          by the sum of squares of elements n+1:m in that column.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDB
!!$*> \verbatim
!!$*>          LDB is INTEGER
!!$*>          The leading dimension of the array B. LDB >= max(1,max(M,N)).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] S
!!$*> \verbatim
!!$*>          S is DOUBLE PRECISION array, dimension (min(M,N))
!!$*>          The singular values of A in decreasing order.
!!$*>          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
!!$*> \endverbatim
!!$*>
!!$*> \param[in] RCOND
!!$*> \verbatim
!!$*>          RCOND is DOUBLE PRECISION
!!$*>          RCOND is used to determine the effective rank of A.
!!$*>          Singular values S(i) <= RCOND*S(1) are treated as zero.
!!$*>          If RCOND < 0, machine precision is used instead.
!!$*> \endverbatim
!!$*>
!!$*> \param[out] RANK
!!$*> \verbatim
!!$*>          RANK is INTEGER
!!$*>          The effective rank of A, i.e., the number of singular values
!!$*>          which are greater than RCOND*S(1).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] WORK
!!$*> \verbatim
!!$*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!$*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LWORK
!!$*> \verbatim
!!$*>          LWORK is INTEGER
!!$*>          The dimension of the array WORK. LWORK must be at least 1.
!!$*>          The exact minimum amount of workspace needed depends on M,
!!$*>          N and NRHS. As long as LWORK is at least
!!$*>              12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,
!!$*>          if M is greater than or equal to N or
!!$*>              12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,
!!$*>          if M is less than N, the code will execute correctly.
!!$*>          SMLSIZ is returned by ILAENV and is equal to the maximum
!!$*>          size of the subproblems at the bottom of the computation
!!$*>          tree (usually about 25), and
!!$*>             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
!!$*>          For good performance, LWORK should generally be larger.
!!$*>
!!$*>          If LWORK = -1, then a workspace query is assumed; the routine
!!$*>          only calculates the optimal size of the WORK array, returns
!!$*>          this value as the first entry of the WORK array, and no error
!!$*>          message related to LWORK is issued by XERBLA.
!!$*> \endverbatim
!!$*>
!!$*> \param[out] IWORK
!!$*> \verbatim
!!$*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!!$*>          LIWORK >= max(1, 3 * MINMN * NLVL + 11 * MINMN),
!!$*>          where MINMN = MIN( M,N ).
!!$*>          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.
!!$*> \endverbatim
!!$*>
!!$*> \param[out] INFO
!!$*> \verbatim
!!$*>          INFO is INTEGER
!!$*>          = 0:  successful exit
!!$*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!!$*>          > 0:  the algorithm for computing the SVD failed to converge;
!!$*>                if INFO = i, i off-diagonal elements of an intermediate
!!$*>                bidiagonal form did not converge to zero.
!!$*> \endverbatim

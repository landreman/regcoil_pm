! Adapted from gs2 by MJL, 20190526.

module regcoil_Gaussian_quadrature
  ! <doc>
  !  Utilities for Gaussian quadrature.
  !  This module provides subroutines to obtain zeros and weights of
  !  Gauss-Legendre and Gauss-Laguerre quadrature rules.
  ! </doc>

  use stel_kinds

  implicit none

  private

  public :: get_legendre_grids

  logical, parameter :: debug=.false.
  logical :: weight_roundoff_correction=.false.
  integer, parameter :: error_unit = 10

contains
  subroutine get_legendre_grids (x1, x2, zero, wgt)
    ! <doc>
    !  returns Legendre zeros and weights in the given interval [x1,x2].
    !  The order is determined from the size of the array 'zero'.
    ! </doc>
    use stel_constants, only: pi
    real(dp), intent (in) :: x1, x2
    real(dp), dimension (:), intent (out) :: zero, wgt
    integer :: i, nn, nh
    double precision :: xold, xnew, pold, pnew
    double precision, dimension (:), allocatable :: zz, ww

    nn = size(zero)
    nh = (nn+1)/2
    allocate (zz(nh))
    allocate (ww(nh))

    ! search zero from 1 to 0 using chebyshev grid points
    ! this is O(nn^2) operations
    xold = cos(pi / (2*(nn+1)))
    pold = legendre_p(nn,xold)
    do i=1, nh
       xnew = cos(pi*(2*i+1)/(2.0*(nn+1)))
       pnew = legendre_p(nn,xnew)
       call find_zero_bisect_newton (nn, xold, xnew, pold, pnew, zz(i))
       xold = xnew
       pold = pnew
    end do
    ! invert them to give zeros in (-1,0]
    zz(1:nn/2) = -zz(1:nn/2)

    ! weight from formula
!    ww = dble(2.0) / (dble(1.0) - zz**2) / legendre_pp(nn,zz,dble(0.0))**2
    ww = dble(2.0) / (dble(1.0) - zz**2) / legendre_pp(nn,zz)**2

    ! rescale (x2 may be smaller than x1)
    zero(1:nh) = real(zz * (x2-x1) + (x1+x2),dp) / 2.0
    zero(nh+1:nn) = real(zz(nn/2:1:-1) * (x1-x2) + (x1+x2),dp) / 2.0
    wgt(1:nh) = real(ww * abs(x2-x1) / 2.0,dp)
    wgt(nh+1:nn) = wgt(nn/2:1:-1)

    deallocate (zz)
    deallocate (ww)

    ! roundoff correction
!!$    if (abs(sum(wgt)/abs(x2-x1)) - 1.0 > epsilon(wgt)) then
!!$       print *, 'roundoff correction occurred'
    if (weight_roundoff_correction) then
       if (mod(nn,2)==0) then
          wgt(nh) = abs(x2-x1) / 2.0 - sum(wgt(1:nh-1))
          wgt(nh+1) = wgt(nh)
       else
          wgt(nh) = abs(x2-x1) - sum(wgt(1:nh-1)) * 2.0
       end if
    end if

    call check_legendre_zero (x1, x2, zero)
    call check_legendre_weights (abs(x2-x1), wgt)

    if (debug) then
       print *, 'get_legendre_grids: sum of weights = ', sum(wgt)
       print *, 'get_legendre_grids: section length = ', abs(x2-x1)
    end if

  end subroutine get_legendre_grids

  subroutine find_zero_bisect_newton (n, xold, xnew, pold, pnew, zz)
    implicit none
    integer, intent (in) :: n
    double precision, intent (in) :: xold, xnew, pold, pnew
!    real, intent (in) :: eps
    double precision, intent (out) :: zz
    integer :: i, maxit=100
    real :: eps
    double precision :: x1, x2, p1, p2, pz

    ! <doc>
    !  eps is declared as real on purpose.
    !  We don't require higher order precision for convergence test
    !  because of a performance reason. The following definition means
    !  eps is a geometric mean of the machine-epsilons in real and double
    !  precisions. 
    !  (note that real/double are promoted to double/double or double/quad 
    !  depending on the compiler.)
    !  
    !  [same applies to eps in find_zero below.]
    ! </doc>

    eps = sqrt(epsilon(eps)*epsilon(x1)) * 2.0
    i=0
    x1 = xold
    x2 = xnew
    p1 = pold
    p2 = pnew

    if (debug) write(*,'(4f10.5)') x1, p1, x2, p2

    ! bisection
    do i=1, 5
       zz = (x1+x2) * dble(.5)
       pz = legendre_p(n,zz)
       if (abs(pz) <= epsilon(pz)) return
       if (pz*p1 < 0.0) then
          p2=pz ; x2=zz
       else
          p1=pz ; x1=zz
       end if
       if (debug) write(*,'(4f10.5)') x1, p1, x2, p2
    end do

    if (debug) print*, 'finished bisection'

    ! newton-raphson
    if (zz==x1) x1 = x2
!    do while (abs(zz/x1-1.0) > eps)
    do i=1, maxit
       x1 = zz
       p1 = legendre_p(n,x1)
!       zz = x1 - p1 / legendre_pp(n,x1,p1)
       zz = x1 - p1 / legendre_pp(n,x1)
       pz = legendre_p(n,zz)
       if (debug) write (*,'(4f10.5)') zz, pz, x1, p1
       if (min(abs(zz/x1-1.0), abs(pz)) < eps) exit
    end do

    if (i==maxit+1) write (error_unit,*) &
         & 'WARNING: too many iterations in get_legendre_grids'

    if (debug) stop "Abort in gauss_quad:find_zero_bisect_newton as debug=.true."

  end subroutine find_zero_bisect_newton

  elemental function legendre_p (n, x)
    implicit none
    integer, intent (in) :: n
    double precision, intent (in) :: x
    integer :: k
    double precision :: p, p1, p2, legendre_p

    select case (n)
    case (0)
       legendre_p = dble(1.0)
    case (1)
       legendre_p = x
    case default
       p1 = x
       p2 = dble(1.0)
       do k=2, n
          p = ((2*k-1)*x*p1 - (k-1)*p2) / k
          p2 = p1
          p1 = p
       end do
       legendre_p = p
    end select

  end function legendre_p

!  elemental function legendre_pp (n, x, p1)
  elemental function legendre_pp (n, x)
    implicit none
    integer, intent (in) :: n
    double precision, intent (in) :: x
!    double precision, intent (in), optional :: p1
    double precision :: legendre_pp

!    if (present(p1)) then
!       legendre_pp = n * ( x * p1 - legendre_p(n-1,x) ) &
!            / (x**2 - dble(1.0))
!    else
    legendre_pp = n * ( x * legendre_p(n,x) - legendre_p(n-1,x) ) &
         / (x**2 - dble(1.0))
!    end if

  end function legendre_pp

  subroutine check_legendre_zero (x0, x1, zero)
    implicit none
    real(dp), intent (in) :: x0, x1
    real(dp), dimension (:), intent (in) :: zero
    logical :: error=.false.
    integer :: nn, nh
    real(dp) :: xx, xmin, xmax
    real(dp), dimension (:), allocatable :: zz

    nn = size(zero)
    nh = (nn+1)/2
    error = .false.
    xmin = min(x0, x1)
    xmax = max(x0, x1)
    allocate (zz(nn))
    zz = zero
    if (zz(1) > zz(nn)) zz(1:nn) = zero(nn:1:-1)
    if (zz(1) < xmin .or. zz(nn) > xmax) then
       write (error_unit,*) 'ERROR in legendre: grid out of range'
       error = .true.
    end if

    if (nn==1) then
       if (abs(2.0*zz(1)/(xmin+xmax) - 1.0) > epsilon(0.0)) then
          write (error_unit, '("ERROR in legendre: zz(1)= ", f20.15)') zz(1)
          error = .true.
       end if
    else
       ! check if distances at the edge: This suffices for nn<=3
       if (zz(2)-zz(1) <= zz(1)-xmin .or. zz(nn)-zz(nn-1) <= xmax-zz(nn)) then
          write (error_unit,*) 'ERROR in legendre: wrong distance at edge'
          error = .true.
       end if
       ! check distances at the center: The above and this suffices for nn=4
       if (mod(nn,2)==0 .and. nn>=4) then
          if ( zz(nh+1)-zz(nh) <= zz(nh)-zz(nh-1) .or. &
               zz(nh+1)-zz(nh) <= zz(nh+2)-zz(nh+1) ) then
             write (error_unit,*) &
                  & 'ERROR in legendre: wrong distance at center'
             error = .true.
          end if
       end if
       if (nn >= 5) then
          ! check if distances are increasing toward center
          ! lower half
          if (any(zz(3:nh)-zz(2:nh-1) <= zz(2:nh-1)-zz(1:nh-2))) then
             write (error_unit,*) 'ERROR in legendre: distance decreasing toward center'
             error = .true.
          end if
          ! upper half
          ! The separate use of nh and nn/2 are intentionally so that they
          ! both work for even and odd cases
          if ( any(zz(nn/2+2:nn-1)-zz(nn/2+1:nn-2) &
               & <= zz(nn/2+3:nn)-zz(nn/2+2:nn-1)) ) then
             write (error_unit,*) 'ERROR in legendre: distance decreasing toward center'
             error = .true.
          end if
       end if

    end if

    ! check if legendre_p(n, zero(i)) are close enough to zero
    if (debug) then
       xx = maxval(abs(real( legendre_p(nn, &
            dble((zz(:)-xmin)/(xmax-xmin)*2.0-1.0)), dp )))
       if (xx/nn**2 > epsilon(xx)) then
          write (error_unit,*) 'WARNING in legendre: maxval(n,zz(:))= ', xx
          ! Do not stop as it is mostly a minor issue
       end if
    end if

    if (error) stop 'STOP in check_legendre_zero'

  end subroutine check_legendre_zero

!  subroutine check_legendre_weights (norm, wgt, eps)
  subroutine check_legendre_weights (norm, wgt)
    implicit none
    real(dp), intent (in) :: norm!, eps
    real(dp), dimension (:), intent (in) :: wgt
    logical :: error=.false.
    integer :: n, nh
    real(dp) :: s

    n = size(wgt)
    error = .false.
    nh = (n+1)/2

    ! check if weights are all positive
    if (any(wgt < 0.)) then
       write (error_unit,*) 'ERROR in legendre: weights got negative'
       error = .true.
    end if

    if (n>=2) then
       ! check symmetry of weights
       if ( any( abs(wgt(n:n+1-n/2:-1)/wgt(1:n/2) - 1.0) > epsilon(wgt) ) ) then
          write (error_unit,*) 'WARNING in legendre: symmetry of weights broken'
          error = .true.
       end if

       ! check if weights are increasing toward center
       if (n>=3) then
          if (any(wgt(2:nh) <= wgt(1:nh-1))) then
             write (error_unit,*) 'ERROR in legendre: weights decreasing toward center'
             error = .true.
          end if
       end if
    end if

    ! check if their sum is close enough to normalized value
    if (debug) then
!       s = sum(wgt)
       s = sum(dble(wgt)) ! ignore roundoff error arising from 8-byte summation
       if (abs(norm/s-1.0) > epsilon(s)) then
          write (error_unit,*) 'WARNING in legendre: weights summation incorrect:', &
               & size(wgt), s/norm-1.0
          ! Do not stop as it is mostly a minor issue
       end if
    end if

    if (error) stop 'STOP in check_legendre_weights'

  end subroutine check_legendre_weights

end module regcoil_Gaussian_quadrature

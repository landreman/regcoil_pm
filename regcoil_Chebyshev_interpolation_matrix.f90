subroutine regcoil_Chebyshev_interpolation_matrix(N, M, xk, x, matrix)
  ! Interpolates from Chebyshev points xk to a new grid x.
  ! Based on the Matlab function chebint included with DMSuite.
  !
  ! Inputs:
  ! N = number of Chebyshev points.
  ! M = number of points on the new grid.
  ! xk(N) = Chebyshev points.  (This subroutine will not check that they are indeed Chebyshev points!!)
  ! x(M) = grid points onto which we will interpolate.
  !
  ! Outputs:
  ! matrix(M,N) = interpolation matrix

  use stel_kinds

  implicit none

  integer, intent(in) :: N, M
  real(dp), intent(in) :: xk(N), x(M)
  real(dp), intent(out) :: matrix(M,N)

  integer :: i, j
  real(dp), allocatable :: w(:), D(:,:), denominator(:), xScaled(:), xkScaled(:)
  real(dp) :: temp, xMin, xMax, xMid

  if (N == 1) then
     matrix = 1.0_dp
     return
  end if

  if (N<1) then
     print *,"Error in regcoil_Chebyshev_interpolation_matrix.f90! N must be at least 1."
     stop
  end if

  if (M<1) then
     print *,"Error! M must be at least 1."
     stop
  end if

  allocate(xScaled(M))
  allocate(xkScaled(N))
  xMin = minval(xk)
  xMax = maxval(xk)
  xMid = (xMin+xMax)*(0.5d+0)
  xScaled = -2*(x-xMid)/(xMax-xMin)
  xkScaled = -2*(xk-xMid)/(xMax-xMin)

!!$  print *,"xScaled:",xScaled
!!$  print *,"xkScaled:",xkScaled

  allocate(w(N))
  do i=1,N
     w(i) = (-1) ** (i-1)
  end do
  w(1) = w(1) / (2.0d+0)
  w(N) = w(N) / (2.0d+0)
  
  allocate(D(M,N))
  do i=1,M
     do j=1,N
        temp = xScaled(i) - xkScaled(j)
        if (abs(temp)<1.0d-15) then
           D(i,j) = 1/(temp+1.0d-15)
        else
           D(i,j) = 1/temp
        end if
     end do
  end do

  allocate(denominator(M))
  denominator = matmul(D,w)

!!$  print *,"w:",w
!!$  print *,"D:"
!!$  do i=1,M
!!$     print *,D(i,:)
!!$  end do
!!$  print *,"Denominator:",denominator

  do i=1,M
     do j=1,N
        matrix(i,j) = D(i,j) * w(j) / denominator(i)
     end do
  end do

  deallocate(w, D, denominator, xScaled, xkScaled)


end subroutine Regcoil_Chebyshev_interpolation_matrix


!-------------------------------------------------------------------------------
! algebra.f90
!
! Contains subroutines for basic linear algraic functions used in other modules
!-------------------------------------------------------------------------------
module algebra

use magpie_globals, only: dp

implicit none

contains

!-------------------------------------------------------------------------------
! inverse_3x3(m11, m12, m13, m21, m22, m23, m31, m32, m33, 
!             i11, i12, i13, i21, i22, i23, i31, i32, i33)
!
! Computes the inverse of a 3x3 matrix.
! Input variables are the components m_ij of the matrix to invert.
! Output variables are the components i_ij of the inverse.
!
! Source: Wolfram MathWorld
!-------------------------------------------------------------------------------
subroutine inverse_3x3(m11, m12, m13, m21, m22, m23, m31, m32, m33, &
                       i11, i12, i13, i21, i22, i23, i31, i32, i33)

    implicit none

    real(dp), intent(IN)  :: m11, m12, m13, m21, m22, m23, m31, m32, m33
    real(dp), intent(OUT) :: i11, i12, i13, i21, i22, i23, i31, i32, i33
    real(dp) :: det

    det = m11*(m22*m33-m23*m32) - m12*(m21*m33-m23*m31) + m13*(m21*m32-m22*m31)

    i11 = (m22*m33 - m23*m32) / det
    i12 = (m13*m32 - m12*m33) / det
    i13 = (m12*m23 - m13*m22) / det
    i21 = (m23*m31 - m21*m33) / det
    i22 = (m11*m33 - m13*m31) / det
    i23 = (m13*m21 - m11*m23) / det
    i31 = (m21*m32 - m22*m31) / det
    i32 = (m12*m31 - m11*m32) / det
    i33 = (m11*m22 - m12*m21) / det

end subroutine inverse_3x3

!-------------------------------------------------------------------------------
! inverse_2x2(m11, m12, m21, m22, i11, i12, i21, i22)
!
! Computes the inverse of a 2x2 matrix.
! Input variables are the components m_ij of the matrix to invert.
! Output variables are the components i_ij of the inverse.
!
! Source: Wolfram MathWorld
!-------------------------------------------------------------------------------
subroutine inverse_2x2(m11, m12, m21, m22, i11, i12, i21, i22)

    implicit none

    real(dp), intent(IN)  :: m11, m12, m21, m22
    real(dp), intent(OUT) :: i11, i12, i21, i22
    real(dp) :: det

    det = m11*m22 - m12*m21

    i11 =  m22 / det
    i12 = -m12 / det
    i21 = -m21 / det
    i22 =  m11 / det

end subroutine inverse_2x2

!-------------------------------------------------------------------------------
! product_3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33, 
!             b11, b12, b13, b21, b22, b23, b31, b32, b33,
!             p11, p12, p13, p21, p22, p23, p31, p32, p33)
!
! Computes the product p of two 3x3 matrices a and b
!-------------------------------------------------------------------------------
subroutine product_3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33, &
                       b11, b12, b13, b21, b22, b23, b31, b32, b33, &
                       p11, p12, p13, p21, p22, p23, p31, p32, p33)

    implicit none

    real(dp), intent(IN)  :: a11, a12, a13, a21, a22, a23, a31, a32, a33
    real(dp), intent(IN)  :: b11, b12, b13, b21, b22, b23, b31, b32, b33
    real(dp), intent(OUT) :: p11, p12, p13, p21, p22, p23, p31, p32, p33

    p11 = a11*b11 + a12*b21 + a13*b31
    p12 = a11*b12 + a12*b22 + a13*b32
    p13 = a11*b13 + a12*b23 + a13*b33
    p21 = a21*b11 + a22*b21 + a23*b31
    p22 = a21*b12 + a22*b22 + a23*b32
    p23 = a21*b13 + a22*b23 + a23*b33
    p31 = a31*b11 + a32*b21 + a33*b31
    p32 = a31*b12 + a32*b22 + a33*b32
    p33 = a31*b13 + a32*b23 + a33*b33

end subroutine product_3x3

!-------------------------------------------------------------------------------
! linear_interpolate(n_data, x_data, y_data, n_query, x_query, y_interp)
!
! Basic 1d linear interpolator. 
! 
! *** Note: x_data are assumed to be in ascending order!!
!
! Input parameters:
!     integer :: n_data           -> number of (x,y) pairs in the input data
!     real(dp) :: x_data, y_data  -> x and y values of the input data
!     integer :: n_query          -> number of query points
!
! Output parameters:
!     real(dp) :: y_interp        -> interpolation of y at the query points
!-------------------------------------------------------------------------------
subroutine linear_interpolate(n_data, x_data, y_data, n_query, x_query, &
                              y_interp)

    implicit none

    integer, intent(IN)  :: n_data, n_query
    real(dp), intent(IN)  :: x_data(:), y_data(:), x_query(:)
    real(dp), intent(OUT) :: y_interp(:)
    integer :: i, j
    real(dp) :: dydx

    do i = 1, n_query
        do j = 2, n_data
            if (x_query(i) < x_data(j) .or. j == n_data) then
                dydx = (y_data(j)-y_data(j-1))/(x_data(j)-x_data(j-1))
                y_interp(i) = y_data(j-1) + dydx*(x_query(i) - x_data(j-1))
                exit
            end if
        end do
    end do

end subroutine linear_interpolate

end module algebra

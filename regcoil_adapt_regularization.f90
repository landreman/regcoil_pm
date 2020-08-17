subroutine regcoil_adapt_regularization(isaved)

  use regcoil_variables
  use stel_constants
  use magpie_globals, only: stell_symm

  implicit none

  integer, intent(IN) :: isaved
  real(dp) :: moment_rel, reg_matrix_factor
  integer :: izeta, j_rZetaZ, j, offset

  print "(a)", "    Updating the regularization factors"

  ! Calculate adaptation factors based on overshoots/undershoots of max moment
  do izeta = 1, nzeta_coil

     moment_rel = abs_M(1,izeta,1,isaved)*qhex_arr_base(izeta)%vol &
                     / qhex_max_moment(izeta)
     if (moment_rel > 1.0d0 .or. adaptation_factor(izeta) > 1.0d0) &
        adaptation_factor(izeta) = &
           adaptation_factor(izeta) &
              * (1.d0 + adaptation_alpha*(moment_rel**2 - 1.d0))

  end do

  ! Update the regularization matrix
  if (stell_symm) then
      reg_matrix_factor = 2.0d0
  else
      reg_matrix_factor = 1.0d0
  end if

  do j_RZetaZ = 1, 3

     offset = (j_RZetaZ-1)*nzeta_coil

     do j = 1, nzeta_coil
        select case (regularization_type)
        case (1)
           matrix_regularization(offset+j,offset+j) = &
              reg_matrix_factor * adaptation_factor(j) * qhex_arr(j)%vol &
                 * d(1,j) ** regularization_d_exponent
        case (2)
           matrix_regularization(offset+j,offset+j) = &
              reg_matrix_factor * adaptation_factor(j) * nzeta_coil_inv &
                 * ( qhex_arr(j)%vol / qhex_max_moment(j) ) ** 2
        end select
     end do

  end do

end subroutine regcoil_adapt_regularization



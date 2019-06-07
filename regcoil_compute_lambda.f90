subroutine regcoil_compute_lambda()

  use regcoil_variables, only: nlambda, lambda_min, lambda_max, lambda, lambda_option, lambda_option_single, lambda_single, lambda_option_scan, verbose, nd
  use stel_kinds

  implicit none

  integer :: j

  if (allocated(lambda)) deallocate(lambda)

  if (trim(lambda_option) == lambda_option_single) then
     nlambda = 1
     allocate(lambda(nd))
     lambda = lambda_single
     if (verbose) print "(a,es10.3)", " Using a single value of lambda: ",lambda_single
     return
  end if

  allocate(lambda(nlambda))
  
  lambda(1) = 0
  do j = 1,nlambda-1
     lambda(j+1) = lambda_min * exp((log(lambda_max/lambda_min)*(j-1))/(nlambda-2))
  end do

  if (trim(lambda_option)==lambda_option_scan .and. verbose) then
     print *,"We will use the following values of the regularization weight lambda:"
     print "(*(es10.3))",lambda
  end if

end subroutine regcoil_compute_lambda

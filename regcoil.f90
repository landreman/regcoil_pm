! Main program

program regcoil

  use regcoil_variables
  use regcoil_init_plasma_mod

  implicit none

  integer :: tic, toc, countrate

  print *,"This is REGCOIL_PM,"
  print *,"a regularized least-squares method for designing permanent magnets for stellarators."
  call system_clock(tic,countrate)

  call regcoil_read_input()
  call regcoil_validate_input()
  call regcoil_compute_lambda()

  ! Define the position vector and normal vector at each grid point for the surfaces:
  call regcoil_init_plasma()
  call regcoil_init_coil_surface()

  ! Initialize some of the vectors and matrices needed:
  call regcoil_read_bnorm()
  call regcoil_build_matrices()
  call regcoil_prepare_solve()

  select case (trim(lambda_option))
  case (lambda_option_single,lambda_option_scan)
     call regcoil_lambda_scan()
  case (lambda_option_search)
     call regcoil_auto_regularization_solve()
  case default
     print *,"Invalid lambda_option:",lambda_option
     stop
  end select

  call system_clock(toc)
  total_time = real(toc-tic)/countrate

  call regcoil_write_output()

  if (write_mgrid) call regcoil_write_mgrid()
 
  print *,"REGCOIL_PM complete. Total time=",total_time,"sec."
  print *,"You can run regcoilPlot ",trim(output_filename)," to plot results."

end program regcoil

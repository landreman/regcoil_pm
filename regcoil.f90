! Main program

program regcoil

  use regcoil_variables
  use regcoil_init_plasma_mod

  implicit none

  integer :: tic, toc, countrate

  print "(a)","This is REGCOIL_PM,"
  print "(a)","a regularized least-squares method for designing permanent magnets for stellarators."
  call system_clock(tic,countrate)

  call regcoil_read_input()
  call regcoil_validate_input()
  call regcoil_compute_lambda()

  ! Define the position vector and normal vector at each grid point for the surfaces:
  call regcoil_init_plasma()

  ! Define the magnets
  select case (trim(magnet_type))
  case ('continuous')
     call regcoil_init_coil_surface()
  case ('qhex')
     call magpie_init_qhex_set()
  case default
     print *, "Invalid magnet_type option: ", magnet_type
  end select

  ! Initialize some of the vectors and matrices needed:
  call regcoil_init_ports()
  call regcoil_read_bnorm()

  select case (trim(magnet_type))
  case ('continuous') 
     call regcoil_init_basis_functions()
  case ('qhex')
     num_basis_functions = nzeta_coil
     system_size = 3 * num_basis_functions
  end select

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

  call regcoil_evaluate_outer_surface()
  call regcoil_write_output()
  if (trim(magnet_type) == 'qhex') call magpie_write_output()

  if (write_mgrid) call regcoil_write_mgrid()
 
  print *,"REGCOIL_PM complete. Total time=",total_time,"sec."
  print *,"You can run regcoilPlot ",trim(output_filename)," to plot results."

end program regcoil

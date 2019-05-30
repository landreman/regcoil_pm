subroutine regcoil_write_output

  use regcoil_variables
  use ezcdf

  implicit none

  integer :: ierr, ncid

  ! Same convention as in VMEC:
  ! Prefix vn_ indicates the variable name used in the .nc file.

  ! Scalars:
  character(len=*), parameter :: &
       vn_lambda_option = "lambda_option", &
       vn_nfp = "nfp", &
       vn_geometry_option_plasma = "geometry_option_plasma", &
       vn_geometry_option_coil = "geometry_option_coil", &
       vn_ntheta_plasma = "ntheta_plasma", &
       vn_nzeta_plasma = "nzeta_plasma", &
       vn_nzetal_plasma = "nzetal_plasma", &
       vn_ntheta_coil = "ntheta_coil", &
       vn_nzeta_coil = "nzeta_coil", &
       vn_nzetal_coil = "nzetal_coil", &
       vn_ns_magnetization = "ns_magnetization", &
       vn_ns_integration = "ns_integration", &
       vn_system_size = "system_size", &
       vn_a_plasma  = "a_plasma", &
       vn_a_coil  = "a_coil", &
       vn_R0_plasma  = "R0_plasma", &
       vn_R0_coil  = "R0_coil", &
       vn_mpol_magnetization = "mpol_magnetization", &
       vn_ntor_magnetization = "ntor_magnetization", &
       vn_mnmax_magnetization = "mnmax_magnetization", &
       vn_mnmax_plasma = "mnmax_plasma", &
       vn_mnmax_coil = "mnmax_coil", &
       vn_num_basis_functions = "num_basis_functions", &
       vn_symmetry_option = "symmetry_option", &
       vn_area_plasma = "area_plasma", &
       vn_area_coil = "area_coil", &
       vn_volume_plasma = "volume_plasma", &
       vn_volume_coil = "volume_coil", &
       vn_curpol = "curpol", &
       vn_nlambda = "nlambda", &
       vn_nsaved = "nsaved", &
       vn_nd = "nd", &
       vn_total_time = "total_time", &
       vn_exit_code = "exit_code", &
       vn_chi2_B_target = "chi2_B_target", &
       vn_sign_normal = "sign_normal", &
       vn_d_initial = "d_initial"

  ! Arrays with dimension 1
  character(len=*), parameter :: &
       vn_theta_plasma = "theta_plasma", &
       vn_zeta_plasma = "zeta_plasma", &
       vn_zetal_plasma = "zetal_plasma", &
       vn_theta_coil = "theta_coil", &
       vn_zeta_coil = "zeta_coil", &
       vn_zetal_coil = "zetal_coil", &
       vn_s_magnetization = "s_magnetization", &
       vn_s_integration = "s_integration", &
       vn_s_weights = "s_weights", &
       vn_xm_magnetization = "xm_magnetization", &
       vn_xn_magnetization = "xn_magnetization", &
       vn_xm_plasma = "xm_plasma", &
       vn_xn_plasma = "xn_plasma", &
       vn_xm_coil = "xm_coil", &
       vn_xn_coil = "xn_coil", &
       vn_rmnc_plasma = "rmnc_plasma", &
       vn_rmns_plasma = "rmns_plasma", &
       vn_zmnc_plasma = "zmnc_plasma", &
       vn_zmns_plasma = "zmns_plasma", &
       vn_rmnc_coil = "rmnc_coil", &
       vn_rmns_coil = "rmns_coil", &
       vn_zmnc_coil = "zmnc_coil", &
       vn_zmns_coil = "zmns_coil", &
       vn_RHS_B = "RHS_B", &
       vn_lambda = "lambda", &
       vn_chi2_B = "chi2_B", &
       vn_chi2_M = "chi2_M", &
       vn_max_Bnormal = "max_Bnormal", &
       vn_max_M = "max_M", &
       vn_min_M = "min_M"

  ! Arrays with dimension 2
  character(len=*), parameter :: &
       vn_norm_normal_plasma  = "norm_normal_plasma", &
       vn_norm_normal_coil  = "norm_normal_coil", &
       vn_Bnormal_from_TF_and_plasma_current = "Bnormal_from_TF_and_plasma_current", &
       vn_mean_curvature_coil = "mean_curvature_coil", &
       vn_matrix_B = "matrix_B", &
       vn_matrix_regularization = "matrix_regularization", &
       vn_interpolate_magnetization_to_integration = "interpolate_magnetization_to_integration"

  ! Arrays with dimension 3
  character(len=*), parameter :: &
       vn_r_plasma  = "r_plasma", &
       vn_r_coil  = "r_coil", &
       vn_drdtheta_plasma  = "drdtheta_plasma", &
       vn_drdtheta_coil  = "drdtheta_coil", &
       vn_drdzeta_plasma  = "drdzeta_plasma", &
       vn_drdzeta_coil  = "drdzeta_coil", &
       vn_normal_plasma = "normal_plasma", &
       vn_normal_coil = "normal_coil", &
       vn_Bnormal_total = "Bnormal_total", &
       vn_Jacobian_coil = "Jacobian_coil"

  ! Arrays with dimension 4
  character(len=*), parameter :: &
       vn_g  = "g", &
       vn_abs_M = "abs_M", &
       vn_magnetization_vector_mn = "magnetization_vector_mn"

  ! Arrays with dimension 5
  character(len=*), parameter :: &
       vn_magnetization_vector  = "magnetization_vector"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now create variables that name the dimensions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Arrays with dimension 1:
  character(len=*), parameter, dimension(1) :: &
       ntheta_plasma_dim = (/'ntheta_plasma'/), &
       nzeta_plasma_dim = (/'nzeta_plasma'/), &
       nzetal_plasma_dim = (/'nzetal_plasma'/), &
       ntheta_coil_dim = (/'ntheta_coil'/), &
       nzeta_coil_dim = (/'nzeta_coil'/), &
       nzetal_coil_dim = (/'nzetal_coil'/), &
       mnmax_magnetization_dim = (/'mnmax_magnetization'/), &
       mnmax_plasma_dim = (/'mnmax_plasma'/), &
       mnmax_coil_dim = (/'mnmax_coil'/), &
       nthetanzeta_plasma_dim = (/'ntheta_nzeta_plasma'/), &
       num_basis_functions_dim = (/'num_basis_functions'/), &
       nlambda_dim = (/'nlambda'/), &
       nsaved_dim = (/'nsaved'/), &
       system_size_dim = (/'system_size'/), &
       ns_magnetization_dim = (/'ns_magnetization'/), &
       ns_integration_dim = (/'ns_integration'/)

  ! Arrays with dimension 2:
  ! The form of the array declarations here is inspired by
  ! http://stackoverflow.com/questions/21552430/gfortran-does-not-allow-character-arrays-with-varying-component-lengths
  character(len=*), parameter, dimension(2) :: &
       ntheta_nzeta_plasma_dim = (/ character(len=50) :: 'ntheta_plasma','nzeta_plasma'/), &
       ntheta_nzeta_coil_dim = (/ character(len=50) :: 'ntheta_coil','nzeta_coil'/), &
       nthetanzeta_plasma_nthetanzeta_coil_dim = (/ character(len=50) :: 'ntheta_nzeta_plasma','ntheta_nzeta_coil'/), &
       nthetanzeta_plasma_basis_dim = (/ character(len=50) :: 'ntheta_nzeta_plasma','num_basis_functions'/), &
       basis_basis_dim = (/ character(len=50) :: 'num_basis_functions','num_basis_functions'/), &
       basis_nsaved_dim = (/ character(len=50) :: 'num_basis_functions','nsaved'/), &
       system_size_system_size_dim = (/ character(len=50) :: 'system_size','system_size'/), &
       ns_integration_ns_magnetization_dim = (/ character(len=50) :: 'ns_integration', 'ns_magnetization' /)

  ! Arrays with dimension 3:
  character(len=*), parameter, dimension(3) :: &
       xyz_ntheta_nzetal_plasma_dim = (/ character(len=50) :: 'xyz','ntheta_plasma','nzetal_plasma'/), &
       xyz_ntheta_nzetal_coil_dim = (/ character(len=50) :: 'xyz','ntheta_coil','nzetal_coil'/), &
       ntheta_nzeta_coil_nsaved_dim = (/ character(len=50) :: 'ntheta_coil','nzeta_coil','nsaved'/), &
       ntheta_nzeta_plasma_nsaved_dim = (/ character(len=50) :: 'ntheta_plasma','nzeta_plasma','nsaved'/), &
       ntheta_nzeta_coil_ns_integration_dim = (/ character(len=50) :: 'ntheta_coil','nzeta_coil','ns_integration'/)

  ! Arrays with dimension 4:
  character(len=*), parameter, dimension(4) :: &
       ntheta_nzeta_plasma_basis_ns_magnetization_RZetaZ_dim = (/ character(len=50) :: 'ntheta_nzeta_plasma','num_basis_functions','ns_magnetization','RZetaZ'/), &
       ntheta_nzeta_coil_ns_magnetization_nsaved_dim = (/ character(len=50) :: 'ntheta_coil','nzeta_coil', 'ns_magnetization','nsaved' /), &
       basis_ns_magnetization_RZetaZ_nsaved_dim = (/ character(len=50) :: 'num_basis_functions','ns_magnetization','RZetaZ','nsaved' /) 

  ! Arrays with dimension 5:
  character(len=*), parameter, dimension(5) :: &
       ntheta_nzeta_coil_ns_magnetization_RZetaZ_nsaved_dim = (/ character(len=50) :: 'ntheta_coil','nzeta_coil', 'ns_magnetization','RZetaZ','nsaved' /)

  character(len=*), parameter :: input_parameter_text = ' See the user manual documentation for the input parameter of the same name.'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call cdf_open(ncid,output_filename,'w',ierr)
  IF (ierr .ne. 0) then
     print *,"Error opening output file ",output_filename
     stop
  end IF

  ! Scalars

  call cdf_define(ncid, vn_lambda_option, lambda_option)
  call cdf_setatt(ncid, vn_lambda_option, "Option determining how the values of the regularization parameter lambda were selected. " // input_parameter_text)

  call cdf_define(ncid, vn_nfp, nfp)
  call cdf_setatt(ncid, vn_nfp, 'Number of field periods, i.e. the number of identical toroidal segments, 5 for W7-X, 4 for HSX, etc. ' // &
       'Equivalent to the VMEC variable of the same name.' // input_parameter_text)

  call cdf_define(ncid, vn_geometry_option_plasma, geometry_option_plasma)
  call cdf_setatt(ncid, vn_geometry_option_plasma, 'Method used to define the geometry of the plasma surface.' // input_parameter_text)

  call cdf_define(ncid, vn_geometry_option_coil, geometry_option_coil)
  call cdf_setatt(ncid, vn_geometry_option_coil, 'Method used to define the geometry of the inner magnetization surface.' // input_parameter_text)

  call cdf_define(ncid, vn_ntheta_plasma, ntheta_plasma)
  call cdf_setatt(ncid, vn_ntheta_plasma, 'Number of grid points used in the poloidal angle theta on the plasma surface.' // input_parameter_text)

  call cdf_define(ncid, vn_nzeta_plasma, nzeta_plasma)
  call cdf_setatt(ncid, vn_nzeta_plasma, 'Number of grid points used in the toridal angle zeta per identical toroidal period, on the plasma surface.' // input_parameter_text)

  call cdf_define(ncid, vn_nzetal_plasma, nzetal_plasma)
  call cdf_setatt(ncid, vn_nzetal_plasma, 'Number of grid points used in the toridal angle, including all the nfp identical toroidal periods, on the plasma surface.')

  call cdf_define(ncid, vn_ntheta_coil, ntheta_coil)
  call cdf_setatt(ncid, vn_ntheta_coil, 'Number of grid points used in the poloidal angle theta on the inner magnetization surface.' // input_parameter_text)

  call cdf_define(ncid, vn_nzeta_coil, nzeta_coil)
  call cdf_setatt(ncid, vn_nzeta_coil, 'Number of grid points used in the toridal angle zeta per identical toroidal period, on the inner magnetization surface.' // input_parameter_text)

  call cdf_define(ncid, vn_nzetal_coil, nzetal_coil)
  call cdf_setatt(ncid, vn_nzetal_coil, 'Number of grid points used in the toridal angle, including all the nfp identical toroidal periods, on the inner magnetization surface.')

  call cdf_define(ncid, vn_ns_magnetization, ns_magnetization)
  call cdf_setatt(ncid, vn_ns_magnetization, 'Number of grid points in the radial (s) direction used to discretize the magnetization.')

  call cdf_define(ncid, vn_ns_integration, ns_integration)
  call cdf_setatt(ncid, vn_ns_integration, 'Number of grid points in the radial (s) direction used for integration over the magnetization volume.')

  call cdf_define(ncid, vn_system_size, system_size)
  call cdf_setatt(ncid, vn_system_size, 'Size of the right-hand-side and solution vector for the normal equations, i.e. the number of discrete degrees of freedom in the magnetization.')

  call cdf_define(ncid, vn_a_plasma, a_plasma)
  call cdf_define(ncid, vn_a_coil, a_coil)
  call cdf_define(ncid, vn_R0_plasma, R0_plasma)
  call cdf_define(ncid, vn_R0_coil, R0_coil)

  call cdf_define(ncid, vn_mpol_magnetization, mpol_magnetization)
  call cdf_setatt(ncid, vn_mpol_magnetization, 'The maximum poloidal mode number retained in the magnetization.' // input_parameter_text)

  call cdf_define(ncid, vn_ntor_magnetization, ntor_magnetization)
  call cdf_setatt(ncid, vn_ntor_magnetization, 'ntor_magnetization * nfp is the maximum toroidal mode number retained in the magnetization.' // input_parameter_text)

  call cdf_define(ncid, vn_mnmax_magnetization, mnmax_magnetization)
  call cdf_setatt(ncid, vn_mnmax_magnetization, 'Number of unique (m,n) pairs for the Fourier modes of the magnetization. ' // &
       'Equal to mpol_magnetization*(ntor_magnetization*2+1) + ntor_magnetization.')

  call cdf_define(ncid, vn_mnmax_plasma, mnmax_plasma)
  call cdf_setatt(ncid, vn_mnmax_plasma, 'Number of unique (m,n) pairs for the Fourier modes retained in the plasma surface.')

  call cdf_define(ncid, vn_mnmax_coil, mnmax_coil)
  call cdf_setatt(ncid, vn_mnmax_coil, 'Number of unique (m,n) pairs for the Fourier modes retained in the inner magnetization surface.')

  call cdf_define(ncid, vn_num_basis_functions, num_basis_functions)
  call cdf_setatt(ncid, vn_num_basis_functions, 'Number of cos(m*theta-n*zeta) and/or sin(m*theta-n*zeta) Fourier modes ' // &
       'of the magnetization. Equal to mnmax_magnetization * 1 or 2, depending on symmetry_option.')

  call cdf_define(ncid, vn_symmetry_option, symmetry_option)
  call cdf_setatt(ncid, vn_symmetry_option, '1 = The magnetization was forced to be stellarator-symmetric. ' // &
       '2 = The magnetization was forced to be stellarator-antisymmetric. ' // &
       '3 = The magnetization was allowed to have both stellarator-symmetric and antisymmetric components')

  call cdf_define(ncid, vn_area_plasma, area_plasma)
  call cdf_setatt(ncid, vn_area_plasma, 'Area of the plasma surface in meters^2')

  call cdf_define(ncid, vn_area_coil, area_coil)
  call cdf_setatt(ncid, vn_area_coil, 'Area of the inner magnetization surface in meters^2')

  call cdf_define(ncid, vn_volume_plasma, volume_plasma)
  call cdf_setatt(ncid, vn_volume_plasma, 'Volume of the plasma surface in meters^3')

  call cdf_define(ncid, vn_volume_coil, volume_coil)
  call cdf_setatt(ncid, vn_volume_coil, 'Volume of the inner magnetization surface in meters^3. Note that this quantity is NOT the volume of the magnetization region!')

  call cdf_define(ncid, vn_nlambda, nlambda)
  call cdf_setatt(ncid, vn_nlambda, 'Number of values of the regularization parameter lambda examined.')

  call cdf_define(ncid, vn_nd, nd)
  call cdf_setatt(ncid, vn_nd, 'Number of iterations performed to determine the thickness d.')

  call cdf_define(ncid, vn_nsaved, nsaved)
  call cdf_setatt(ncid, vn_nsaved, 'Number of configurations saved in the output file. This value equals nd if lambda_option=single, otherwise nsaved=nlambda.')

  call cdf_define(ncid, vn_total_time, total_time)
  call cdf_setatt(ncid, vn_total_time, 'Total time it took regcoil to run, in seconds.')

  call cdf_define(ncid, vn_exit_code, exit_code)
  call cdf_setatt(ncid, vn_exit_code, "Only meaningful when lambda_option=lambda_option_search so a lambda search is performed. " // &
       "exit_code = 0 means the lambda search was successful. " // &
       "exit_code = -1 means the lambda search did not converge to the requested tolerance within nlambda iterations. " // &
       "exit_code = -2 means the current_density_target you have set is not achievable because it is too low. " // &
       "exit_code = -3 means the current_density_target you have set is not achievable because it is too high.")

  if (trim(lambda_option)==lambda_option_search) then
     call cdf_define(ncid, vn_chi2_B_target, chi2_B_target)
     call cdf_setatt(ncid, vn_chi2_B_target, 'The value of chi^2_B at the final value of regularization parameter lambda resulting from the lambda search. ' // &
          'Units = Tesla^2 meters^2.')
  end if

  call cdf_define(ncid, vn_sign_normal, sign_normal)
  call cdf_setatt(ncid, vn_sign_normal, '+1 if the magnetization region is obtained by moving in the normal_coil direction from the boundary surface, or ' // &
       '-1 if the magnetization region is obtained by moving opposite to the normal_coil direction from the boundary surface.')

  call cdf_define(ncid, vn_d_initial, d_initial)
  call cdf_setatt(ncid, vn_d_initial, 'Initial thickness of the magnetization region, in meters. ' // input_parameter_text)

  ! Arrays with dimension 1

  call cdf_define(ncid, vn_theta_plasma, theta_plasma, dimname=ntheta_plasma_dim)
  call cdf_setatt(ncid, vn_theta_plasma, 'Grid points of the poloidal angle on the plasma surface.')

  call cdf_define(ncid, vn_zeta_plasma, zeta_plasma, dimname=nzeta_plasma_dim)
  call cdf_setatt(ncid, vn_theta_plasma, 'Grid points of the toroidal angle on the plasma surface. ' // &
       'Only the first of the nfp identical toroidal periods is included')

  call cdf_define(ncid, vn_zetal_plasma, zetal_plasma, dimname=nzetal_plasma_dim)
  call cdf_setatt(ncid, vn_theta_plasma, 'Grid points of the toroidal angle on the plasma surface, including all nfp toroidal periods.')

  call cdf_define(ncid, vn_theta_coil, theta_coil, dimname=ntheta_coil_dim)
  call cdf_setatt(ncid, vn_theta_coil, 'Grid points of the poloidal angle on the inner magnetization surface.')

  call cdf_define(ncid, vn_zeta_coil, zeta_coil, dimname=nzeta_coil_dim)
  call cdf_setatt(ncid, vn_theta_coil, 'Grid points of the toroidal angle on the inner magnetization surface. ' // &
       'Only the first of the nfp identical toroidal periods is included')

  call cdf_define(ncid, vn_zetal_coil, zetal_coil, dimname=nzetal_coil_dim)
  call cdf_setatt(ncid, vn_theta_coil, 'Grid points of the toroidal angle on the inner magnetization surface, including all nfp toroidal periods.')

  call cdf_define(ncid, vn_s_magnetization, s_magnetization, dimname=ns_magnetization_dim)
  call cdf_setatt(ncid, vn_s_magnetization, 'Grid points in the radial (s) direction in the discrete representation of the magnetization.')

  call cdf_define(ncid, vn_s_integration, s_integration, dimname=ns_integration_dim)
  call cdf_setatt(ncid, vn_s_integration, 'Grid points in the radial (s) direction used for integration over the magnetization region.')

  call cdf_define(ncid, vn_s_weights, s_weights, dimname=ns_integration_dim)
  call cdf_setatt(ncid, vn_s_weights, 'Integration weights used for integration over the radial (s) direction of the magnetization region.')

  call cdf_define(ncid, vn_xm_magnetization, xm_magnetization, dimname=mnmax_magnetization_dim)
  call cdf_setatt(ncid, vn_xm_magnetization, 'Values of poloidal mode number m used in the Fourier representation of the magnetization.')
  call cdf_define(ncid, vn_xn_magnetization, xn_magnetization, dimname=mnmax_magnetization_dim)
  call cdf_setatt(ncid, vn_xm_magnetization, 'Values of toroidal mode number n used in the Fourier representation of the magnetization.')

  call cdf_define(ncid, vn_xm_plasma, xm_plasma, dimname=mnmax_plasma_dim)
  call cdf_setatt(ncid, vn_xm_plasma, 'Values of poloidal mode number m used in the Fourier representation of the plasma surface.')
  call cdf_define(ncid, vn_xn_plasma, xn_plasma, dimname=mnmax_plasma_dim)
  call cdf_setatt(ncid, vn_xm_plasma, 'Values of toroidal mode number n used in the Fourier representation of the plasma surface.')

  call cdf_define(ncid, vn_xm_coil, xm_coil, dimname=mnmax_coil_dim)
  call cdf_setatt(ncid, vn_xm_coil, 'Values of poloidal mode number m used in the Fourier representation of the inner magnetization surface.')
  call cdf_define(ncid, vn_xn_coil, xn_coil, dimname=mnmax_coil_dim)
  call cdf_setatt(ncid, vn_xm_coil, 'Values of toroidal mode number n used in the Fourier representation of the inner magnetization surface.')

  call cdf_define(ncid, vn_rmnc_plasma, rmnc_plasma, dimname=mnmax_plasma_dim)
  call cdf_setatt(ncid, vn_rmnc_plasma, 'Amplitudes of the cosine(m*theta-n*zeta) terms in a Fourier expansion of the cylindrical coordinate ' // &
       'R(theta,zeta) defining the plasma surface. The corresponding mode numbers (m,n) are stored in xm_plasma and xn_plasma.')
  if (lasym) then
     call cdf_define(ncid, vn_rmns_plasma, rmns_plasma, dimname=mnmax_plasma_dim)
     call cdf_setatt(ncid, vn_rmns_plasma, 'Amplitudes of the sine(m*theta-n*zeta) terms in a Fourier expansion of the cylindrical coordinate ' // &
          'R(theta,zeta) defining the plasma surface. The corresponding mode numbers (m,n) are stored in xm_plasma and xn_plasma.')
     call cdf_define(ncid, vn_zmnc_plasma, zmnc_plasma, dimname=mnmax_plasma_dim)
     call cdf_setatt(ncid, vn_zmnc_plasma, 'Amplitudes of the cosine(m*theta-n*zeta) terms in a Fourier expansion of the coordinate ' // &
          'Z(theta,zeta) defining the plasma surface. The corresponding mode numbers (m,n) are stored in xm_plasma and xn_plasma.')
  end if
  call cdf_define(ncid, vn_zmns_plasma, zmns_plasma, dimname=mnmax_plasma_dim)
  call cdf_setatt(ncid, vn_zmns_plasma, 'Amplitudes of the sine(m*theta-n*zeta) terms in a Fourier expansion of the coordinate ' // &
       'Z(theta,zeta) defining the plasma surface. The corresponding mode numbers (m,n) are stored in xm_plasma and xn_plasma.')

  call cdf_define(ncid, vn_rmnc_coil, rmnc_coil, dimname=mnmax_coil_dim)
  call cdf_setatt(ncid, vn_rmnc_coil, 'Amplitudes of the cosine(m*theta-n*zeta) terms in a Fourier expansion of the cylindrical coordinate ' // &
       'R(theta,zeta) defining the coil surface. The corresponding mode numbers (m,n) are stored in xm_coil and xn_coil.')
  if (lasym) then
     call cdf_define(ncid, vn_rmns_coil, rmns_coil, dimname=mnmax_coil_dim)
     call cdf_setatt(ncid, vn_rmns_coil, 'Amplitudes of the sine(m*theta-n*zeta) terms in a Fourier expansion of the cylindrical coordinate ' // &
          'R(theta,zeta) defining the inner magnetization surface. The corresponding mode numbers (m,n) are stored in xm_coil and xn_coil.')
     call cdf_define(ncid, vn_zmnc_coil, zmnc_coil, dimname=mnmax_coil_dim)
     call cdf_setatt(ncid, vn_zmnc_coil, 'Amplitudes of the cosine(m*theta-n*zeta) terms in a Fourier expansion of the coordinate ' // &
          'Z(theta,zeta) defining the inner magnetization surface. The corresponding mode numbers (m,n) are stored in xm_coil and xn_coil.')
  end if
  call cdf_define(ncid, vn_zmns_coil, zmns_coil, dimname=mnmax_coil_dim)
  call cdf_setatt(ncid, vn_zmns_coil, 'Amplitudes of the sine(m*theta-n*zeta) terms in a Fourier expansion of the coordinate ' // &
       'Z(theta,zeta) defining the inner magnetization surface. The corresponding mode numbers (m,n) are stored in xm_coil and xn_coil.')

  call cdf_define(ncid, vn_RHS_B, RHS_B, dimname=system_size_dim)

  call cdf_define(ncid, vn_lambda, lambda(1:nsaved), dimname=nsaved_dim)
  call cdf_setatt(ncid, vn_lambda, 'Values of the regularization parameter that were used, in SI units (Tesla^2 / Ampere^2)')

  call cdf_define(ncid, vn_chi2_B, chi2_B(1:nsaved), dimname=nsaved_dim)
  call cdf_setatt(ncid, vn_chi2_B, 'Values of chi^2_B (the area integral over the plasma surface of |B_normal|^2) that resulted for each value of lambda, in SI units (Tesla^2 meter^2)')

  call cdf_define(ncid, vn_chi2_M, chi2_M(1:nsaved), dimname=nsaved_dim)
  call cdf_setatt(ncid, vn_chi2_M, 'Values of chi^2_M (the volume integral over the magnetization region of magnetization squared times thickness d) that resulted for each value of lambda, in SI units (Amperes^2 meters^2).')

  call cdf_define(ncid, vn_max_Bnormal, max_Bnormal(1:nsaved), dimname=nsaved_dim)
  call cdf_setatt(ncid, vn_max_Bnormal, 'Maximum (over the plasma surface) magnetic field normal to the target plasma shape that resulted for each value of lambda, in Tesla.')

  call cdf_define(ncid, vn_max_M, max_M(1:nsaved), dimname=nsaved_dim) ! We only write elements 1:nsaved in case of a lambda search.
  call cdf_setatt(ncid, vn_max_M, 'Maximum (over the magnetization region) magnitude of the magnetization that resulted for each value of lambda, in SI units (Amperes / meter).')

  call cdf_define(ncid, vn_min_M, min_M(1:nsaved), dimname=nsaved_dim) ! We only write elements 1:nsaved in case of a lambda search.
  call cdf_setatt(ncid, vn_min_M, 'Minimum (over the magnetization region) magnitude of the magnetization that resulted for each value of lambda, in SI units (Amperes / meter).')

  ! Arrays with dimension 2

  call cdf_define(ncid, vn_norm_normal_plasma,  norm_normal_plasma,  dimname=ntheta_nzeta_plasma_dim)
  call cdf_setatt(ncid, vn_norm_normal_plasma, '|N|, where N = (d r / d zeta) cross (d r / d theta) is a non-unit-length normal vector ' // &
       'and r is the posiiton vector, for the plasma surface. This quantity is the Jacobian appearing in area integrals: ' // &
       'int d^2a = int dtheta int dzeta |N|. Units = meters^2.')

  call cdf_define(ncid, vn_norm_normal_coil,  norm_normal_coil,  dimname=ntheta_nzeta_coil_dim)
  call cdf_setatt(ncid, vn_norm_normal_coil, '|N|, where N = (d r / d zeta) cross (d r / d theta) is a non-unit-length normal vector ' // &
       'and r is the posiiton vector, for the inner magnetization surface. This quantity is the Jacobian appearing in area integrals: ' // &
       'int d^2a = int dtheta int dzeta |N|. Units = meters^2.')

  call cdf_define(ncid, vn_Bnormal_from_TF_and_plasma_current, Bnormal_from_TF_and_plasma_current, dimname=ntheta_nzeta_plasma_dim)
  call cdf_setatt(ncid, vn_Bnormal_from_TF_and_plasma_current, 'Contribution to the magnetic field normal to the plasma surface from currents inside the plasma ' // &
       'and any electromagnetic coils, such as planar TF coils that generate the net toroidal flux. Units = Tesla.')

  call cdf_define(ncid, vn_mean_curvature_coil,  mean_curvature_coil,  dimname=ntheta_nzeta_coil_dim)
  call cdf_setatt(ncid, vn_mean_curvature_coil,  "Mean curvature (i.e. average of the two principal curvatures) of the inner magnetization surface. Units = 1/meters.")

  if (save_level < 3) then
     call cdf_define(ncid, vn_matrix_B,  matrix_B,  dimname=system_size_system_size_dim)
     call cdf_define(ncid, vn_matrix_regularization,  matrix_regularization,  dimname=system_size_system_size_dim)
  end if

  call cdf_define(ncid, vn_interpolate_magnetization_to_integration,  interpolate_magnetization_to_integration,  dimname=ns_integration_ns_magnetization_dim)
  call cdf_setatt(ncid, vn_interpolate_magnetization_to_integration,  "Matrix to transfer quantities from the s_magnetization radial grid to the s_integration radial grid. Dimensionless.")

  ! Arrays with dimension 3

  call cdf_define(ncid, vn_r_plasma,  r_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)
  call cdf_setatt(ncid, vn_r_plasma, 'Position vector describing the plasma boundary surface, in Cartesian components, with units of meters. ' // &
       'The dimension of size 3 corresponds to the (x,y,z) components.')

  call cdf_define(ncid, vn_r_coil,  r_coil,  dimname=xyz_ntheta_nzetal_coil_dim)
  call cdf_setatt(ncid, vn_r_coil,  'Position vector describing the plasma boundary surface, in Cartesian components, with units of meters. ' // &
       'The dimension of size 3 corresponds to the (x,y,z) components.')

  if (save_level < 3) then
     call cdf_define(ncid, vn_drdtheta_plasma,  drdtheta_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)
     call cdf_define(ncid, vn_drdtheta_coil,  drdtheta_coil,  dimname=xyz_ntheta_nzetal_coil_dim)
     
     call cdf_define(ncid, vn_drdzeta_plasma,  drdzeta_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)
     call cdf_define(ncid, vn_drdzeta_coil,  drdzeta_coil,  dimname=xyz_ntheta_nzetal_coil_dim)

     call cdf_define(ncid, vn_normal_plasma,  normal_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)
     call cdf_define(ncid, vn_normal_coil,  normal_coil,  dimname=xyz_ntheta_nzetal_coil_dim)

  end if

  call cdf_define(ncid, vn_Bnormal_total, Bnormal_total(:,:,1:nsaved), dimname=ntheta_nzeta_plasma_nsaved_dim)
  call cdf_setatt(ncid, vn_Bnormal_total, 'Residual magnetic field normal to the plasma surface, in units of Tesla, ' // &
       'for each value of the regularization parameter lambda considered.')

  call cdf_define(ncid, vn_Jacobian_coil, Jacobian_coil, dimname=ntheta_nzeta_coil_ns_integration_dim)
  call cdf_setatt(ncid, vn_Jacobian_coil, "Jacobian of the (s,theta,zeta) coordinates describing the magnetization region. Units = meters^3.")
 
  ! Arrays with dimension 4

  if (save_level < 1) then
     call cdf_define(ncid, vn_g, g,  dimname=ntheta_nzeta_plasma_basis_ns_magnetization_RZetaZ_dim)
     call cdf_setatt(ncid, vn_g, 'Matrix relating the magnetization to its contribution to B_normal on the plasma surface.')
  end if

  call cdf_define(ncid, vn_abs_M, abs_M, dimname=ntheta_nzeta_coil_ns_magnetization_nsaved_dim)
  call cdf_setatt(ncid, vn_abs_M, 'Magnitude of the magnetization vector, |M|. Units = Amperes / meter.')

  call cdf_define(ncid, vn_magnetization_vector_mn, magnetization_vector_mn, dimname=basis_ns_magnetization_RZetaZ_nsaved_dim)
  call cdf_setatt(ncid, vn_magnetization_vector_mn, 'The magnetization vector, M, with its dependence on (theta,zeta) described in Fourier space. ' // &
       'The dimension of size 3 corresponds to the components of M along the cylindrical unit basis vectors, ' // &
       'with these basis vectors evaluated at toroidal angle zeta, not necessarily at the position where M is evaluated. Units = Amperes / meter.')

  ! Arrays with dimension 5

  call cdf_define(ncid, vn_magnetization_vector, magnetization_vector, dimname=ntheta_nzeta_coil_ns_magnetization_RZetaZ_nsaved_dim)
  call cdf_setatt(ncid, vn_magnetization_vector, 'The magnetization vector, M, with its dependence on (theta,zeta) described in real space rather than Fourier space. ' // &
       'The dimension of size 3 corresponds to the components of M along the cylindrical unit basis vectors, ' // &
       'with these basis vectors evaluated at toroidal angle zeta, not necessarily at the position where M is evaluated. Units = Amperes / meter.')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! Done with cdf_define calls. Now write the data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  ! Scalars

  call cdf_write(ncid, vn_lambda_option, lambda_option)
  call cdf_write(ncid, vn_nfp, nfp)
  call cdf_write(ncid, vn_geometry_option_plasma, geometry_option_plasma)
  call cdf_write(ncid, vn_geometry_option_coil, geometry_option_coil)
  call cdf_write(ncid, vn_ntheta_plasma, ntheta_plasma)
  call cdf_write(ncid, vn_nzeta_plasma, nzeta_plasma)
  call cdf_write(ncid, vn_nzetal_plasma, nzetal_plasma)
  call cdf_write(ncid, vn_ntheta_coil, ntheta_coil)
  call cdf_write(ncid, vn_nzeta_coil, nzeta_coil)
  call cdf_write(ncid, vn_nzetal_coil, nzetal_coil)
  call cdf_write(ncid, vn_ns_magnetization, ns_magnetization)
  call cdf_write(ncid, vn_ns_integration, ns_integration)
  call cdf_write(ncid, vn_system_size, system_size)
  call cdf_write(ncid, vn_a_plasma, a_plasma)
  call cdf_write(ncid, vn_a_coil, a_coil)
  call cdf_write(ncid, vn_R0_plasma, R0_plasma)
  call cdf_write(ncid, vn_R0_coil, R0_coil)
  call cdf_write(ncid, vn_mpol_magnetization, mpol_magnetization)
  call cdf_write(ncid, vn_ntor_magnetization, ntor_magnetization)
  call cdf_write(ncid, vn_mnmax_magnetization, mnmax_magnetization)
  call cdf_write(ncid, vn_mnmax_plasma, mnmax_plasma)
  call cdf_write(ncid, vn_mnmax_coil, mnmax_coil)
  call cdf_write(ncid, vn_num_basis_functions, num_basis_functions)
  call cdf_write(ncid, vn_symmetry_option, symmetry_option)
  call cdf_write(ncid, vn_area_plasma, area_plasma)
  call cdf_write(ncid, vn_area_coil, area_coil)
  call cdf_write(ncid, vn_volume_plasma, volume_plasma)
  call cdf_write(ncid, vn_volume_coil, volume_coil)
  call cdf_write(ncid, vn_nlambda, nlambda)
  call cdf_write(ncid, vn_nd, nd)
  call cdf_write(ncid, vn_nsaved, nsaved)
  call cdf_write(ncid, vn_total_time, total_time)
  call cdf_write(ncid, vn_exit_code, exit_code)
  if (trim(lambda_option)==lambda_option_search) call cdf_write(ncid, vn_chi2_B_target, chi2_B_target)
  call cdf_write(ncid, vn_sign_normal, sign_normal)
  call cdf_write(ncid, vn_d_initial, d_initial)

  ! Arrays with dimension 1

  call cdf_write(ncid, vn_theta_plasma, theta_plasma)
  call cdf_write(ncid, vn_zeta_plasma, zeta_plasma)
  call cdf_write(ncid, vn_zetal_plasma, zetal_plasma)
  call cdf_write(ncid, vn_theta_coil, theta_coil)
  call cdf_write(ncid, vn_zeta_coil, zeta_coil)
  call cdf_write(ncid, vn_zetal_coil, zetal_coil)
  call cdf_write(ncid, vn_s_magnetization, s_magnetization)
  call cdf_write(ncid, vn_s_integration, s_integration)
  call cdf_write(ncid, vn_s_weights, s_weights)
  call cdf_write(ncid, vn_xm_magnetization, xm_magnetization)
  call cdf_write(ncid, vn_xn_magnetization, xn_magnetization)
  call cdf_write(ncid, vn_xm_plasma, xm_plasma)
  call cdf_write(ncid, vn_xn_plasma, xn_plasma)
  call cdf_write(ncid, vn_xm_coil, xm_coil)
  call cdf_write(ncid, vn_xn_coil, xn_coil)
  call cdf_write(ncid, vn_rmnc_plasma, rmnc_plasma)
  if (lasym) then
     call cdf_write(ncid, vn_rmns_plasma, rmns_plasma)
     call cdf_write(ncid, vn_zmnc_plasma, zmnc_plasma)
  end if
  call cdf_write(ncid, vn_zmns_plasma, zmns_plasma)
  call cdf_write(ncid, vn_rmnc_coil, rmnc_coil)
  if (lasym) then
     call cdf_write(ncid, vn_rmns_coil, rmns_coil)
     call cdf_write(ncid, vn_zmnc_coil, zmnc_coil)
  end if
  call cdf_write(ncid, vn_zmns_coil, zmns_coil)
  call cdf_write(ncid, vn_RHS_B, RHS_B)
  call cdf_write(ncid, vn_lambda, lambda(1:nsaved))
  call cdf_write(ncid, vn_chi2_B, chi2_B(1:nsaved))
  call cdf_write(ncid, vn_chi2_M, chi2_M(1:nsaved))
  call cdf_write(ncid, vn_max_Bnormal, max_Bnormal(1:nsaved))
  call cdf_write(ncid, vn_max_M, max_M(1:nsaved))
  call cdf_write(ncid, vn_min_M, min_M(1:nsaved))

  ! Arrays with dimension 2

  call cdf_write(ncid, vn_norm_normal_plasma,  norm_normal_plasma)
  call cdf_write(ncid, vn_norm_normal_coil,  norm_normal_coil)
  call cdf_write(ncid, vn_Bnormal_from_TF_and_plasma_current, Bnormal_from_TF_and_plasma_current)
  call cdf_write(ncid, vn_mean_curvature_coil, mean_curvature_coil)
  if (save_level < 3) then
     call cdf_write(ncid, vn_matrix_B, matrix_B)
     call cdf_write(ncid, vn_matrix_regularization, matrix_regularization)
  end if
  call cdf_write(ncid, vn_interpolate_magnetization_to_integration,  interpolate_magnetization_to_integration)

  ! Arrays with dimension 3

  call cdf_write(ncid, vn_r_plasma, r_plasma)
  call cdf_write(ncid, vn_r_coil, r_coil)

  if (save_level < 3) then
     call cdf_write(ncid, vn_drdtheta_plasma, drdtheta_plasma)
     call cdf_write(ncid, vn_drdtheta_coil, drdtheta_coil)

     call cdf_write(ncid, vn_drdzeta_plasma, drdzeta_plasma)
     call cdf_write(ncid, vn_drdzeta_coil, drdzeta_coil)

     call cdf_write(ncid, vn_normal_plasma, normal_plasma)
     call cdf_write(ncid, vn_normal_coil, normal_coil)

  end if

  call cdf_write(ncid, vn_Bnormal_total, Bnormal_total(:,:,1:nsaved))
  call cdf_write(ncid, vn_Jacobian_coil, Jacobian_coil)

  ! Arrays with dimension 4

  if (save_level<1) then
     call cdf_write(ncid, vn_g,  g)
  end if

  call cdf_write(ncid, vn_abs_M, abs_M)
  call cdf_write(ncid, vn_magnetization_vector_mn, magnetization_vector_mn)

  ! Arrays with dimension 5

  call cdf_write(ncid, vn_magnetization_vector, magnetization_vector)

  ! Finish up:
  call cdf_close(ncid)

end subroutine regcoil_write_output

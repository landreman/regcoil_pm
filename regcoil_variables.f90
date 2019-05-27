module regcoil_variables

  use stel_kinds

  implicit none

  logical :: verbose = .true.

  character(len=*), parameter :: &
       lambda_option_single = "single", &
       lambda_option_scan = "scan", &
       lambda_option_search = "search"
  character(len=200) :: lambda_option = lambda_option_scan

  integer :: ntheta_plasma=64, nzeta_plasma=64, nzetal_plasma
  integer :: ntheta_coil=64, nzeta_coil=64, nzetal_coil

  integer :: geometry_option_plasma = 0
  integer :: geometry_option_coil = 0

  real(dp) :: R0_plasma = 10.0, R0_coil = 10.0
  real(dp) :: a_plasma = 0.5, a_coil = 1.0
  real(dp) :: separation=0.2

  character(len=200) :: wout_filename=""
  character(len=200) :: shape_filename_plasma=""
  character(len=200) :: nescin_filename="nescin.out"
  character(len=200) :: nescout_filename=""
  character(len=200) :: efit_filename=""
  character(len=200) :: output_filename

  real(dp), dimension(:), allocatable :: theta_plasma, zeta_plasma, zetal_plasma
  real(dp), dimension(:,:,:), allocatable :: r_plasma, drdtheta_plasma, drdzeta_plasma, normal_plasma

  real(dp), dimension(:,:,:,:), allocatable :: g, inductance

  real(dp), dimension(:,:), allocatable :: Bnormal_from_TF_and_plasma_current
  real(dp), dimension(:,:), allocatable :: matrix_B, matrix_regularization
  real(dp), dimension(:), allocatable :: RHS_B, RHS_regularization
  real(dp), dimension(:,:,:), allocatable :: Bnormal_total
  real(dp), dimension(:), allocatable :: chi2_B, chi2_M, max_Bnormal, max_M, min_M
  real(dp), dimension(:,:,:,:), allocatable :: M_R, M_zeta, M_Z, abs_M
  real(dp), dimension(:,:,:), allocatable :: M_R_mn, M_zeta_mn, M_Z_mn

  real(dp), dimension(:), allocatable :: theta_coil, zeta_coil, zetal_coil
  real(dp), dimension(:,:,:), allocatable :: r_coil, drdtheta_coil, drdzeta_coil, normal_coil
  real(dp), dimension(:,:,:), allocatable :: d2rdtheta2_coil, d2rdthetadzeta_coil, d2rdzeta2_coil

  real(dp), dimension(:,:), allocatable :: norm_normal_plasma, norm_normal_coil
  real(dp), dimension(:,:), allocatable :: basis_functions

  real(dp) :: dtheta_plasma, dzeta_plasma, dtheta_coil, dzeta_coil

  integer :: mpol_magnetization=12
  integer :: ntor_magnetization=12
  integer :: mnmax_plasma, mnmax_coil, mnmax_magnetization
  integer :: num_basis_functions, system_size
  integer, dimension(:), allocatable :: xm_plasma, xn_plasma, xm_coil, xn_coil, xm_magnetization, xn_magnetization
  real(dp), dimension(:), allocatable :: rmns_plasma, zmnc_plasma, rmnc_plasma, zmns_plasma
  real(dp), dimension(:), allocatable :: rmns_coil, zmnc_coil, rmnc_coil, zmns_coil
  integer :: nfp
  logical :: lasym
  integer :: max_mpol_coil = 24, max_ntor_coil = 24 ! These variables are upper limits on the # of Fourier modes used to describe a uniform-offset coil surface.
  integer :: mpol_coil_filter = 24, ntor_coil_filter = 24

  integer :: save_level = 3
  integer :: nfp_imposed = 1

  integer :: symmetry_option = 3
  real(dp) :: total_time

  integer :: efit_num_modes = 10
  real(dp) :: efit_psiN = 0.98
  real(dp) :: constant_arclength_tolerance = 1.0e-6

  real(dp) :: mpol_transform_refinement=5, ntor_transform_refinement=1
  real(dp) :: area_plasma, area_coil, volume_plasma, volume_coil

  logical :: load_bnorm = .true.
  character(len=200) :: bnorm_filename=""
  real(dp) :: curpol = 1  ! number which multiplies data in bnorm file.
  integer :: nbf ! number of Fourier harmonics in FOCUS format boundary.
  integer, dimension(:), allocatable :: bfn, bfm
  real(dp), dimension(:), allocatable :: bfs, bfc

  integer :: nlambda = 4
  real(dp) :: lambda_min = 1.0d-19, lambda_max = 1.0d-13
  real(dp) :: lambda_single = 0
  real(dp), dimension(:), allocatable :: lambda

  real(dp), dimension(:,:), allocatable :: matrix
  real(dp), dimension(:), allocatable :: RHS, solution

  ! Variables needed by LAPACK:
  integer :: LAPACK_INFO, LAPACK_LWORK
  real(dp), dimension(:), allocatable :: LAPACK_WORK
  integer, dimension(:), allocatable :: LAPACK_IPIV

  real(dp) :: target_value = 8.0d+6
  real(dp) :: lambda_search_tolerance = 1.0d-5
  integer :: exit_code = 0
  real(dp) :: chi2_B_target = 0

  real(dp), dimension(:,:,:), allocatable :: Jacobian_coil
  real(dp), dimension(:,:), allocatable :: mean_curvature_coil
  real(dp), dimension(:,:), allocatable :: d
  integer :: ns_magnetization = 1
  integer :: ns_integration = 2
  real(dp) :: d_initial = 0.01d+0
  real(dp), dimension(:), allocatable :: s_integration, s_weights, s_magnetization
  real(dp), dimension(:,:), allocatable :: interpolate_magnetization_to_integration

  character(len=*), parameter :: &
       target_option_max_M = "max_M", &
       target_option_rms_M = "rms_M", &
       target_option_chi2_M = "chi2_M", &
       target_option_max_Bnormal = "max_Bnormal", &
       target_option_rms_Bnormal = "rms_Bnormal", &
       target_option_chi2_B = "chi2_B"
  character(len=200) :: target_option = target_option_max_M

  character(len=*), parameter :: &
       s_integration_option_Gaussian = "Gaussian", &
       s_integration_option_Chebyshev = "Chebyshev"
  character(len=200) :: s_integration_option = s_integration_option_Gaussian

  namelist / regcoil_nml / ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, &
       geometry_option_plasma, geometry_option_coil, &
       R0_plasma, R0_coil, a_plasma, a_coil, &
       separation, wout_filename, &
       save_level, nfp_imposed, symmetry_option, &
       mpol_magnetization, ntor_magnetization, mpol_coil_filter, ntor_coil_filter, &
       nescin_filename, efit_filename, efit_psiN, efit_num_modes, &
       mpol_transform_refinement, ntor_transform_refinement, max_mpol_coil, max_ntor_coil, &
       load_bnorm, bnorm_filename, &
       shape_filename_plasma, nlambda, lambda_min, lambda_max, lambda_option, verbose, nescout_filename, &
       target_option, target_value, lambda_search_tolerance, &
       ns_magnetization, ns_integration, d_initial, s_integration_option, lambda_single

end module regcoil_variables


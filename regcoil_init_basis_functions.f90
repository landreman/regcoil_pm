subroutine regcoil_init_basis_functions()

  use regcoil_variables
  use regcoil_init_Fourier_modes_mod
  
  implicit none

  integer :: itheta_coil, izeta_coil
  integer :: index_coil, imn, iflag, offset
  integer :: tic, toc, countrate
  real(dp) :: angle, sinangle, cosangle


  call system_clock(tic,countrate)
  if (verbose) print *,"Initializing basis functions"

  ! Initialize Fourier arrays
  call regcoil_init_Fourier_modes(mpol_magnetization, ntor_magnetization, mnmax_magnetization, xm_magnetization, xn_magnetization, .false.)
  xn_magnetization = xn_magnetization * nfp

  allocate(dmnc(mnmax_magnetization))
  allocate(dmns(mnmax_magnetization))
  dmnc = 0
  dmns = 0
  
  select case (symmetry_option)
  case (1,2)
     num_basis_functions = mnmax_magnetization
  case (3)
     num_basis_functions = mnmax_magnetization * 2
  case default
     print *,"Error! Invalid setting for symmetry_option:",symmetry_option
     stop
  end select
  if (include_constant_basis_function) num_basis_functions = num_basis_functions + 1

  system_size = 3 * ns_magnetization * num_basis_functions
  
  if (allocated(basis_functions_R)) deallocate(basis_functions_R)
  allocate(basis_functions_R(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 1!'

  if (allocated(basis_functions_zeta_Z)) deallocate(basis_functions_zeta_Z)
  allocate(basis_functions_zeta_Z(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 1!'

  offset = 0
  if (include_constant_basis_function) then
     ! First basis function is the constant function:
     ! (For stellarator symmetry, this basis function always has amplitude=0 for the zeta and Z components of M, but we include it in the code anyway just so the zeta and Z blocks have the same size as the R block.)
     basis_functions_R(:, 1) = 1
     basis_functions_zeta_Z(:, 1) = 1
     offset = 1
  end if

  ! Remaining basis functions:
  ! These loops could be made faster
  ! by using the sum-angle trig identities and pretabulating the trig functions.
  ! But these loops are not the rate-limiting step, so I'll use the more transparent direct method here.
  select case (symmetry_option)
  case (1,2)
     ! Assume stellarator symmetry.
     do izeta_coil = 1, nzeta_coil
        do itheta_coil = 1, ntheta_coil
           index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
           do imn = 1, mnmax_magnetization
              angle = xm_magnetization(imn)*theta_coil(itheta_coil)-xn_magnetization(imn)*zeta_coil(izeta_coil)
              sinangle = sin(angle)
              cosangle = cos(angle)
              basis_functions_R(index_coil, imn + offset) = sinangle
              basis_functions_zeta_Z(index_coil, imn + offset) = cosangle
           end do
        end do
     end do

  case (3)
     ! Don't assume stellarator symmetry.

     do izeta_coil = 1, nzeta_coil
        do itheta_coil = 1, ntheta_coil
           index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
           do imn = 1, mnmax_magnetization
              angle = xm_magnetization(imn)*theta_coil(itheta_coil)-xn_magnetization(imn)*zeta_coil(izeta_coil)
              sinangle = sin(angle)
              cosangle = cos(angle)
              basis_functions_R(index_coil, imn + offset) = sinangle
              basis_functions_zeta_Z(index_coil, imn + offset) = sinangle
              basis_functions_R(index_coil, imn + offset + mnmax_magnetization) = cosangle
              basis_functions_zeta_Z(index_coil, imn + offset + mnmax_magnetization) = cosangle
           end do
        end do
     end do
  case default
     stop "Should not get here."
  end select
  
  call system_clock(toc)
  if (verbose) print *,"Done. Took",real(toc-tic)/countrate,"sec."


end subroutine regcoil_init_basis_functions

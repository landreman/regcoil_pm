subroutine regcoil_validate_input

  use regcoil_variables
  use safe_open_mod

  implicit none

  integer :: iunit = 7, istat, j
  character(300) :: myline
  character(*), parameter :: matchString = "---- Phi(m,n) for"
  character(len=*), parameter :: line="******************************************************************"
  real(dp) :: typical_target_min, typical_target_max

  if (ntheta_plasma < 1) then
     stop "Error! ntheta_plasma must be >= 1."
  end if

  if (ntheta_coil < 1) then
     stop "Error! ntheta_coil must be >= 1."
  end if


  if (nzeta_plasma < 1) then
     stop "Error! nzeta_plasma must be >= 1."
  end if

  if (nzeta_coil < 1) then
     stop "Error! nzeta_coil must be >= 1."
  end if


  if (mpol_magnetization < 0) then
     stop "Error! mpol_magnetization must be >= 0."
  end if


  if (ntor_magnetization < 0) then
     stop "Error! ntor_magnetization must be >= 0."
  end if


  if (save_level < 0) then
     stop "Error! save_level must be >= 0."
  end if

  if (save_level > 3) then
     stop "Error! save_level must be <= 3."
  end if


  if (symmetry_option < 1) then
     stop "Error! symmetry_option must be >= 1."
  end if

  if (symmetry_option > 3) then
     stop "Error! symmetry_option must be <= 3."
  end if


  if (geometry_option_plasma < 0) then
     stop "Error! geometry_option_plasma must be >= 0."
  end if

  if (geometry_option_plasma > 7) then                  ! czhu changed this to option 7 for FOCUS
     stop "Error! geometry_option_plasma must be <= 7." 
  end if


  if (geometry_option_coil < 0) then
     stop "Error! geometry_option_coil must be >= 0."
  end if

  if (geometry_option_coil > 4) then
     stop "Error! geometry_option_coil must be <= 4."
  end if



  if (separation < 0) then
     stop "Error! separation must be >= 0."
  end if



  if (load_bnorm) then
     ! To use bnorm data, we must have a VMEC file to set the overall normalization
     select case (geometry_option_plasma)
     case (2,3,4,7)
        ! Yes, we have a VMEC file available.
     case (0,1,5)
        stop "Error! If load_bnorm=.t., the plasma surface must come from a vmec wout file."
     case default
        stop "Error! Invalid geometry_option_plasma"
     end select
  end if

  if (nlambda < 1) then
     stop "nlambda must be at least 1."
  end if

  if (lambda_min <= 0) then
     stop "lambda_min must be greater than 0."
  end if

  if (lambda_max < lambda_min) then
     stop "lambda_max must be >= lambda_min."
  end if

!!$  if (general_option<1) then
!!$     stop "general_option must be at least 1."
!!$  end if
!!$  if (general_option>5) then
!!$     stop "general_option must be no more than 5."
!!$  end if

  ! General_option=2 does not make sense for permanent magnets
!!$  if (general_option==2) then
!!$     ! Replace nlambda with the number of current potentials saved in the nescout file.
!!$     if (verbose) print *,"Opening nescout file",nescout_filename
!!$     call safe_open(iunit, istat, trim(nescout_filename), 'old', 'formatted')
!!$     if (istat .ne. 0) then
!!$        stop 'Error opening nescout file'
!!$     endif
!!$     j = 0
!!$     do
!!$        read (iunit,"(a)",iostat=istat) myline
!!$        if (istat<0) exit
!!$        if (myline(:len(matchString)) == matchString) then
!!$           j = j + 1
!!$        end if
!!$     end do
!!$     if (verbose) print *,"Detected",j,"current potentials in the nescout file."
!!$     nlambda = j
!!$  end if

  if (target_value<=0) then
     stop "target_value must be positive."
  end if

!!$  if (general_option == 4 ) then
!!$     print *,"It is recommended that you run with general_option=5 instead of 4"
!!$     print *,"to verify that this value of target_value is attainable."
!!$  end if

  !if (general_option==4 .or. general_option==5) then
  if (trim(lambda_option)==lambda_option_search) then
     select case (trim(target_option))
     case (target_option_max_M)
        typical_target_min = 1e5
        typical_target_max = 3e8
     case (target_option_chi2_M)
        typical_target_min = 1e14
        typical_target_max = 1e17
     case (target_option_max_Bnormal,target_option_rms_Bnormal)
        typical_target_min = 1e-5
        typical_target_max = 5
     case (target_option_chi2_B)
        typical_target_min = 1e-6
        typical_target_max = 100
     case default
        print *,"Invalid target_option: ",target_option
        stop
     end select

     if (target_value < typical_target_min) then
        print "(a)",line
        print "(a)","Warning! The value of target_value you have set is surprisingly small."
        print "(a)",line
     end if
     if (target_value > typical_target_max) then
        print "(a)",line
        print "(a)","Warning! The value of target_value you have set is surprisingly large."
        print "(a)",line
     end if
  end if

  if (ns_integration < ns_magnetization) then
     stop "ns_integration must be >= ns_magnetization"
  end if

  select case (trim(lambda_option))
  case (lambda_option_single)
  case (lambda_option_scan)
  case (lambda_option_search)
  case default
     print *,"Error! Unrecognized lambda_option: ",trim(lambda_option)
     stop
  end select

  select case (trim(s_integration_option))
  case (s_integration_option_Gaussian)
  case (s_integration_option_Chebyshev)
  case default
     print *,"Error! Unrecognized s_integration_option: ",trim(s_integration_option)
     stop
  end select

  select case (trim(d_option))
  case (d_option_uniform)
     nd = 1
  case (d_option_Picard)
  case (d_option_Anderson_old)
  case (d_option_Anderson)
  case default
     print *,"Error! Unrecognized d_option: ",trim(d_option)
     stop
  end select

  !if (symmetry_option .ne. 3) stop "Error! Presently only symmetry_option=3 works."

  if ((sign_normal.ne.1) .and. (sign_normal.ne.-1)) stop "Error! sign_normal must be 1 or -1."

  if (nd < 1) stop "nd must be at least 1."

  if (target_mu0_M < 0) stop "target_mu0_M must be positive."

  if (trim(magnet_type) == 'qhex') then
     if (ns_magnetization /= 1) then
        print *, "Warning: ns_magnetization > 1 is not supported for qhex geometry. Overriding input value."
        ns_magnetization = 1
     end if
     if (ns_integration /= 1) then
        print *, "Warning: ns_integration > 1 is not supported for qhex geometry. Overriding input/default value."
        ns_integration = 1
     end if
     if (filter_d) then
        print *, "Warning: Fourier filtering of d is not supported for qhex geometry. Setting filter_d to .false."
        filter_d = .false.
     end if
  end if

end subroutine regcoil_validate_input

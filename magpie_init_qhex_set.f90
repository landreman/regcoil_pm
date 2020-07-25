subroutine magpie_init_qhex_set()

  use regcoil_variables
  
  use stel_kinds
  use stel_constants
  use omp_lib

  use magpie_globals,  only: len_fname, len_mag_label, stell_symm, tor_symm, &
                             limiting_surf_file, limiting_surf_phi_sign
  use input,           only: read_namelists, read_surface_file
  use surface_calc,    only: surface
  use build_qhex,      only: build_qhex_set_wrap
  use intersections,   only: intersections_to_check
  use repetition,      only: geometry_repetition_qhex
  
  implicit none
  
  type(surface) :: limiting_surf
  integer :: iflag, i, tic, toc, countrate

  call system_clock(tic,countrate)
  if (verbose) print *,"Initializing set of hexhedral magnets."

  ! Determine which intersections (ports, etc.) should be checked, if any
  call intersections_to_check()

  ! Import the surface around which the hexahedra will be built
  call read_surface_file(limiting_surf_file, nfp, limiting_surf_phi_sign, &
                         stell_symm, limiting_surf)

  ! Construct the qhex array
  call build_qhex_set_wrap(limiting_surf, nzeta_coil, qhex_arr_base)
  nzeta_coil_inv = 1.0 / real(nzeta_coil, dp)

  ! Maximum allowable dipole moment
  allocate(qhex_max_moment(nzeta_coil))
  do i = 1, nzeta_coil
     qhex_max_moment(i) = qhex_arr_base(i)%vol * target_mu0_M / mu0
  end do

  ! Extend the array to wrap around the torus
  call geometry_repetition_qhex('torus', nfp, stell_symm, tor_symm, &
           nzeta_coil, qhex_arr_base, nzetal_coil, qhex_arr)

  ntheta_coil   = 1 

  if (allocated(theta_coil)) deallocate(theta_coil)
  allocate(theta_coil(ntheta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 1'

  if (allocated(zeta_coil)) deallocate(zeta_coil)
  allocate(zeta_coil(nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 2'

  if (allocated(zetal_coil)) deallocate(zetal_coil)
  allocate(zetal_coil(nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 3'

  do i = 1,ntheta_coil
     theta_coil(i) = twopi*(i-1.0_dp)/ntheta_coil
  end do

  do i = 1,nzeta_coil
     zeta_coil(i) = atan2(qhex_arr_base(i)%oy, qhex_arr_base(i)%ox)
  end do

  do i = 1,nzetal_coil
     zetal_coil(i) = atan2(qhex_arr(i)%oy, qhex_arr(i)%ox)
  end do

  ! Should be unity (if ever used); typically they will be multiplied by 
  ! factors that already include the volume in m^3 of the hexahedron
  dtheta_coil = 1._dp
  dzeta_coil  = 1._dp

  ! d = the "height" property of each hexahedron
  allocate(d(ntheta_coil, nzeta_coil))
  allocate(max_d_before_singularity(ntheta_coil, nzeta_coil))
  do i = 1, nzeta_coil
     d(1,i) = qhex_arr_base(i)%height
  end do
  print '(a,es10.3,a,es10.3)', "Max d = ", maxval(d), "; Min d = ", minval(d)

  ! Forbid d from exceeding the initialized height of each hexahedron
  max_d_before_singularity = d

  ! Other stuff normally initialized in regcoil_evaluate_coil_surface
  allocate(interpolate_magnetization_to_integration(ns_magnetization, ns_integration))
  interpolate_magnetization_to_integration = 1._dp

  allocate(s_magnetization_weights(ns_magnetization), s_magnetization(ns_magnetization))
  s_magnetization_weights = 1._dp
  s_magnetization = 1._dp

  allocate(s_integration(ns_integration), s_weights(ns_integration))
  s_integration = 1._dp
  s_weights = 1._dp


end subroutine magpie_init_qhex_set

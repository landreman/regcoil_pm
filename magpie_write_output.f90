subroutine magpie_write_output

  use regcoil_variables
  use magpie_globals,    only: output_base_name, repeat_to_fill, stell_symm, tor_symm, &
                               M_max, focus_stell_symm, focus_tor_symm, focus_no_symm, &
                               len_mag_label, rho_initial, rho_free, &
                               axis_free, corners_file, focus_file
  use qhex_properties,   only: quad_hexahedron
  use output,            only: write_corners_file, write_focus_file
  use magnet_properties, only: magnet, qhex_magnet
  use repetition,        only: geometry_repetition_qhex

  implicit none

  integer :: nQhexTot, i
  type(quad_hexahedron), dimension(:), allocatable :: qhex_arr_out
  type(magnet), dimension(:), allocatable :: magnets
  logical, dimension(:), allocatable :: temp
  character(len=len_mag_label) :: label
  real(dp) :: Mx, My, Mz, m_theta, m_phi
  integer :: focus_symm_key

  ! Qhex corners file (geometric data for the hexahedra)
  call geometry_repetition_qhex(repeat_to_fill, nfp, stell_symm, tor_symm, &
                                nzeta_coil, qhex_arr_base, nQhexTot, qhex_arr_out)
  allocate(temp(nQhexTot))
  if (corners_file) &
     call write_corners_file(qhex_arr_out, nQhexTot, .false., .false., temp, output_base_name)
  deallocate(temp)

  ! Focus file for the magnets
  allocate(magnets(nzeta_coil))

  if (stell_symm) then
     focus_symm_key = focus_stell_symm
  else if (tor_symm) then
     focus_symm_key = focus_tor_symm
  else
     focus_symm_key = focus_no_symm
  end if

  do i = 1, nzeta_coil
     write(label, fmt='(A, I10.10)') 'pm_', i
     magnets(i) = qhex_magnet(qhex_arr_base(i), M_max, focus_symm_key, rho_free, axis_free, &
                              rho_initial, label)

     magnets(i)%rho_initial = abs_M(1,i,1,nsaved) / M_max

     ! Cartesian components of the magnetization vector
     Mx = magnetization_vector(1,i,1,1,nsaved) * cos(zeta_coil(i)) &
           - magnetization_vector(1,i,1,2,nsaved) * sin(zeta_coil(i))
     My = magnetization_vector(1,i,1,1,nsaved) * sin(zeta_coil(i)) &
           + magnetization_vector(1,i,1,2,nsaved) * cos(zeta_coil(i))
     Mz = magnetization_vector(1,i,1,3,nsaved)

     ! Polar and azimuthal angle of the magnetization vector
     m_phi   = atan2(My, Mx)
     m_theta = atan2(sqrt(Mx*Mx + My*My), Mz)

     magnets(i)%phi = m_phi
     magnets(i)%theta = m_theta
     magnets(i)%nx = Mx / abs_M(1,i,1,nsaved)
     magnets(i)%ny = My / abs_M(1,i,1,nsaved)
     magnets(i)%nz = Mz / abs_M(1,i,1,nsaved)
  end do

  if (focus_file) &
     call write_focus_file(magnets, nzeta_coil, output_base_name) 

  deallocate(magnets)

end subroutine magpie_write_output


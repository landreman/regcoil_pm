!-------------------------------------------------------------------------------
! Construction of a qhex magnet set
!-------------------------------------------------------------------------------
module build_qhex

implicit none

contains

subroutine build_qhex_set(reference_surf, limiting_surf, &
               grid, nQhexes, qhexes, nOverlaps, overlaps_out, &
               nMagnets, magnets)

    use magpie_globals,                                                       &
        only: len_mag_label, dp, nfp,                                         &
              qhex_nTheta, qhex_nPhi, qhex_unif_poloidal,                     &
              qhex_incl_upper, qhex_adj_interval, gap_rad, qhex_nSplitLo,     &
              qhex_max_ht, qhex_max_ht_up, qhex_max_ht_lo, qhex_nSplitUp,     &
              M_max, stell_symm, tor_symm, focus_no_symm, focus_tor_symm,     &
              focus_stell_symm, rho_initial, rho_free, axis_free,             &
              output_base_name, focus_file, corners_file, base_grid_file,     &
              incl_err_corners, incl_ovl_corners, max_ovl_check_interval     
    use surface_calc,      only: surface
    use base_grid,         only: base_grid_struct, build_base_grid
    use qhex_properties,   only: quad_hexahedron, qhex_array
    use qhex_array,        only: build_qhex_array, split_qhex
    use magnet_properties, only: magnet, qhex_magnet
    use qhex_overlap,      only: object_overlap

    implicit none

    type(surface), intent(IN) :: limiting_surf, reference_surf
    integer, intent(OUT) :: nQhexes, nOverlaps, nMagnets
    type(quad_hexahedron), dimension(:), allocatable, intent(OUT) :: qhexes
    type(magnet), dimension(:), allocatable, intent(OUT) :: magnets
    logical, dimension(:), allocatable, intent(OUT) :: overlaps_out
    type(base_grid_struct), intent(OUT) :: grid
    integer :: nSplitQhexesUp, nSplitQhexesLo, nSplitQhexes
    integer :: i, j, i_pol, i_tor, idx, n, nErrQhexes, focus_symm_key
    type(qhex_array)       :: qhArr
    type(quad_hexahedron), allocatable, dimension(:) :: all_qhexes
    type(quad_hexahedron), allocatable, dimension(:) :: qh_subs_lo, qh_subs_up
    character(len=len_mag_label) :: label
    real(dp) :: max_ht_lo, max_ht_up
    logical, allocatable, dimension(:) :: overlaps 

    ! Construct the base grid 
    call build_base_grid(reference_surf, qhex_nTheta+1, qhex_nPhi+1, &
                         qhex_unif_poloidal, grid)

    ! Construct the qhex array according to surface input
    max_ht_lo = qhex_max_ht
    max_ht_up = qhex_max_ht
    if (qhex_max_ht_lo> 0.0) max_ht_lo = qhex_max_ht_lo
    if (qhex_max_ht_up> 0.0) max_ht_up = qhex_max_ht_up
    call build_qhex_array(reference_surf, limiting_surf, grid, &
             qhex_incl_upper, max_ht_lo, max_ht_up, qhex_adj_interval, qhArr)

    ! Split the qhexes in the grid and put them in an 1D array
    nSplitQhexesLo = qhArr%nLower * qhex_nSplitLo
    nSplitQhexesUp = qhArr%nUpper * qhex_nSplitUp
    nSplitQhexes = nSplitQhexesUp + nSplitQhexesLo
    write(*,fmt='(A)') '    Subdividing hexahedra to final magnet geometry'
    if (qhArr%incl_upper) then
        write(*,fmt='(A, I0, A, I0, A, I0, A)') '        (', nSplitQhexes, &
            ' total; ', nSplitQhexesLo, ' lower and ', nSplitQhexesUp, ' upper)'
    else
        write(*,fmt='(A, I0, A)') '        (', nSplitQhexes, ' total)'
    end if
    nErrQhexes = 0
    allocate(all_qhexes(nSplitQhexes))
    allocate(qh_subs_lo(qhex_nSplitLo), qh_subs_up(qhex_nSplitUp))
    do i = 1, qhArr%nLower
        i_pol = mod(i-1, qhArr%nTheta) + 1
        i_tor = (i-1)/qhArr%nTheta + 1
        call split_qhex(qhArr%pairs(i_pol, i_tor)%lo, qhex_nSplitLo, gap_rad, &
                        qh_subs_lo)
        do j = 1, qhex_nSplitLo
            idx = (i-1)*qhex_nSplitLo + j
            all_qhexes(idx) = qh_subs_lo(j)
            if (all_qhexes(idx)%err) nErrQhexes = nErrQhexes + 1
        end do
        if (qhArr%incl_upper) then
            call split_qhex(qhArr%pairs(i_pol, i_tor)%up, qhex_nSplitUp, &
                            gap_rad, qh_subs_up)
            do j = 1, qhex_nSplitUp
                idx = (i-1)*qhex_nSplitUp + j + nSplitQhexesLo
                all_qhexes(idx) = qh_subs_up(j)
                if (all_qhexes(idx)%err) nErrQhexes = nErrQhexes + 1
            end do
        end if
    end do

    ! Check the qhexes for intersections with applicable ports & other obstacles
    write(*,fmt='(A)') '    Checking for overlaps with ports and other objects'
    allocate(overlaps(nSplitQhexes))
    nOverlaps = 0
    do i = 1, nSplitQhexes
        call object_overlap(all_qhexes(i), max_ovl_check_interval, overlaps(i))
        if ((.not. all_qhexes(i)%err) .and. overlaps(i)) then
            nOverlaps = nOverlaps + 1
        end if
    end do
    write(*,fmt='(A, I0, A)') '        ', nOverlaps, ' magnets eliminated ' &
                              // 'due to overlaps'

    ! Make an array of magnets from the non-erroneous, non-overlapping qhexes
    nMagnets = nSplitQhexes - nErrQhexes - nOverlaps
    if (stell_symm) then
        focus_symm_key = focus_stell_symm
    else if (tor_symm) then
        focus_symm_key = focus_tor_symm
    else
        focus_symm_key = focus_no_symm
    end if
    allocate(magnets(nMagnets))
    n = 0
    do i = 1, nSplitQhexes
        if ( (.not. all_qhexes(i)%err) .and. (.not. overlaps(i)) ) then
            n = n + 1
            write(label, fmt='(A, I10.10)') 'pm_', i
            magnets(n) = qhex_magnet(all_qhexes(i), M_max, focus_symm_key, &
                                 rho_free, axis_free, rho_initial, label)
        end if
    end do

    ! Prepare the output qhex and overlaps array
    nQhexes = nMagnets
    if (incl_ovl_corners) nQhexes = nQhexes + nOverlaps
    if (incl_err_corners) nQhexes = nQhexes + nErrQhexes
    allocate(qhexes(nQhexes), overlaps_out(nQhexes))
    n = 0
    do i = 1, nSplitQhexes
        if ((incl_ovl_corners .or. (.not. overlaps(i))) .and. &
            (incl_err_corners .or. (.not. all_qhexes(i)%err))) then
            n = n + 1
            overlaps_out(n) = overlaps(i)
            qhexes(n) = all_qhexes(i)
        end if
    end do

    ! Deallocate arrays
    deallocate(all_qhexes, overlaps)

end subroutine build_qhex_set

subroutine build_qhex_set_wrap(limiting_surf, nQhexes, qhexes)

    use surface_calc,      only: surface
    use base_grid,         only: base_grid_struct
    use qhex_properties,   only: quad_hexahedron
    use magnet_properties, only: magnet

    implicit none

    type(surface), intent(IN) :: limiting_surf
    integer, intent(OUT) :: nQhexes
    type(quad_hexahedron), dimension(:), allocatable, intent(OUT) :: qhexes
    integer :: nOverlaps, nMagnets
    type(base_grid_struct) :: grid
    type(magnet), dimension(:), allocatable :: magnets
    logical, dimension(:), allocatable :: overlaps

    call build_qhex_set(limiting_surf, limiting_surf, grid, nQhexes, &
                             qhexes, nOverlaps, overlaps, nMagnets, magnets)

    deallocate(overlaps, magnets)

end subroutine build_qhex_set_wrap

end module build_qhex


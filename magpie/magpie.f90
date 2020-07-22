!-------------------------------------------------------------------------------
! magpie.f90
!
! Constructs arrangements of quadrilateral-faced hexahedra that represent the
! geometry of potential arrays of permanent magnets for stellarators.
!
! Documentation: kchammond.github.io/MAGPIE/
!-------------------------------------------------------------------------------
program magpie

    use magpie_globals,                                                       &
        only: version, dp, len_fname,                                         &
              limiting_surf_file, limiting_surf_phi_sign,                     &
              reference_surf_file, reference_surf_phi_sign,                   &
              nfp, magnet_type, stell_symm, output_base_name,                 &
              focus_file, corners_file, cbrick_file, base_grid_file,          &
              magnet_vtk, wframe_vtk, incl_ovl_cbricks,                       &
              incl_err_corners, incl_ovl_corners, repeat_to_fill,             &
              stell_symm, tor_symm, incl_ovl_trecs
    use input,             only: read_namelists, read_surface_file           
    use surface_calc,      only: surface
    use intersections,     only: intersections_to_check
    use base_grid,         only: base_grid_struct
    use qhex_properties,   only: quad_hexahedron
    use cbrick_properties, only: cbrick
    use trec_properties,   only: trec
    use magnet_properties, only: magnet
    use repetition,        only: geometry_repetition_qhex_ovl
    use build_qhex,        only: build_qhex_set
    use output,            only: write_base_grid_file, write_corners_file,    &
                                 write_focus_file, write_cbrick_file,         &
                                 write_cbrick_vtk, write_cbrick_wf_vtk,       &
                                 write_trec_corners_file

    implicit none

    character(len=len_fname) :: infile
    integer :: cmd_arg_stat, nQhexes, nQhexes_rep, nMagnets, nOverlaps
    integer :: nCbricks, nTrecs
    type(surface) :: limiting_surf, reference_surf
    logical       :: has_limiting_surf, has_reference_surf
    type(base_grid_struct) :: bg
    type(quad_hexahedron), dimension(:), allocatable :: qhexes, qhexes_rep
    type(cbrick), dimension(:), allocatable :: cbricks
    type(trec), dimension(:), allocatable :: trecs
    type(magnet), dimension(:), allocatable :: magnets
    logical, dimension(:), allocatable :: overlaps, overlaps_rep

    ! Obtain the command line argument (input file name)
    call get_command_argument(1, infile, status=cmd_arg_stat)
    if (cmd_arg_stat /= 0) stop 'Usage: xmagpie [infile]'

    ! Print title and version
    write(*,fmt='(A,A,A)') '--------------- MAGPIE version ', trim(version), &
                           ' ---------------'

    ! Read namelists from input file
    call read_namelists(infile)

    ! Determine which intersections should be checked, if any
    call intersections_to_check()

    ! Read data on limiting and/or reference surfaces as applicable
    has_limiting_surf = .false.
    has_reference_surf = .false.
    if (trim(limiting_surf_file) /= '') then
        call read_surface_file(limiting_surf_file, nfp, &
                 limiting_surf_phi_sign, stell_symm, limiting_surf)
        has_limiting_surf = .true.
    end if
    if (trim(reference_surf_file) /= '') then
        call read_surface_file(reference_surf_file, nfp, &
                 reference_surf_phi_sign, stell_symm, reference_surf)
        has_reference_surf = .true.
    end if

    ! If only one surface is provided, use for both ref and lim surf
    if (has_reference_surf .and. .not. has_limiting_surf) then
        limiting_surf = reference_surf
    else if (has_limiting_surf .and. .not. has_reference_surf) then
        reference_surf = limiting_surf
    else if (.not. has_reference_surf .and. .not. has_limiting_surf) then
        write(*,*) 'Error: no reference or limiting surface data provided'
        stop
    end if

    ! Proceed according to the magnet type
    select case (trim(magnet_type))
        case ('qhex')

            call build_qhex_set(reference_surf, limiting_surf, bg, &
                     nQhexes, qhexes, nOverlaps, overlaps, nMagnets, magnets)

            ! Handle repetition of the magnets for geometric output files
            call geometry_repetition_qhex_ovl( &
                             repeat_to_fill, nfp, stell_symm, tor_symm,  &
                             nQhexes,     qhexes,     overlaps, &
                             nQhexes_rep, qhexes_rep, overlaps_rep)

        case ('cbrick')
            call build_cbrick_set(reference_surf, nCbricks, cbricks, &
                                  nMagnets, magnets, nOverlaps, overlaps)
        case ('trec')
            call build_trec_set(reference_surf, nTrecs, trecs, &
                                nMagnets, magnets, nOverlaps, overlaps)
        case default
            write(*,*) 'Error: unrecognized magnet_type'
            stop
    end select

    ! Write output files
    select case (trim(magnet_type))
        case ('qhex')
            if (corners_file) &
                call write_corners_file(qhexes_rep, nQhexes_rep,       &
                         incl_err_corners, incl_ovl_corners, overlaps, &
                         output_base_name)
            if (base_grid_file) &
                call write_base_grid_file(bg, output_base_name)
        case ('cbrick')
            if (cbrick_file) &
                call write_cbrick_file(cbricks, nCbricks, incl_ovl_cbricks, &
                                       overlaps, output_base_name)
            if (magnet_vtk)  &
                call write_cbrick_vtk(cbricks, nCbricks, incl_ovl_cbricks, &
                                      overlaps, output_base_name)
            if (wframe_vtk) &
                 call write_cbrick_wf_vtk(cbricks, nCbricks, incl_ovl_cbricks, &
                                          overlaps, output_base_name)

        case ('trec')
            if (corners_file) &
                call write_trec_corners_file(trecs, nTrecs, incl_ovl_trecs, &
                                             overlaps, output_base_name)
    end select
    if (focus_file) call write_focus_file(magnets, nMagnets, output_base_name)

    ! Deallocate
    deallocate(overlaps, magnets)
    if (trim(magnet_type) == 'qhex'  ) deallocate(qhexes, qhexes_rep)
    if (trim(magnet_type) == 'cbrick') deallocate(cbricks)
    if (trim(magnet_type) == 'trec'  ) deallocate(trecs)

    ! Print summary
    write(*,fmt='(A)') '--------------- MAGPIE finished ---------------'
    
end program magpie


  

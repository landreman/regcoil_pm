!-------------------------------------------------------------------------------
! Construction of a trec magnet set
!-------------------------------------------------------------------------------
subroutine build_trec_set(lim_surf, nTrecs_out, trecs_out, &
                          nMagnets, magnets, nColliding, colliding_out)

    use magpie_globals,    &
        only: dp, len_mag_label, &
              trec_rmin, trec_rmax, trec_zmin, trec_zmax, trec_pmin,          &
              trec_pmax, trec_dr, trec_dz, trec_nPhi, gap_trr, gap_trz,       &
              gap_trp, gap_lim, trec_lpMin, trec_nTheta, radial_extent, pi,   &
              M_max, stell_symm, tor_symm, focus_no_symm, focus_tor_symm,     &
              focus_stell_symm, rho_initial, rho_free, axis_free,             &
              output_base_name, focus_file, corners_file,                     &
              incl_ovl_trecs, max_ovl_check_interval,                         &
              trec_ax_init
    use surface_calc,      only: surface, surf_r, surf_z
    use surface_solve,     only: surf_perp_intersect_2d
    use trec_properties,   only: trec, trec_grid_params, new_trec,            &
                                 trec_in_polygon_5plane,                      &
                                 trec_from_polygon_5plane,                    &
                                 estimate_trec_theta, trec_set_perpendicular, &
                                 trec_obj_overlap, iIB, iIF, iCM, iOB, iOF
    use magnet_properties, only: magnet, trec_magnet
    use input,             only: to_lowercase
    use output,            only: write_focus_file, write_trec_corners_file

    implicit none

    type(surface), intent(IN) :: lim_surf
    integer, intent(OUT) :: nTrecs_out, nMagnets, nColliding
    type(trec), dimension(:), allocatable :: trecs_out
    logical, dimension(:), allocatable, intent(OUT) :: colliding_out
    type(trec), dimension(:,:,:), allocatable :: trArr
    logical, dimension(:,:,:), allocatable :: in_bounds
    type(trec), dimension(:), allocatable :: trecs 
    logical, dimension(:), allocatable :: colliding 
    type(magnet), dimension(:), allocatable :: magnets
    real(dp) :: rTri, r0_grid, z0_grid, trec_dphi, r_cm, z_cm, lr, lp, lz
    real(dp), dimension(:), allocatable :: phi_cm
    integer :: nr, nz, i, j, k, l, n, nRects, focus_symm_key
    real(dp) :: dTheta, phi_ib, phi_ob, phi_if, phi_of, theta_cm, minDist
    real(dp), dimension(trec_nTheta) :: theta
    real(dp), dimension(trec_nTheta,5) :: rSect, zSect
    character(len=len_mag_label) :: label

    write(*,fmt='(A)') '    Constructing grid of trecs within allotted volume'

    ! Determine radial and vertical dimensions of the grid of potential bricks
    trec_dphi = (trec_pmax - trec_pmin)/real(trec_nPhi,dp)
    call trec_grid_params(trec_dr, trec_dz, trec_dphi, 0.5*gap_trp,         &
                          trec_lpMin, trec_rmin, trec_rmax, trec_zmin,      &
                          trec_zmax, rTri, r0_grid, z0_grid, nr, nz)

    allocate(trArr(nr,nz,trec_nphi), in_bounds(nr,nz,trec_nphi), &
             phi_cm(trec_nphi))

    ! Toroidal angles at which the bricks may be centered
    do k = 1, trec_nPhi
        phi_cm(k) = trec_pmin + (k-0.5)*trec_dphi
    end do


    ! Populate the array with potential bricks
    do i = 1, nr
        ! Note that rects are "pushed" against the inner side of the cell
        ! in the radial dimension rather than being centered
        r_cm = r0_grid + (i-0.5)*trec_dr - 0.5*gap_trr
        lr   = trec_dr - gap_trr
        lp   = 2 * ( (r_cm - 0.5*lr - rTri) * sin(0.5*trec_dphi) )
        do j = 1, nz
            ! By contrast, rects are centered within cells in the z dimension
            z_cm = z0_grid + (j-0.5)*trec_dz + 0.5*gap_trz
            lz   = trec_dz - gap_trz
            do k = 1, trec_nPhi
                trArr(i,j,k) = new_trec(r_cm, phi_cm(k), z_cm, lr, lp, lz)
            end do
        end do
    end do

    write(*,fmt='(A)') '    Identifying prisms within acceptable bounds'

    ! Check whether each brick in the grid is within the acceptable bounds
    nRects = 0
    dTheta = 2*pi/(trec_nTheta-1)
    do i = 1, trec_nTheta
        theta(i) = (i-1)*dTheta
    end do
    do k = 1, trec_nPhi

        ! Xsects of limiting surface at phi angle of centroids for tor. sector k
        do l = 1, trec_nTheta
            rSect(l,iCM) = surf_r(lim_surf, theta(l), phi_cm(k))
            zSect(l,iCM) = surf_z(lim_surf, theta(l), phi_cm(k))
        end do

        do i = 1, nr

            ! Limiting surf xsects at toroidal angles of corners
            phi_ib = atan2(trArr(i,1,k)%yib, trArr(i,1,k)%xib)
            phi_ob = atan2(trArr(i,1,k)%yob, trArr(i,1,k)%xob)
            phi_if = atan2(trArr(i,1,k)%yif, trArr(i,1,k)%xif)
            phi_of = atan2(trArr(i,1,k)%yof, trArr(i,1,k)%xof)
            do l = 1, trec_nTheta
                rSect(l,iIB) = surf_r(lim_surf, theta(l), phi_ib)
                zSect(l,iIB) = surf_z(lim_surf, theta(l), phi_ib)
                rSect(l,iIF) = surf_r(lim_surf, theta(l), phi_if)
                zSect(l,iIF) = surf_z(lim_surf, theta(l), phi_if)
                rSect(l,iOB) = surf_r(lim_surf, theta(l), phi_ob)
                zSect(l,iOB) = surf_z(lim_surf, theta(l), phi_ob)
                rSect(l,iOF) = surf_r(lim_surf, theta(l), phi_of)
                zSect(l,iOF) = surf_z(lim_surf, theta(l), phi_of)
            end do

            do j = 1, nz

                ! Eliminate prisms with corners inside the limiting surface
                if (trec_in_polygon_5plane(trArr(i,j,k), trec_nTheta, &
                                           rSect, zSect)) then
                    in_bounds(i,j,k) = .false.

                ! Eliminate prisms with corners too close/far outside lim. surf.
                else
                    minDist = trec_from_polygon_5plane(trArr(i,j,k), &
                                  trec_nTheta, rSect, zSect)
                    if (minDist > gap_lim .and. minDist < radial_extent) then

                        in_bounds(i,j,k) = .true.

                        ! For prisms to be retained, set axis perp. to surface
                        ! if applicable
                        if (to_lowercase(trec_ax_init) == 'normal') then
                            theta_cm = &
                                estimate_trec_theta(trArr(i,j,k), trec_nTheta, &
                                    theta, rSect(:,iCM), zSect(:,iCM))
                            call trec_set_perpendicular(trArr(i,j,k), &
                                     lim_surf, minDist, theta_cm, phi_cm(k))
                        end if

                        nRects = nRects + 1
                    else
                        in_bounds(i,j,k) = .false.
                    end if
                end if

            end do
        end do

    end do

    write(*,fmt='(A,I0,A)') '        ', nRects, ' prisms retained'

    ! Determine which of the in-bounds prisms overlaps ports or other objects
    write(*,fmt='(A)') &
        '    Checking for collisions with ports and other objects'
    allocate(trecs(nRects), colliding(nRects))
    n = 0
    nMagnets = 0
    do i = 1, nr
        do j = 1, nz
            do k = 1, trec_nphi
                if (in_bounds(i,j,k)) then
                    n = n + 1
                    trecs(n) = trArr(i,j,k)
                    call trec_obj_overlap(trecs(n), &
                                           max_ovl_check_interval, colliding(n))
                    if (.not. colliding(n)) nMagnets = nMagnets + 1
                end if
            end do
        end do
    end do
    nColliding = nRects - nMagnets
    write(*,fmt='(A,I0,A)') '        ', nColliding, &
                            ' removed due to collisions' 

    ! Populate the magnet array with bricks that are in-bounds and non-colliding
    if (stell_symm) then
        focus_symm_key = focus_stell_symm
    else if (tor_symm) then
        focus_symm_key = focus_tor_symm
    else
        focus_symm_key = focus_no_symm
    end if
    allocate(magnets(nMagnets))
    n = 0
    do i = 1, nRects
        if (.not. colliding(i)) then
            n = n + 1
            write(label, fmt='(A,I10.10)') 'pm_', n
            magnets(n) = &
                trec_magnet(trecs(i), M_max, focus_symm_key, &
                              rho_free, axis_free, rho_initial, label)
        end if
    end do

    ! Prepare the output trecs and overlaps array
    if (incl_ovl_trecs) then
        nTrecs_out = nRects
    else
        nTrecs_out = nRects - nColliding
    end if
    allocate(trecs_out(nTrecs_out), colliding_out(nTrecs_out))
    n = 0
    do i = 1, nRects
        if (incl_ovl_trecs .or. (.not. colliding(i))) then
            n = n + 1
            colliding_out(n) = colliding(i)
            trecs_out(n) = trecs(i)
        end if
    end do

    ! Deallocate arrays
    deallocate(trArr, in_bounds, phi_cm, trecs, colliding)

end subroutine build_trec_set



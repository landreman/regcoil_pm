!-------------------------------------------------------------------------------
! Construction of a cbrick magnet set
!-------------------------------------------------------------------------------
subroutine build_cbrick_set(lim_surf, nBricks_out, cbricks_out, &
                            nMagnets, magnets, nColliding, colliding_out)

    use magpie_globals,    &
        only: dp, len_mag_label, &
              cbrick_rmin, cbrick_rmax, cbrick_zmin, cbrick_zmax, cbrick_pmin, &
              cbrick_pmax, cbrick_dr, cbrick_dz, cbrick_nPhi, gap_brr,         &
              gap_brz, gap_brp, gap_lim, cbrick_nTheta, radial_extent, pi,     &
              M_max, stell_symm, tor_symm, focus_no_symm, focus_tor_symm,      &
              focus_stell_symm, rho_initial, rho_free, axis_free,              &
              incl_ovl_cbricks, max_ovl_check_interval, cbrick_ax_init
    use surface_calc,      only: surface, surf_r, surf_z
    use surface_solve,     only: surf_perp_intersect_2d
    use cbrick_properties, only: cbrick, cbrick_grid_params, new_cbrick,       &
                                 corner_in_polygon, corners_from_polygon,      &
                                 estimate_cbrick_theta, set_perpendicular,     &
                                 cbrick_obj_overlap
    use magnet_properties, only: magnet, cbrick_magnet
    use input,             only: to_lowercase
    use output,            only: write_focus_file, write_cbrick_file,          &
                                 write_cbrick_vtk, write_cbrick_wf_vtk

    implicit none

    type(surface), intent(IN) :: lim_surf
    integer, intent(OUT) :: nBricks_out, nMagnets, nColliding
    type(cbrick), dimension(:), allocatable, intent(OUT) :: cbricks_out
    logical, dimension(:), allocatable, intent(OUT) :: colliding_out
    type(cbrick), dimension(:,:,:), allocatable :: cbArr
    logical, dimension(:,:,:), allocatable :: in_bounds
    type(cbrick), dimension(:), allocatable :: cbricks
    logical, dimension(:), allocatable :: colliding
    type(magnet), dimension(:), allocatable, intent(OUT) :: magnets
    real(dp) :: r0_grid, z0_grid, cbrick_dphi, r0, z0, r1, z1, phi0, phi1
    integer :: nr, nz, i, j, k, n, nBricks, focus_symm_key
    real(dp) :: dTheta, phi_back, phi_front, phi_cm, theta_cm
    real(dp), dimension(cbrick_nTheta) :: rSectB, zSectB, rSectF, zSectF
    real(dp), dimension(cbrick_nTheta) :: rSectM, zSectM, theta
    real(dp), dimension(8) :: dists
    character(len=len_mag_label) :: label

    write(*,fmt='(A)') '    Constructing grid of cbricks within allotted volume'

    ! Determine radial and vertical dimensions of the grid of potential bricks
    call cbrick_grid_params(cbrick_dr, cbrick_dz, cbrick_rmin, cbrick_rmax, &
                            cbrick_zmin, cbrick_zmax, r0_grid, z0_grid, nr, nz)

    allocate(cbArr(nr,nz,cbrick_nphi), in_bounds(nr,nz,cbrick_nphi))

    ! Populate the array with potential bricks
    cbrick_dphi = (cbrick_pmax - cbrick_pmin)/real(cbrick_nPhi,dp)
    do i = 1, nr
        r0 = r0_grid + (i-1)*cbrick_dr + 0.5*gap_brr
        r1 = r0_grid +     i*cbrick_dr - 0.5*gap_brr
        do j = 1, nz
            z0 = z0_grid + (j-1)*cbrick_dz + 0.5*gap_brz
            z1 = z0_grid +     j*cbrick_dz - 0.5*gap_brz
            do k = 1, cbrick_nPhi
                phi0 = cbrick_pmin + (k-1)*cbrick_dphi + 0.5*gap_brp
                phi1 = cbrick_pmin +     k*cbrick_dphi - 0.5*gap_brp
                cbArr(i,j,k) = new_cbrick(r0, r1, z0, z1, phi0, phi1)
            end do
        end do
    end do

    write(*,fmt='(A)') '    Identifying bricks within acceptable bounds'

    ! Check whether each brick in the grid is within the acceptable bounds
    nBricks = 0
    dTheta = 2*pi/(cbrick_nTheta-1)
    do i = 1, cbrick_nTheta
        theta(i) = (i-1)*dTheta
    end do
    do k = 1, cbrick_nPhi

        phi_back  = cbArr(1,1,k)%phib
        phi_front = cbArr(1,1,k)%phif
        phi_cm    = cbArr(1,1,k)%phi
        
        ! Obtain xsects of limiting surface at front & back tor. angle of brick
        do i = 1, cbrick_nTheta
            rSectB(i)  = surf_r(lim_surf, theta(i), phi_back)
            zSectB(i)  = surf_z(lim_surf, theta(i), phi_back)
            rSectF(i) = surf_r(lim_surf, theta(i), phi_front)
            zSectF(i) = surf_z(lim_surf, theta(i), phi_front)
            rSectM(i) = surf_r(lim_surf, theta(i), phi_cm)
            zSectM(i) = surf_z(lim_surf, theta(i), phi_cm)
        end do

        do i = 1, nr
            do j = 1, nz

                ! Eliminate bricks with corners inside the limiting surface
                if (corner_in_polygon(cbArr(i,j,k), &
                                      cbrick_nTheta, rSectB, zSectB) .or. &
                    corner_in_polygon(cbArr(i,j,k), &
                                      cbrick_nTheta, rSectF, zSectF)     ) then
                    in_bounds(i,j,k) = .false.

                ! Eliminate bricks with corners too close/far outside lim. surf.
                else
                    dists(1:4) = &
                        corners_from_polygon(cbArr(i,j,k), cbrick_nTheta, &
                                             rSectB, zSectB)
                    dists(5:8) = &
                        corners_from_polygon(cbArr(i,j,k), cbrick_nTheta, &
                                             rSectF, zSectF)
                    if (minval(dists) > gap_lim .and. &
                        minval(dists) < radial_extent) then
                        in_bounds(i,j,k) = .true.

                        ! For bricks to be retained, set axis perp. to surface
                        ! if applicable
                        if (to_lowercase(cbrick_ax_init) == 'normal') then
                            theta_cm = &
                                estimate_cbrick_theta(cbArr(i,j,k), &
                                    cbrick_nTheta, theta, rSectM, zSectM)
                            call set_perpendicular(cbArr(i,j,k), lim_surf, &
                                    dists(1), theta_cm, phi_cm)
                        end if

                        nBricks = nBricks + 1
                    else
                        in_bounds(i,j,k) = .false.
                    end if
                end if

            end do
        end do

    end do

    write(*,fmt='(A,I0,A)') '        ', nBricks, ' bricks retained'

    ! Determine which of the in-bounds bricks overlaps ports or other objects
    write(*,fmt='(A)') &
        '    Checking for collisions with ports and other objects'
    allocate(cbricks(nBricks), colliding(nBricks))
    n = 0
    nMagnets = 0
    do i = 1, nr
        do j = 1, nz
            do k = 1, cbrick_nphi
                if (in_bounds(i,j,k)) then
                    n = n + 1
                    cbricks(n) = cbArr(i,j,k)
                    call cbrick_obj_overlap(cbricks(n), &
                                           max_ovl_check_interval, colliding(n))
                    if (.not. colliding(n)) nMagnets = nMagnets + 1
                end if
            end do
        end do
    end do
    nColliding = nBricks - nMagnets

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
    do i = 1, nBricks
        if (.not. colliding(i)) then
            n = n + 1
            write(label, fmt='(A,I10.10)') 'pm_', n
            magnets(n) = &
                cbrick_magnet(cbricks(i), M_max, focus_symm_key, &
                              rho_free, axis_free, rho_initial, label)
        end if
    end do

    ! Prepare the output cbricks and overlaps array
    if (incl_ovl_cbricks) then
        nBricks_out = nBricks
    else
        nBricks_out = nBricks - nColliding
    end if
    allocate(cbricks_out(nBricks_out), colliding_out(nBricks_out))
    n = 0
    do i = 1, nBricks
        if (incl_ovl_cbricks .or. (.not. colliding(i))) then
            n = n + 1
            colliding_out(n) = colliding(i)
            cbricks_out(n) = cbricks(i)
        end if
    end do

    ! Deallocate arrays
    deallocate(cbArr, in_bounds, cbricks, colliding)

end subroutine build_cbrick_set



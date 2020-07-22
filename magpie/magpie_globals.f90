module magpie_globals

    implicit none

    integer, parameter :: dp = 8
    real(dp), parameter :: pi = 3.1415926535897931

    ! Default lengths of certain strings
    integer, parameter :: len_fname     = 200
    integer, parameter :: len_line      = 500
    integer, parameter :: len_mag_label = 50  ! max length of magnet name/ID 
    integer, parameter :: len_suffix    = 20
    integer, parameter :: len_version   = 20

    ! Version number
    character(len=len_version) :: version = '0.0.0'

    ! Sizes of input arrays
    integer, parameter :: ntg = 500         ! Max. num. theta vals., base grid
    integer, parameter :: npg = 500         ! Max. num. phi vals., base grid
    integer, parameter :: nt = ntg - 1      ! Max theta dim. for magnet array
    integer, parameter :: np = npg - 1      ! Max phi dim. for magnet array

    ! Integer unit numbers for different files
    integer, parameter :: unit_in      = 100
    integer, parameter :: unit_surf    = 101
    integer, parameter :: unit_focus   = 102
    integer, parameter :: unit_corners = 103
    integer, parameter :: unit_bg      = 104
    integer, parameter :: unit_port    = 105 
    integer, parameter :: unit_cbrick  = 106
    integer, parameter :: unit_vtk     = 107
    integer, parameter :: unit_wf_vtk  = 108

    integer, parameter :: maxIter = 20

    real(dp), parameter :: dist_tol = 1.0e-7

    ! Parameters for geometric output files
    integer, parameter :: arc_min_pts_per_deg = 1
    integer, parameter :: arcs_per_cbrick = 4

    ! FOCUS-specific parameters
    integer  :: focus_type = 2       ! "coil" type for permanent magnet
    real(dp) :: rho_initial = 0.0    ! initialization value for rho
    logical  :: rho_free = .true.    ! True if rho is optimizable
    logical  :: axis_free = .false.  ! True if moment direction is optimizable
    integer, parameter :: focus_no_symm = 0
    integer, parameter :: focus_tor_symm = 1
    integer, parameter :: focus_stell_symm = 2

    ! Magnet type
    character(len=len_suffix) :: magnet_type = ''  ! qhex, cbrick, or trec

    ! Required gap spacing
    real(dp) :: gap_lim = 0.0
    real(dp) :: gap_qhx = 0.0
    real(dp) :: gap_rad = 0.0

    ! Parameters for curved-brick magnets
    real(dp) :: gap_brr = 0.0     ! Gaps between bricks, radial dimension (m)
    real(dp) :: gap_brz = 0.0     ! Gaps between bricks, vertical dimension (m)
    real(dp) :: gap_brp = 0.0     ! Gaps between bricks, toroidal dimension (m)
    real(dp) :: cbrick_rmin = 0.0 ! Minimum major radius for bricks (m)
    real(dp) :: cbrick_rmax = 0.0 ! Maximum major radius for bricks (m)
    real(dp) :: cbrick_zmin = 0.0 ! Minimum z-coordinate (m) for bricks
    real(dp) :: cbrick_zmax = 0.0 ! Maximum z-coordinate (m) for bricks
    real(dp) :: cbrick_pmin = 0.0 ! Minimum tor angle (rad) for bricks
    real(dp) :: cbrick_pmax = 0.0 ! Maximum tor angle (rad) for bricks
    real(dp) :: cbrick_dr = 0.0   ! Radial dim., grid cell for each brick (m)
    real(dp) :: cbrick_dz = 0.0   ! Vertical dim., grid cell for each brick (m)
    integer :: cbrick_nPhi = 0    ! Number of toroidal partitions in brick grid
    character(len=len_suffix) :: cbrick_ax_init = 'normal'
              ! Initial polarization axis of each brick ('radial' or 'normal')

    ! Parameters for trapezoidally-enclosed rectangular prisms
    real(dp) :: gap_trr = 0.0    ! Gaps between prisms, radial dimension (m)
    real(dp) :: gap_trz = 0.0    ! Gaps between prisms, vertical dimension (m)
    real(dp) :: gap_trp = 0.0    ! Gaps between prisms, toroidal dimension (m)
    real(dp) :: trec_rmin = 0.0  ! Minimum major radius for prisms (m)
    real(dp) :: trec_rmax = 0.0  ! Maximum major radius for prisms (m)
    real(dp) :: trec_zmin = 0.0  ! Minimum z-coordinate (m) for prisms
    real(dp) :: trec_zmax = 0.0  ! Maximum z-coordinate (m) for prisms
    real(dp) :: trec_pmin = 0.0  ! Minimum tor angle (rad) for prisms
    real(dp) :: trec_pmax = 0.0  ! Maximum tor angle (rad) for prisms
    real(dp) :: trec_dr = 0.0    ! Radial dim., grid cell for each prism (m)
    real(dp) :: trec_dz = 0.0    ! Vertical dim., grid cell for each prism (m)
    real(dp) :: trec_lpMin = 0.0 ! Minimium allowable toroidal dimension (m)
    integer :: trec_nPhi = 0     ! Number of toroidal partitions in prism grid
    character(len=len_suffix) :: trec_ax_init = 'normal'
              ! Initial polarization axis of each prism ('radial' or 'normal')

    ! Number of poloidal points for bounds checking
    integer, parameter :: cbrick_nTheta = 360
    integer, parameter :: trec_nTheta   = 360

    ! Precision for finite-interval checking
    real(dp) :: face_adj_interval = 0.001
    real(dp) :: max_ovl_check_interval = 0.005 ! Max spacing between test points
                                               ! qhex faces for checking 
                                               ! overlaps with ports

    ! Parameters for the base grid
    logical :: qhex_cust_vert_theta = .false. ! Use user-input qhex_vert_theta
    logical :: qhex_cust_vert_phi = .false.   ! Use user-input qhex_vert_phi
    logical :: qhex_cust_vert_sep = .false.   ! Use user-input qhex_vert_sep
    real(dp) :: qhex_grid_theta0 = 0.0        ! Start point in theta 
    real(dp) :: qhex_grid_theta1 = 2.0*pi     ! End point in theta
    real(dp) :: qhex_grid_phi0 = 0.0          ! Start point in phi
    real(dp) :: qhex_grid_phi1 = 2.0*pi       ! End point in phi
    real(dp) :: radial_extent= 0.0          ! Radial separation from ref. surf
    real(dp), dimension(ntg,npg) :: qhex_vert_theta ! Custom vertex theta array
    real(dp), dimension(ntg,npg) :: qhex_vert_phi   ! Custom vertex phi array
    real(dp), dimension(ntg,npg) :: qhex_vert_sep   ! Custom vertex separ. array

    ! Parameters for constructing the quadrilaterally-faced hexahedra
    integer :: qhex_nBoxPol = 2 ! # pol. base grid points for boxcar averaging
    integer :: qhex_nBoxTor = 2 ! # tor. base grid points for boxcar averaging


    ! Magnet array parameters
    integer :: nfp                 ! Number of field periods for the config.
    logical :: stell_symm = .true. ! True if stellarator symmetry is assumed
    logical :: tor_symm = .true.   ! True if tor. symm is assumed
    integer :: qhex_nTheta         ! Theta dimension of magnet array
    integer :: qhex_nPhi           ! Phi dim. of magnet array per field period
    real(dp) :: qhex_max_ht = 0.0  ! Maximum height of an individual qhex magnet
    real(dp) :: qhex_max_ht_up = 0.0 ! supersedes qhex_max_ht, upper qhex only
    real(dp) :: qhex_max_ht_lo = 0.0 ! supersedes qhex_max_ht, lower qhex only
    integer :: qhex_nSplitLo = 1    ! # subdivisions to apply to each lower qhex
    integer :: qhex_nSplitUp = 1    ! # subdivisions to apply to each upper qhex
    logical :: qhex_unif_poloidal = .true.  ! If uniform dtheta/ds should be
                                            ! enforced on grid (recommended)
    logical :: qhex_incl_upper    = .false. ! Attempt to fill in concave regions
                                            ! an "upper" layer of magnets
    logical :: qhex_poloidally_closed = .false.
    logical :: qhex_toroidally_closed = .false.
    real(dp) :: qhex_adj_interval = 0.001  ! Interval for dimension adjustment
                                           ! (to resolve overlaps)
    real(dp) :: M_max = 0.0        ! Max. magnetization (dipole momemt/volume)

    ! Surface information
    character(len=len_fname) :: limiting_surf_file = ''
    character(len=len_fname) :: reference_surf_file = ''
    integer :: limiting_surf_phi_sign = 1
    integer :: reference_surf_phi_sign = 1

    ! Output file options
    logical :: focus_file   = .true.      ! Write a focus input file
    logical :: corners_file = .true.      ! Write a corners file (qhex only)
    logical :: base_grid_file = .true.    ! Write a base grid file (qhex only)
    logical :: cbrick_file = .true.       ! Write a cbrick file (cbrick only)
    logical :: magnet_vtk  = .false.      ! Write a vtk file with magnet geom
    logical :: wframe_vtk  = .false.      ! Write a wireframe vtk file
    logical :: incl_err_corners = .false. ! Include erroneous qhexes in
                                          ! the corners file
    logical :: incl_ovl_corners = .false. ! Include qhexes with overlaps
                                          ! in the corners file
    logical :: incl_ovl_cbricks = .false. ! Include cbricks with overlaps in the
                                          ! cbrick file
    logical :: incl_ovl_trecs   = .false. ! Include trecs with overlaps in the
                                          ! cbrick file
    character(len=len_suffix) :: repeat_to_fill = 'none'
    character(len=len_suffix), parameter :: suffix_focus   = '.focus'
    character(len=len_suffix), parameter :: suffix_corners = '_corners.txt'
    character(len=len_suffix), parameter :: suffix_bg      = '_BaseGrid.txt'
    character(len=len_suffix), parameter :: suffix_cbrick  = '_cbricks.txt'
    character(len=len_suffix), parameter :: suffix_mag_vtk = '_magnets.vtk'
    character(len=len_suffix), parameter :: suffix_wf_vtk  = '_wframe.vtk'
    character(len=len_fname) :: output_base_name = 'magpie_out'

    ! Locations of files with port parameter data
    character(len=len_fname) :: dir_ncsx_param = "port_geometry/"
    character(len=len_fname) :: file_param_ncsx_cir = "ncsx_circular_port_params.csv"
    character(len=len_fname) :: file_param_ncsx_nbp = "ncsx_nb_port_params.csv"
    character(len=len_fname) :: file_param_ncsx_p12 = "ncsx_port_12_params.csv"
    character(len=len_fname) :: file_param_ncsx_p04 = "ncsx_port_4_params.csv"
    character(len=len_fname) :: file_param_ncsx_dom = "ncsx_dome_cylindrical_params.csv"
    ! Minimum spacing between magnets and port edges
    real(dp) :: port_gap = 0.0

    ! To be evaluated after checking individual port inclusion variables
    logical :: check_ncsx_ports  = .false.

    ! Ports to be included (none by default)
    logical :: incl_port_ncsx_nb = .false.
    logical :: incl_port_ncsx_02 = .false.
    logical :: incl_port_ncsx_04 = .false.
    logical :: incl_port_ncsx_05 = .false.
    logical :: incl_port_ncsx_06 = .false.
    logical :: incl_port_ncsx_07 = .false.
    logical :: incl_port_ncsx_08 = .false.
    logical :: incl_port_ncsx_09 = .false.
    logical :: incl_port_ncsx_10 = .false.
    logical :: incl_port_ncsx_11 = .false.
    logical :: incl_port_ncsx_12 = .false.
    logical :: incl_port_ncsx_15 = .false.
    logical :: incl_port_ncsx_17 = .false.
    logical :: incl_port_ncsx_18 = .false.
    logical :: incl_port_ncsx_dome = .false.  ! Note: in practice, exclusion of 
                                          ! the dome would imply the exclusion 
                                          ! of ports 17 and 18 as well.
    namelist /magpie/ &
        magnet_type,                 &
        radial_extent,             &
        nfp,                         &
        stell_symm,                  &
        tor_symm,                    &
        M_max,                       &
        rho_free,                    &
        rho_initial,                 &
        axis_free,                   &
        qhex_nTheta,                 &
        qhex_nPhi,                   &
        qhex_unif_poloidal,          &
        qhex_incl_upper,             &
        qhex_nSplitLo,               &
        qhex_nSplitUp,               &
        qhex_poloidally_closed,      &
        qhex_toroidally_closed,      &
        qhex_adj_interval,           &
        qhex_cust_vert_theta,        &
        qhex_cust_vert_phi,          &
        qhex_cust_vert_sep,          &
        qhex_vert_theta,             &
        qhex_vert_phi,               &
        qhex_vert_sep,               &
        qhex_grid_theta0,            &
        qhex_grid_theta1,            &
        qhex_grid_phi0,              &
        qhex_grid_phi1,              &
        qhex_max_ht,                 &
        qhex_max_ht_up,              &
        qhex_max_ht_lo,              &
        qhex_nBoxPol,                &
        qhex_nBoxTor,                &
        limiting_surf_file,          &
        reference_surf_file,         &
        limiting_surf_phi_sign,      &
        reference_surf_phi_sign,     &
        gap_lim,                     &
        gap_qhx,                     &
        gap_rad,                     &
        cbrick_rmin,                 &
        cbrick_rmax,                 &
        cbrick_zmin,                 &
        cbrick_zmax,                 &
        cbrick_pmin,                 &
        cbrick_pmax,                 &
        cbrick_dr,                   &
        cbrick_dz,                   &
        cbrick_nPhi,                 &
        cbrick_ax_init,              &
        gap_brr,                     &
        gap_brz,                     &
        gap_brp,                     &
        trec_rmin,                   &
        trec_rmax,                   &
        trec_zmin,                   &
        trec_zmax,                   &
        trec_pmin,                   &
        trec_pmax,                   &
        trec_dr,                     &
        trec_dz,                     &
        trec_nPhi,                   &
        trec_lpMin,                  &
        trec_ax_init,                &
        gap_trr,                     &
        gap_trz,                     &
        gap_trp,                     &
        output_base_name,            &
        incl_err_corners,            &
        incl_ovl_corners,            &
        incl_ovl_cbricks,            &
        incl_ovl_trecs,              &
        repeat_to_fill,              &
        focus_file,                  &
        corners_file,                &
        base_grid_file,              &
        cbrick_file,                 &
        magnet_vtk,                  &
        wframe_vtk,                  &
        dir_ncsx_param,              &
        file_param_ncsx_cir,         &
        file_param_ncsx_nbp,         &
        file_param_ncsx_p04,         &
        file_param_ncsx_p12,         &
        file_param_ncsx_dom,         &
        port_gap,                    &
        incl_port_ncsx_nb,           &
        incl_port_ncsx_02,           &
        incl_port_ncsx_04,           &
        incl_port_ncsx_05,           &
        incl_port_ncsx_06,           &
        incl_port_ncsx_07,           &
        incl_port_ncsx_08,           &
        incl_port_ncsx_09,           &
        incl_port_ncsx_10,           &
        incl_port_ncsx_11,           &
        incl_port_ncsx_12,           &
        incl_port_ncsx_15,           &
        incl_port_ncsx_17,           &
        incl_port_ncsx_18,           &
        incl_port_ncsx_dome,         &
        max_ovl_check_interval

end module magpie_globals


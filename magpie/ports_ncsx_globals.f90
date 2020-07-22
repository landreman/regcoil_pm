!-------------------------------------------------------------------------------
! ports_ncsx_globals.f90
! 
! Global variables used by the subroutines that determine whether points
! are inside NCSX ports
!-------------------------------------------------------------------------------
module ports_ncsx_globals

    use magpie_globals, only: dp

    implicit none

    real(dp), parameter :: in2m = 0.0254  ! number of meters in one inch
    real(dp), parameter :: pi = 3.141592653589793
    real(dp), parameter :: deg2rad = pi/180. ! degrees in one radian

    ! If false, skip all ncsx_ports subroutines (e.g., if no namelist found)
    logical :: ncsx_ports_on = .false.

    ! Checks if the values in this module have been initialized
    logical :: ncsx_ports_initialized = .false.

    ! Required clearance between magnet locations and port walls (meters)
    real(dp) :: port_gap = 0.

    ! Number of field periods in NCSX
    integer :: nfp_ncsx = 3

    ! Number of ports included in the model
    integer :: nPortsNCSX = 15

    ! Circular port parameters
    integer, parameter          :: nCirPrt = 11        ! number of circ. ports
    integer, dimension(nCirPrt) :: id_cir              ! ID number of port
    real(dp), dimension(nCirPrt) :: xo_cir, yo_cir, zo_cir ! origin of port axis
    real(dp), dimension(nCirPrt) :: ax_cir, ay_cir, az_cir ! axis dir. from origin
    real(dp), dimension(nCirPrt) :: ir_cir                 ! inner radius
    real(dp), dimension(nCirPrt) :: thick_cir       ! wall thickness incl. clearance
    real(dp), dimension(nCirPrt) :: l0_cir, l1_cir  ! port endpts. along axis

    ! Dome parameters (presently modeled as an additional circular port)
    real(dp) :: xo_dom,  yo_dom,  zo_dom   ! origin of dome axis
    real(dp) :: ax_dom,  ay_dom,  az_dom   ! axis dir. from origin
    real(dp) :: xo_dom2, yo_dom2, zo_dom2  ! origin of axis in adjacent half-mod
    real(dp) :: ax_dom2, ay_dom2, az_dom2  ! axis dir. in adjacent half-mod
    real(dp) :: ir_dom                     ! inner radius
    real(dp) :: thick_dom                  ! wall thickness
    real(dp) :: l0_dom, l1_dom             ! endpoints of dome along axis

    ! NB port parameters
    ! Note: there are 6 "ref. circles" @ rounded edges, but only 1-3, 6 are used
    real(dp) :: yo1_pnb, zo1_pnb ! y, z coords, 1st ref. circle (at midplane)
    real(dp) :: yo2_pnb, zo2_pnb ! y, z coords, 2nd ref. circle (above midplane)
    real(dp) :: yo3_pnb, zo3_pnb ! y, z coords, 3rd ref. circle (above midplane)
    real(dp) :: yo6_pnb, zo6_pnb ! y, z coords, 6th ref. circle (below midplane)
    real(dp) :: or1sq_pnb      ! square of outer radius of circle 1 incl. clear.
    real(dp) :: or2sq_pnb      ! square of outer radius of circle 2 incl. clear.
    real(dp) :: or3sq_pnb      ! square of outer radius of circle 3 incl. clear.
    real(dp) :: or6sq_pnb      ! square of outer radius of circle 6 incl. clear.
    real(dp) :: l0_pnb, l1_pnb ! dists of inner & outer endpoints from z-axis
    real(dp) :: zMax_pnb           ! maximum extent in z of the port
    real(dp) :: yMax_pnb           ! maximum extent in y of the port
    real(dp) :: yTop_pnb, zTop_pnb ! intersect btwn top circle and upper line seg
    real(dp) ::  yUp_pnb,  zUp_pnb ! intersect btwn mid circle and upper line seg
    real(dp) ::  yLo_pnb,  zLo_pnb ! intersect btwn mid circle and lower line seg
    real(dp) :: yBot_pnb, zBot_pnb ! intersect btwn bot circle and lower line seg
    real(dp) :: slopeUp_pnb        ! slope, line seg connecting top to mid circle
    real(dp) :: slopeLo_pnb        ! slope, line seg connecting mid to bot circle

    ! Port 4 parameters 
    ! Note: "a" refers to the CW-facing side; "b" refers to the CCW-facing side
    real(dp) :: phi_p04     ! azimuthal angle of the port's reference axis
    real(dp) :: theta_a_p04 ! angle btwn. a-side inner wall and ref. axis
    real(dp) :: theta_b_p04 ! neg. angle btwn. b-side inner wall and ref. axis
    real(dp) :: w1a_p04     ! distance btwn. axis and a-side wall @ narrow xsect 
    real(dp) :: w1b_p04     ! distance btwn. axis and b-side wall @ narrow xsect
    real(dp) :: w2a_p04     ! distance btwn. axis and a-side wall @ port entrance
    real(dp) :: w2b_p04     ! distance btwn. axis and b-side wall @ port entrance
    real(dp) :: slope_a_p04 ! slope of outer wall segment on a-side
    real(dp) :: slope_b_p04 ! slope of outer wall segment on b-side
    real(dp) :: h_p04       ! height (z-dimension) of port
    real(dp) :: zMax_p04    ! maximum value of z attained (0.5*h + thickness)
    real(dp) :: z_circ_p04  ! height of centers of circ. corners above midplane
    real(dp) :: r_circ_p04  ! radius of rounded corners (in plane perp. to axis)
    real(dp) :: orc_p04     ! outer radius of rounded corners
    real(dp) :: orc2_p04       ! outer radius^2, rounded corners 
    real(dp) :: l0_p04, l1_p04 ! distances of inner and outer endpoints from z-axis
    real(dp) :: l_narrow_p04 ! distance of narrow cross-section from z-axis
    real(dp) :: thick_p04    ! wall thickness including clearance
    real(dp) :: thick_ai_p04 ! wall thickness perp. to port axis, a side, inner
    real(dp) :: thick_ao_p04 ! wall thickness perp. to port axis, a side, outer
    real(dp) :: thick_bi_p04 ! wall thickness perp. to port axis, b side, inner
    real(dp) :: thick_bo_p04 ! wall thickness perp. to port axis, b side, outer
    real(dp) :: rhat_x_p04   ! unit vector along the port axis, x component
    real(dp) :: rhat_y_p04  ! unit vector along the port axis, y component
    real(dp) :: phat_x_p04  ! unit vector perpendicular to port axis, x component
    real(dp) :: phat_y_p04  ! unit vector perpendicular to port axis, y component

    ! Port 12 parameters
    real(dp) :: xMin_p12    ! Minimum possible value of the x coordinate
    real(dp) :: xMax_p12    ! Maximum possible value of the x coordinate
    real(dp) :: yMax_p12    ! Maximum possible value of the y coordinate
    real(dp) :: xo1_p12     ! x coord. of the center of the inboard circle
    real(dp) :: xo2_p12     ! x coord. of the center of the outboard circle
    real(dp) :: xInb_p12    ! x coord., intersect. of line seg w/ inboard circ.
    real(dp) :: yInb_p12    ! y coord., intersect. of line seg w/ inboard circ.
    real(dp) :: xOut_p12    ! x coord., intersect. of line seg w/ outboard circ.
    real(dp) :: yOut_p12    ! y coord., intersect. of line seg w/ outboard circ.
    real(dp) :: or1_p12     ! outer radius of the inboard circle incl. clearance
    real(dp) :: or2_p12     ! outer radius of the outboard circle incl. clearance
    real(dp) :: or1sq_p12      ! or1_p12**2
    real(dp) :: or2sq_p12      ! or2_p12**2
    real(dp) :: l0_p12, l1_p12 ! port endpoints (as distances from the xy plane)
    real(dp) :: slope_p12      ! slope of line segment

    logical, dimension(nCirPrt) :: incl_circular_ports


end module ports_ncsx_globals


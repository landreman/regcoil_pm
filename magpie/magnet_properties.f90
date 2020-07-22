!-------------------------------------------------------------------------------
! magnets
!
! Defines the magnet structure and subroutines for initializing magnets
!-------------------------------------------------------------------------------
module magnet_properties

use magpie_globals, only: dp, len_mag_label

implicit none

type magnet
    real(dp) :: ox, oy, oz    ! x, y, z coordinates, center of mass (m)
    real(dp) :: nx, ny, nz    ! x, y, z components, unit polarization axis 
    real(dp) :: phi, theta    ! azimuthal & polar angles of polarization axis
    real(dp) :: moment        ! dipole moment (A*m) (typ. set to max allowable)
    integer  :: symm          ! symmetry designation (for FOCUS)
    logical  :: moment_free   ! true if moment may be varied in optimization
    logical  :: axis_free     ! true if polarization axis may be varied
    real(dp) :: rho_initial   ! initialization value of rho for focus input file
    character(len=len_mag_label) :: label ! name/ID for magnet
end type magnet

contains

!-------------------------------------------------------------------------------
! qhex_magnet(qhex, M, symm, moment_free, axis_free, rho_initial, label)
!
! Generates a magnet object with the geometry of a qhex and an arbitrary 
! magnetization. The dipole moment points in the direction of the "stack normal"
! (snx, sny, snz) of the qhex.
!
! Input parameters:
!     type(quad_hexahedron) :: qhex -> the qhex defining the geometry
!     real(dp) :: M                 -> magnetization (A/m**2)
!     integer  :: symm              -> FOCUS symmetry designation
!     logical  :: moment_free       -> true if moment may vary in optimization
!     real(dp) :: rho_initial       -> initial value of rho to write in file
!     character(len=len_mag_label) :: label -> name/ID for magnet
!-------------------------------------------------------------------------------
type(magnet) function qhex_magnet(qhex, M, symm, moment_free, axis_free, &
                                  rho_initial, label)

    use magpie_globals,  only: len_mag_label
    use qhex_properties, only: quad_hexahedron
    
    implicit none

    type(quad_hexahedron), intent(IN) :: qhex
    real(dp), intent(IN) :: M
    integer, intent(IN)  :: symm
    logical, intent(IN)  :: moment_free, axis_free
    real(dp), intent(IN) :: rho_initial
    character(len=len_mag_label), intent(IN) :: label

    qhex_magnet%ox     = qhex%ox
    qhex_magnet%oy     = qhex%oy
    qhex_magnet%oz     = qhex%oz
    qhex_magnet%nx     = qhex%snx
    qhex_magnet%ny     = qhex%sny
    qhex_magnet%nz     = qhex%snz
    qhex_magnet%moment = qhex%vol * M
    qhex_magnet%phi    = atan2(qhex%sny, qhex%snx)
    qhex_magnet%theta  = atan2(sqrt(qhex%snx**2 + qhex%sny**2), qhex%snz)
    qhex_magnet%symm   = symm
    qhex_magnet%moment_free = moment_free
    qhex_magnet%axis_free   = axis_free
    qhex_magnet%rho_initial = rho_initial
    qhex_magnet%label  = trim(label)

end function qhex_magnet

!-------------------------------------------------------------------------------
! cbrick_magnet(cb, M, symm, moment_free, axis_free, rho_initial, label)
!
! Populates a magnet structure based on the geometric data contained in a
! cbrick structure.
!
! Input parameters
!     type(cbrick) :: cb         -> the cbrick structure defining the geometry
!     real(dp) :: M              -> maximum allowable magnetization (A/m)
!     integer :: symm            -> FOCUS/FAMUS symmetry designation
!     logical :: moment_free     -> True if the moment may vary in optimization
!     logical :: axis_free       -> True if the axis may vary in optimization
!     real(dp) :: rho_initial    -> Initial value of rho to write in file
!     character(len=len_mag_label) :: label  -> name/ID for the magnet
!
! Returns
!     type(magnet) :: mag_out
!-------------------------------------------------------------------------------
function cbrick_magnet(cb, M, symm, moment_free, axis_free, rho_initial, &
                       label) result(mag_out)

    use cbrick_properties, only: cbrick

    implicit none

    type(cbrick), intent(IN) :: cb
    real(dp), intent(IN) :: M
    integer, intent(IN) :: symm
    logical, intent(IN) :: moment_free, axis_free
    real(dp), intent(IN) :: rho_initial
    character(len=len_mag_label), intent(IN) :: label
    type(magnet) :: mag_out

    mag_out%ox = cb%cmx
    mag_out%oy = cb%cmy
    mag_out%oz = cb%cmz
    mag_out%nx = cb%axx
    mag_out%ny = cb%axy
    mag_out%nz = cb%axz
    mag_out%moment = cb%vol * M
    mag_out%phi = atan2(cb%axy, cb%axx)
    mag_out%theta = atan2(sqrt(cb%axx**2 + cb%axy**2), cb%axz)
    mag_out%symm = symm
    mag_out%moment_free = moment_free
    mag_out%axis_free = axis_free
    mag_out%rho_initial = rho_initial
    mag_out%label = label

end function cbrick_magnet

!-------------------------------------------------------------------------------
! trec_magnet(tr, M, symm, moment_free, axis_free, rho_initial, label)
!
! Populates a magnet structure based on the geometric data contained in a
! trec structure.
!
! Input parameters
!     type(trec) :: tr           -> the trec structure defining the geometry
!     real(dp) :: M              -> maximum allowable magnetization (A/m)
!     integer :: symm            -> FOCUS/FAMUS symmetry designation
!     logical :: moment_free     -> True if the moment may vary in optimization
!     logical :: axis_free       -> True if the axis may vary in optimization
!     real(dp) :: rho_initial    -> Initial value of rho to write in file
!     character(len=len_mag_label) :: label  -> name/ID for the magnet
!
! Returns
!     type(magnet) :: mag_out
!-------------------------------------------------------------------------------
function trec_magnet(tr, M, symm, moment_free, axis_free, rho_initial, &
                     label) result(mag_out)

    use trec_properties, only: trec

    implicit none

    type(trec), intent(IN) :: tr
    real(dp), intent(IN) :: M
    integer, intent(IN) :: symm
    logical, intent(IN) :: moment_free, axis_free
    real(dp), intent(IN) :: rho_initial
    character(len=len_mag_label), intent(IN) :: label
    type(magnet) :: mag_out

    mag_out%ox = tr%cmx
    mag_out%oy = tr%cmy
    mag_out%oz = tr%cmz
    mag_out%nx = tr%axx
    mag_out%ny = tr%axy
    mag_out%nz = tr%axz
    mag_out%moment = tr%vol * M
    mag_out%phi = atan2(tr%axy, tr%axx)
    mag_out%theta = atan2(sqrt(tr%axx**2 + tr%axy**2), tr%axz)
    mag_out%symm = symm
    mag_out%moment_free = moment_free
    mag_out%axis_free = axis_free
    mag_out%rho_initial = rho_initial
    mag_out%label = label

end function trec_magnet

end module magnet_properties


!-------------------------------------------------------------------------------
! qhex_properties.f90
!
! Subroutines for determining properties of quadrilaterally-faced hexahedra 
! Normal vectors (not unit length) for the two triangles composing each face
!
! The vertices and faces are labeled according to the general layout of the 
! hexahedron relative to the vacuum vessel (or other bounding surface that 
! encloses the toroidal plasma). The "top" face is closest to the surface of 
! the vessel, and its outward-pointing normal vector points toward the vessel 
! surface (in the minor radial direction). The "base" face is opposite the top 
! face and is constrained to be parallel with the top face. The "toroidal 
! front" face is oriented such that its outward pointing normal vector is 
! approximately in the toroidal direction. The "toroidal back" face is opposite 
! the toroidal front face but not necessarily parallel; its outward normal 
! points approximately in the negative toroidal direction. Similarly, the 
! poloidal front and back faces have outward normals that point in the positive 
! and negative poloidal directions, respectively.
!
! The vertices on the top and base faces are each numbered 1 through 4 as shown
! in the diagram below. Looking from the base face toward the top face and 
! vessel, the numbers go counter-clockwise about the base and top faces, 
! beginning at the intersection of the toroidal and poloidal back planes.
!
!            Poloidal
!             front
!          4----------3 
!          |   Base   |
! Toroidal |  (view   | Toroidal    
!   back   |  toward  |  front     
!          |  vessel) |
!          1----------2
!            Poloidal
!              back
!-------------------------------------------------------------------------------
module qhex_properties

use magpie_globals, only: dp

implicit none

integer, parameter :: nHexFaces = 6
integer, parameter :: nQhexCorners = 8
integer, parameter :: nQhexFaceProps = 15
integer, parameter :: iTOP = 1, iBASE = 2  ! IDs, top and base face
integer, parameter :: iTB  = 3, iTF   = 4  ! IDs, tor. back and front faces
integer, parameter :: iPB  = 5, iPF   = 6  ! IDs, pol. back and front faces
integer, parameter :: iNX  =  1, iNY  =  2, iNZ  =  3  ! 
integer, parameter :: iCX1 =  4, iCY1 =  5, iCZ1 =  6  ! IDs for each row of the
integer, parameter :: iCX2 =  7, iCY2 =  8, iCZ2 =  9  ! face_parameters array
integer, parameter :: iCX3 = 10, iCY3 = 11, iCZ3 = 12  ! (Normal + 4 Corners)
integer, parameter :: iCX4 = 13, iCY4 = 14, iCZ4 = 15  !

type quad_hexahedron
    real(dp) :: ntbx, ntby, ntbz ! unit normal vector of toroidal back plane
    real(dp) :: ntfx, ntfy, ntfz ! unit normal vector of toroidal front plane
    real(dp) :: npbx, npby, npbz ! unit normal vector of poloidal back plane
    real(dp) :: npfx, npfy, npfz ! unit normal vector of poloidal front plane
    real(dp) :: snx,  sny,  snz  ! stack normal vector (toward vessel)
    real(dp) :: otbx, otby, otbz ! reference point on toroidal back plane
    real(dp) :: otfx, otfy, otfz ! reference point on toroidal front plane
    real(dp) :: opbx, opby, opbz ! reference point on poloidal back plane
    real(dp) :: opfx, opfy, opfz ! reference point on poloidal front plane
    real(dp) :: otx,  oty,  otz  ! reference point on top plane
    real(dp) :: obx,  oby,  obz  ! reference point ("origin") on bottom plane
    real(dp) :: vx, vy, vz       ! intersection of stack normal with vessel
    real(dp) :: ves_theta        ! theta, intersect., stack norm with vessel
    real(dp) :: ves_phi          ! phi, intersection of norm with vessel
    real(dp) :: rx, ry, rz       ! intersection of stack normal witih ref surf
    real(dp) :: ref_theta        ! theta, intersect., stack norm with ref surf
    real(dp) :: ref_phi          ! phi, intersect., stack norm with ref surf
    real(dp) :: ctx1, cty1, ctz1 ! coords, corner 1, top face (b tor b pol)
    real(dp) :: ctx2, cty2, ctz2 ! coords, corner 2, top face (f tor b pol)
    real(dp) :: ctx3, cty3, ctz3 ! coords, corner 3, top face (f tor f pol)
    real(dp) :: ctx4, cty4, ctz4 ! coords, corner 4, top face (b tor f pol)
    real(dp) :: cbx1, cby1, cbz1 ! coords, corner 1, bot face (b tor b pol)
    real(dp) :: cbx2, cby2, cbz2 ! coords, corner 2, bot face (f tor b pol)
    real(dp) :: cbx3, cby3, cbz3 ! coords, corner 3, bot face (f tor f pol)
    real(dp) :: cbx4, cby4, cbz4 ! coords, corner 4, bot face (b tor f pol)
    real(dp) :: height           ! distance between base and top planes
    real(dp) :: vol              ! volume contained within the trapezoid
    real(dp) :: ox, oy, oz       ! coords, centroid (center of mass)
    logical :: err              ! true if geometry is erroneous/undesirable
end type quad_hexahedron

! Pair of upper and lower qhex structures that share a grid point in an array
type qhex_pair
    logical :: incl_upper
    type(quad_hexahedron) :: lo
    type(quad_hexahedron) :: up
end type qhex_pair

! Wrapper data type for storing the set of qhex structures for a configuration
type qhex_array
    logical :: incl_upper        ! True if upper qhexes are to be included
    integer :: nLower, nUpper    ! Total # of lower and upper qhexes
    integer :: nTotal            ! Total number of qhex in the config
    integer :: nTheta, nPhi      ! Num of lower/upper pairs in theta, phi dim
    integer :: nErrLo, nErrUp    ! Number of erroneous lower and upper qhex
    integer :: nErr              ! Total number of erroneous qhex
    type(qhex_pair), allocatable :: pairs(:,:)   ! Array of qhex pairs
end type qhex_array

contains

!-------------------------------------------------------------------------------
! qhex_corners(qh)
!
! Populates the coordinates of the eight corners of a quadrilaterally-faced
! hexahedron with defined normal vectors and reference points for each of its 
! six faces.
!-------------------------------------------------------------------------------
subroutine qhex_corners(qh)

    use geometry, only: three_plane_intersect

    implicit none

    type(quad_hexahedron) :: qh
    
    ! Top face, corner 1: toroidal back face and poloidal back face
    call three_plane_intersect(                                                &
                         qh%otx,  qh%oty,  qh%otz,  qh%snx,  qh%sny,  qh%snz,  &
                         qh%otbx, qh%otby, qh%otbz, qh%ntbx, qh%ntby, qh%ntbz, &
                         qh%opbx, qh%opby, qh%opbz, qh%npbx, qh%npby, qh%npbz, &
                         qh%ctx1, qh%cty1, qh%ctz1)

    ! Top face, corner 2: toroidal front face and poloidal back face
    call three_plane_intersect(                                                &
                         qh%otx,  qh%oty,  qh%otz,  qh%snx,  qh%sny,  qh%snz,  &
                         qh%otfx, qh%otfy, qh%otfz, qh%ntfx, qh%ntfy, qh%ntfz, &
                         qh%opbx, qh%opby, qh%opbz, qh%npbx, qh%npby, qh%npbz, &
                         qh%ctx2, qh%cty2, qh%ctz2)

    ! Top face, corner 3: toroidal front face and poloidal front face
    call three_plane_intersect(                                                &
                         qh%otx,  qh%oty,  qh%otz,  qh%snx,  qh%sny,  qh%snz,  &
                         qh%otfx, qh%otfy, qh%otfz, qh%ntfx, qh%ntfy, qh%ntfz, &
                         qh%opfx, qh%opfy, qh%opfz, qh%npfx, qh%npfy, qh%npfz, &
                         qh%ctx3, qh%cty3, qh%ctz3)

    ! Top face, corner 4: toroidal back face and poloidal front face
    call three_plane_intersect(                                                &
                         qh%otx,  qh%oty,  qh%otz,  qh%snx,  qh%sny,  qh%snz,  &
                         qh%otbx, qh%otby, qh%otbz, qh%ntbx, qh%ntby, qh%ntbz, &
                         qh%opfx, qh%opfy, qh%opfz, qh%npfx, qh%npfy, qh%npfz, &
                         qh%ctx4, qh%cty4, qh%ctz4)

    ! Bottom face, corner 1: toroidal back face and poloidal back face
    call three_plane_intersect(                                                &
                         qh%obx,  qh%oby,  qh%obz,  qh%snx,  qh%sny,  qh%snz,  &
                         qh%otbx, qh%otby, qh%otbz, qh%ntbx, qh%ntby, qh%ntbz, &
                         qh%opbx, qh%opby, qh%opbz, qh%npbx, qh%npby, qh%npbz, &
                         qh%cbx1, qh%cby1, qh%cbz1)

    ! Top face, corner 2: toroidal front face and poloidal back face
    call three_plane_intersect(                                                &
                         qh%obx,  qh%oby,  qh%obz,  qh%snx,  qh%sny,  qh%snz,  &
                         qh%otfx, qh%otfy, qh%otfz, qh%ntfx, qh%ntfy, qh%ntfz, &
                         qh%opbx, qh%opby, qh%opbz, qh%npbx, qh%npby, qh%npbz, &
                         qh%cbx2, qh%cby2, qh%cbz2)

    ! Top face, corner 3: toroidal front face and poloidal front face
    call three_plane_intersect(                                                &
                         qh%obx,  qh%oby,  qh%obz,  qh%snx,  qh%sny,  qh%snz,  &
                         qh%otfx, qh%otfy, qh%otfz, qh%ntfx, qh%ntfy, qh%ntfz, &
                         qh%opfx, qh%opfy, qh%opfz, qh%npfx, qh%npfy, qh%npfz, &
                         qh%cbx3, qh%cby3, qh%cbz3)

    ! Top face, corner 4: toroidal back face and poloidal front face
    call three_plane_intersect(                                                &
                         qh%obx,  qh%oby,  qh%obz,  qh%snx,  qh%sny,  qh%snz,  &
                         qh%otbx, qh%otby, qh%otbz, qh%ntbx, qh%ntby, qh%ntbz, &
                         qh%opfx, qh%opfy, qh%opfz, qh%npfx, qh%npfy, qh%npfz, &
                         qh%cbx4, qh%cby4, qh%cbz4)

    ! Update the height parameter
    qh%height = (qh%otx - qh%obx)*qh%snx + (qh%oty - qh%oby)*qh%sny &
                   + (qh%otz - qh%obz)*qh%snz

end subroutine qhex_corners

!-------------------------------------------------------------------------------
! qhex_volume_ctr(tp)
!
! Calculates the volume enclosed in a trapezoid, as well as the centroid.
! 
! The trapezoid structure tp is taken as an input. All 8 corner coordinates are
! assumed to be defined within the structure prior to the call to this 
! subroutine. During execution, the following parameters of tp will be updated:
!     tp%vol: the enclosed volume
!     tp%ox:  the x coordinate of the centroid
!     tp%oy:  the y coordinate of the centroid
!     tp%oz:  the z coordinate of the centroid
! 
! Based on the formulas for volume and centroid in lecture notes 
! by R. Nuernberg, Imperial College of London, accessed 2019-11-03 from
! wwwf.imperial.ac.uk/~rn/centroid.pdf
!-------------------------------------------------------------------------------
subroutine qhex_volume_ctr(qh)

    use geometry, only: cross_prod

    implicit none

    type(quad_hexahedron) :: qh
    integer :: i
    real(dp) :: ntx1, nty1, ntz1, ntx2, nty2, ntz2
    real(dp) :: nbx1, nby1, nbz1, nbx2, nby2, nbz2
    real(dp) :: ntbx1, ntby1, ntbz1, ntbx2, ntby2, ntbz2
    real(dp) :: ntfx1, ntfy1, ntfz1, ntfx2, ntfy2, ntfz2
    real(dp) :: npbx1, npby1, npbz1, npbx2, npby2, npbz2
    real(dp) :: npfx1, npfy1, npfz1, npfx2, npfy2, npfz2
    real(dp) :: trip_prod(12), quadsum_x(12), quadsum_y(12), quadsum_z(12)
    
    ! Base face, triangle 1: base2-base1 x base4-base1
    call cross_prod( qh%cbx2 - qh%cbx1, qh%cby2 - qh%cby1, qh%cbz2 - qh%cbz1, &
                     qh%cbx4 - qh%cbx1, qh%cby4 - qh%cby1, qh%cbz4 - qh%cbz1, &
                     nbx1, nby1, nbz1 )
    trip_prod(1) = qh%cbx1 *  nbx1 + qh%cby1 *  nby1 + qh%cbz1 *  nbz1
    quadsum_x(1) = nbx1 * ((qh%cbx1 + qh%cbx2)**2 + (qh%cbx1 + qh%cbx4)**2 &
                                                  + (qh%cbx2 + qh%cbx4)**2  )
    quadsum_y(1) = nby1 * ((qh%cby1 + qh%cby2)**2 + (qh%cby1 + qh%cby4)**2 &
                                                  + (qh%cby2 + qh%cby4)**2  )
    quadsum_z(1) = nbz1 * ((qh%cbz1 + qh%cbz2)**2 + (qh%cbz1 + qh%cbz4)**2 &
                                                  + (qh%cbz2 + qh%cbz4)**2  )

    ! Base face, triangle 2: base4-base3 x base2-base3
    call cross_prod( qh%cbx4 - qh%cbx3, qh%cby4 - qh%cby3, qh%cbz4 - qh%cbz3, &
                     qh%cbx2 - qh%cbx3, qh%cby2 - qh%cby3, qh%cbz2 - qh%cbz3, &
                     nbx2, nby2, nbz2 )
    trip_prod(2) = qh%cbx3 *  nbx2 + qh%cby3 *  nby2 + qh%cbz3 *  nbz2
    quadsum_x(2) = nbx2 * ((qh%cbx2 + qh%cbx3)**2 + (qh%cbx2 + qh%cbx4)**2 &
                                                  + (qh%cbx3 + qh%cbx4)**2  )
    quadsum_y(2) = nby2 * ((qh%cby2 + qh%cby3)**2 + (qh%cby2 + qh%cby4)**2 &
                                                  + (qh%cby3 + qh%cby4)**2  )
    quadsum_z(2) = nbz2 * ((qh%cbz2 + qh%cbz3)**2 + (qh%cbz2 + qh%cbz4)**2 &
                                                  + (qh%cbz3 + qh%cbz4)**2  )

    ! Poloidal back face, triangle 1: top1-base1 x base2-base1
    call cross_prod( qh%ctx1 - qh%cbx1, qh%cty1 - qh%cby1, qh%ctz1 - qh%cbz1, &
                     qh%cbx2 - qh%cbx1, qh%cby2 - qh%cby1, qh%cbz2 - qh%cbz1, &
                     npbx1, npby1, npbz1 )
    trip_prod(3) = qh%cbx1 * npbx1 + qh%cby1 * npby1 + qh%cbz1 * npbz1
    quadsum_x(3) = npbx1 * ((qh%ctx1 + qh%cbx2)**2 + (qh%ctx1 + qh%cbx1)**2 &
                                                   + (qh%cbx2 + qh%cbx1)**2  )
    quadsum_y(3) = npby1 * ((qh%cty1 + qh%cby2)**2 + (qh%cty1 + qh%cby1)**2 &
                                                   + (qh%cby2 + qh%cby1)**2  )
    quadsum_z(3) = npbz1 * ((qh%ctz1 + qh%cbz2)**2 + (qh%ctz1 + qh%cbz1)**2 &
                                                   + (qh%cbz2 + qh%cbz1)**2  )

    ! Poloidal back face, triangle 2: base2-top2 x top1-top2
    call cross_prod( qh%cbx2 - qh%ctx2, qh%cby2 - qh%cty2, qh%cbz2 - qh%ctz2, &
                     qh%ctx1 - qh%ctx2, qh%cty1 - qh%cty2, qh%ctz1 - qh%ctz2, &
                     npbx2, npby2, npbz2 )
    trip_prod(4) = qh%ctx2 * npbx2 + qh%cty2 * npby2 + qh%ctz2 * npbz2
    quadsum_x(4) = npbx2 * ((qh%ctx1 + qh%cbx2)**2 + (qh%ctx1 + qh%ctx2)**2 &
                                                   + (qh%cbx2 + qh%ctx2)**2  )
    quadsum_y(4) = npby2 * ((qh%cty1 + qh%cby2)**2 + (qh%cty1 + qh%cty2)**2 &
                                                   + (qh%cby2 + qh%cty2)**2  )
    quadsum_z(4) = npbz2 * ((qh%ctz1 + qh%cbz2)**2 + (qh%ctz1 + qh%ctz2)**2 &
                                                   + (qh%cbz2 + qh%ctz2)**2  )

    ! Poloidal front face, triangle 1: top3-base3 x base4-base3
    call cross_prod( qh%ctx3 - qh%cbx3, qh%cty3 - qh%cby3, qh%ctz3 - qh%cbz3, &
                     qh%cbx4 - qh%cbx3, qh%cby4 - qh%cby3, qh%cbz4 - qh%cbz3, &
                     npfx1, npfy1, npfz1 )
    trip_prod(5) = qh%cbx3 * npfx1 + qh%cby3 * npfy1 + qh%cbz3 * npfz1
    quadsum_x(5) = npfx1 * ((qh%ctx3 + qh%cbx3)**2 + (qh%ctx3 + qh%cbx4)**2 &
                                                   + (qh%cbx3 + qh%cbx4)**2  )
    quadsum_y(5) = npfy1 * ((qh%cty3 + qh%cby3)**2 + (qh%cty3 + qh%cby4)**2 &
                                                   + (qh%cby3 + qh%cby4)**2  )
    quadsum_z(5) = npfz1 * ((qh%ctz3 + qh%cbz3)**2 + (qh%ctz3 + qh%cbz4)**2 &
                                                   + (qh%cbz3 + qh%cbz4)**2  )

    ! Poloidal front face, triangle 2: base4-top4 x top3-top4
    call cross_prod( qh%cbx4 - qh%ctx4, qh%cby4 - qh%cty4, qh%cbz4 - qh%ctz4, &
                     qh%ctx3 - qh%ctx4, qh%cty3 - qh%cty4, qh%ctz3 - qh%ctz4, &
                     npfx2, npfy2, npfz2 )
    trip_prod(6) = qh%ctx4 * npfx2 + qh%cty4 * npfy2 + qh%ctz4 * npfz2
    quadsum_x(6) = npfx2 * ((qh%ctx3 + qh%ctx4)**2 + (qh%ctx3 + qh%cbx4)**2 &
                                                   + (qh%ctx4 + qh%cbx4)**2  )
    quadsum_y(6) = npfy2 * ((qh%cty3 + qh%cty4)**2 + (qh%cty3 + qh%cby4)**2 &
                                                   + (qh%cty4 + qh%cby4)**2  )
    quadsum_z(6) = npfz2 * ((qh%ctz3 + qh%ctz4)**2 + (qh%ctz3 + qh%cbz4)**2 &
                                                   + (qh%ctz4 + qh%cbz4)**2  )

    ! Toroidal back face, triangle 1: base4-base1 x top1-base1
    call cross_prod( qh%cbx4 - qh%cbx1, qh%cby4 - qh%cby1, qh%cbz4 - qh%cbz1, &
                     qh%ctx1 - qh%cbx1, qh%cty1 - qh%cby1, qh%ctz1 - qh%cbz1, &
                     ntbx1, ntby1, ntbz1 )
    trip_prod(7) = qh%cbx1 * ntbx1 + qh%cby1 * ntby1 + qh%cbz1 * ntbz1
    quadsum_x(7) = ntbx1 * ((qh%cbx1 + qh%ctx1)**2 + (qh%cbx1 + qh%cbx4)**2 &
                                                   + (qh%ctx1 + qh%cbx4)**2  )
    quadsum_y(7) = ntby1 * ((qh%cby1 + qh%cty1)**2 + (qh%cby1 + qh%cby4)**2 &
                                                   + (qh%cty1 + qh%cby4)**2  )
    quadsum_z(7) = ntbz1 * ((qh%cbz1 + qh%ctz1)**2 + (qh%cbz1 + qh%cbz4)**2 &
                                                   + (qh%ctz1 + qh%cbz4)**2  )

    ! Toroidal back face, triangle 2: top1-top4 x base4-top4 
    call cross_prod( qh%ctx1 - qh%ctx4, qh%cty1 - qh%cty4, qh%ctz1 - qh%ctz4, &
                     qh%cbx4 - qh%ctx4, qh%cby4 - qh%cty4, qh%cbz4 - qh%ctz4, &
                     ntbx2, ntby2, ntbz2 )
    trip_prod(8) = qh%ctx4 * ntbx2 + qh%cty4 * ntby2 + qh%ctz4 * ntbz2
    quadsum_x(8) = ntbx2 * ((qh%ctx1 + qh%ctx4)**2 + (qh%ctx1 + qh%cbx4)**2 &
                                                   + (qh%ctx4 + qh%cbx4)**2  )
    quadsum_y(8) = ntby2 * ((qh%cty1 + qh%cty4)**2 + (qh%cty1 + qh%cby4)**2 &
                                                   + (qh%cty4 + qh%cby4)**2  )
    quadsum_z(8) = ntbz2 * ((qh%ctz1 + qh%ctz4)**2 + (qh%ctz1 + qh%cbz4)**2 &
                                                   + (qh%ctz4 + qh%cbz4)**2  )

    ! Toroidal front face, triangle 1: top2-base2 x base3-base2
    call cross_prod( qh%ctx2 - qh%cbx2, qh%cty2 - qh%cby2, qh%ctz2 - qh%cbz2, &
                     qh%cbx3 - qh%cbx2, qh%cby3 - qh%cby2, qh%cbz3 - qh%cbz2, &
                     ntfx1, ntfy1, ntfz1 )
    trip_prod(9) = qh%cbx2 * ntfx1 + qh%cby2 * ntfy1 + qh%cbz2 * ntfz1
    quadsum_x(9) = ntfx1 * ((qh%ctx2 + qh%cbx2)**2 + (qh%ctx2 + qh%cbx3)**2 &
                                                   + (qh%cbx2 + qh%cbx3)**2  )
    quadsum_y(9) = ntfy1 * ((qh%cty2 + qh%cby2)**2 + (qh%cty2 + qh%cby3)**2 &
                                                   + (qh%cby2 + qh%cby3)**2  )
    quadsum_z(9) = ntfz1 * ((qh%ctz2 + qh%cbz2)**2 + (qh%ctz2 + qh%cbz3)**2 &
                                                   + (qh%cbz2 + qh%cbz3)**2  )

    ! Toroidal front face, triangle 2: base3-top3 x top2-top3
    call cross_prod( qh%cbx3 - qh%ctx3, qh%cby3 - qh%cty3, qh%cbz3 - qh%ctz3, &
                     qh%ctx2 - qh%ctx3, qh%cty2 - qh%cty3, qh%ctz2 - qh%ctz3, &
                     ntfx2, ntfy2, ntfz2 )
    trip_prod(10) = qh%ctx3 * ntfx2 + qh%cty3 * ntfy2 + qh%ctz3 * ntfz2
    quadsum_x(10) = ntfx2 * ((qh%cbx3 + qh%ctx2)**2 + (qh%cbx3 + qh%ctx3)**2 &
                                                    + (qh%ctx2 + qh%ctx3)**2  )
    quadsum_y(10) = ntfy2 * ((qh%cby3 + qh%cty2)**2 + (qh%cby3 + qh%cty3)**2 &
                                                    + (qh%cty2 + qh%cty3)**2  )
    quadsum_z(10) = ntfz2 * ((qh%cbz3 + qh%ctz2)**2 + (qh%cbz3 + qh%ctz3)**2 &
                                                    + (qh%ctz2 + qh%ctz3)**2  )

    ! Top face, triangle 1: top4-top1 x top2-top1
    call cross_prod( qh%ctx4 - qh%ctx1, qh%cty4 - qh%cty1, qh%ctz4 - qh%ctz1, &
                     qh%ctx2 - qh%ctx1, qh%cty2 - qh%cty1, qh%ctz2 - qh%ctz1, &
                     ntx1, nty1, ntz1 )
    trip_prod(11) = qh%ctx1 *  ntx1 + qh%cty1 *  nty1 + qh%ctz1 *  ntz1
    quadsum_x(11) = ntx1 * ((qh%ctx1 + qh%ctx2)**2 + (qh%ctx1 + qh%ctx4)**2 &
                                                   + (qh%ctx2 + qh%ctx4)**2  )
    quadsum_y(11) = nty1 * ((qh%cty1 + qh%cty2)**2 + (qh%cty1 + qh%cty4)**2 &
                                                   + (qh%cty2 + qh%cty4)**2  )
    quadsum_z(11) = ntz1 * ((qh%ctz1 + qh%ctz2)**2 + (qh%ctz1 + qh%ctz4)**2 &
                                                   + (qh%ctz2 + qh%ctz4)**2  )

    ! Top face, triangle 2: top2-top3 x top4-top3
    call cross_prod( qh%ctx2 - qh%ctx3, qh%cty2 - qh%cty3, qh%ctz2 - qh%ctz3, &
                     qh%ctx4 - qh%ctx3, qh%cty4 - qh%cty3, qh%ctz4 - qh%ctz3, &
                     ntx2, nty2, ntz2 )
    trip_prod(12) = qh%ctx3 *  ntx2 + qh%cty3 *  nty2 + qh%ctz3 *  ntz2
    quadsum_x(12) = ntx2 * ((qh%ctx4 + qh%ctx2)**2 + (qh%ctx4 + qh%ctx3)**2 &
                                                   + (qh%ctx2 + qh%ctx3)**2  )
    quadsum_y(12) = nty2 * ((qh%cty4 + qh%cty2)**2 + (qh%cty4 + qh%cty3)**2 &
                                                   + (qh%cty2 + qh%cty3)**2  )
    quadsum_z(12) = ntz2 * ((qh%ctz4 + qh%ctz2)**2 + (qh%ctz4 + qh%ctz3)**2 &
                                                   + (qh%ctz2 + qh%ctz3)**2  )

    ! For the volume, add the contributions from each cross product, 
    ! dotted with a vector from the origin to one of the vertices
    qh%vol = 0.0 
    do i = 1, 12
        !write(*, fmt='(A, I2, A, ES11.4)') &
        !      'triple_product ', i, ': ', trip_prod(i)
        qh%vol = qh%vol + trip_prod(i)
    end do
    qh%vol = qh%vol / 6.0

    ! Calculate the centroid based on the quadrature sum contributions
    qh%ox = 0.0
    qh%oy = 0.0
    qh%oz = 0.0
    do i = 1, 12
        !write(*, fmt='(A, I2, A, ES11.4)') &
        !      'quadsum ', i, ': ', quadsum_x(i)
        qh%ox = qh%ox + quadsum_x(i)
        qh%oy = qh%oy + quadsum_y(i)
        qh%oz = qh%oz + quadsum_z(i)
    end do
    qh%ox = qh%ox / (qh%vol * 48.0)
    qh%oy = qh%oy / (qh%vol * 48.0)
    qh%oz = qh%oz / (qh%vol * 48.0)

end subroutine qhex_volume_ctr

!-------------------------------------------------------------------------------
! qhex_face_parameters(qh, fp)
!
! Collects normal vectors and vertex coordinates of each of the six faces of a
! a quadrilateral hexahedron into an ordered 2D array to enable iteration 
! through the hexahedron's faces.
!
! Input parameter:
!     type(quad_hexahedron) :: qh -> quad_hexahedron struct for which the array 
!                                    is to be made
!
! Output parameter:
!     real(dp), dimension(6,15) :: fp
!                             ->  array in which each row corresponds to a face
!                                 of the hexahedron, and the columns contain the
!                                 x, y, and z components of the normal vector
!                                 and the coordinates of the vertices
!-------------------------------------------------------------------------------
subroutine qhex_face_parameters(qh, fp)

    implicit none

    type(quad_hexahedron), intent(IN) :: qh
    real(dp), dimension(nHexFaces,nQhexFaceProps), intent(OUT) :: fp

    ! Top face: row 1 of the array
    fp(iTOP,iNX:iNZ)   = (/ qh%snx,  qh%sny,  qh%snz  /)
    fp(iTOP,iCX1:iCZ1) = (/ qh%ctx1, qh%cty1, qh%ctz1 /)
    fp(iTOP,iCX2:iCZ2) = (/ qh%ctx2, qh%cty2, qh%ctz2 /)
    fp(iTOP,iCX3:iCZ3) = (/ qh%ctx3, qh%cty3, qh%ctz3 /)
    fp(iTOP,iCX4:iCZ4) = (/ qh%ctx4, qh%cty4, qh%ctz4 /)

    ! Base/"bottom" face: row 2 of the array
    fp(iBASE,iNX:iNZ)   = (/ qh%snx,  qh%sny,  qh%snz  /)
    fp(iBASE,iCX1:iCZ1) = (/ qh%cbx1, qh%cby1, qh%cbz1 /)
    fp(iBASE,iCX2:iCZ2) = (/ qh%cbx2, qh%cby2, qh%cbz2 /)
    fp(iBASE,iCX3:iCZ3) = (/ qh%cbx3, qh%cby3, qh%cbz3 /)
    fp(iBASE,iCX4:iCZ4) = (/ qh%cbx4, qh%cby4, qh%cbz4 /)

    ! Toroidal back face: row 3 of the array
    fp(iTB,iNX:iNZ)   = (/ qh%ntbx, qh%ntby, qh%ntbz /)
    fp(iTB,iCX1:iCZ1) = (/ qh%cbx1, qh%cby1, qh%cbz1 /)
    fp(iTB,iCX2:iCZ2) = (/ qh%cbx4, qh%cby4, qh%cbz4 /)
    fp(iTB,iCX3:iCZ3) = (/ qh%ctx4, qh%cty4, qh%ctz4 /)
    fp(iTB,iCX4:iCZ4) = (/ qh%ctx1, qh%cty1, qh%ctz1 /)

    ! Toroidal front face: row 4 of the array
    fp(iTF,iNX:iNZ)   = (/ qh%ntfx, qh%ntfy, qh%ntfz /)
    fp(iTF,iCX1:iCZ1) = (/ qh%cbx3, qh%cby3, qh%cbz3 /)
    fp(iTF,iCX2:iCZ2) = (/ qh%cbx2, qh%cby2, qh%cbz2 /)
    fp(iTF,iCX3:iCZ3) = (/ qh%ctx2, qh%cty2, qh%ctz2 /)
    fp(iTF,iCX4:iCZ4) = (/ qh%ctx3, qh%cty3, qh%ctz3 /)

    ! Poloidal back face: row 5 of the array
    fp(iPB,iNX:iNZ)   = (/ qh%npbx, qh%npby, qh%npbz /)
    fp(iPB,iCX1:iCZ1) = (/ qh%cbx2, qh%cby2, qh%cbz2 /)
    fp(iPB,iCX2:iCZ2) = (/ qh%cbx1, qh%cby1, qh%cbz1 /)
    fp(iPB,iCX3:iCZ3) = (/ qh%ctx1, qh%cty1, qh%ctz1 /)
    fp(iPB,iCX4:iCZ4) = (/ qh%ctx2, qh%cty2, qh%ctz2 /)

    ! Poloidal front face: row 6 of the array
    fp(iPF,iNX:iNZ)   = (/ qh%npfx, qh%npfy, qh%npfz /)
    fp(iPF,iCX1:iCZ1) = (/ qh%cbx4, qh%cby4, qh%cbz4 /)
    fp(iPF,iCX2:iCZ2) = (/ qh%cbx3, qh%cby3, qh%cbz3 /)
    fp(iPF,iCX3:iCZ3) = (/ qh%ctx3, qh%cty3, qh%ctz3 /)
    fp(iPF,iCX4:iCZ4) = (/ qh%ctx4, qh%cty4, qh%ctz4 /)

end subroutine qhex_face_parameters

!-------------------------------------------------------------------------------
! check_qhex(qh)
! 
! Checks a quadrilateral hexahedron for erroneous or otherwise undesirable 
! properties. If one or more such properties are found, the error flag of the 
! input trapezoid is set to true.
!
! List of erroneous properties:
!     - The height of the top face above the base face is negative
!     - Any pair of edges of the top face intersect one another
!     - Any pair of edges of the bottom face intersect one another
!-------------------------------------------------------------------------------
subroutine check_qhex(qh)

    use geometry, only: plane_elev

    implicit none

    type(quad_hexahedron) :: qh
    real(dp) :: height

    qh%err = .false.

    !---------------------------------------------------------------------------
    ! Ensure that the top plane's height along the stack normal is positive
    !---------------------------------------------------------------------------

    height = (qh%otx-qh%obx)*qh%snx + (qh%oty-qh%oby)*qh%sny + &
             (qh%otz-qh%obz)*qh%snz
    if (height < 0.0) qh%err = .true.

    !---------------------------------------------------------------------------
    ! Ensure that lateral faces are not crossing and have not switched places
    ! (note: this test may make the next two tests redundant)
    !---------------------------------------------------------------------------

    ! Toroidal front and back planes (elev. of toroidal back points above 
    ! toroidal front plane must be negative, as the normal vectors point 
    ! outward)
    if (plane_elev(qh%ctx1, qh%cty1, qh%ctz1, &
                   qh%otfx, qh%otfy, qh%otfz, &
                   qh%ntfx, qh%ntfy, qh%ntfz   ) > 0) qh%err = .true.
    if (plane_elev(qh%ctx4, qh%cty4, qh%ctz4, &
                   qh%otfx, qh%otfy, qh%otfz, &
                   qh%ntfx, qh%ntfy, qh%ntfz   ) > 0) qh%err = .true.
    if (plane_elev(qh%cbx1, qh%cby1, qh%cbz1, &
                   qh%otfx, qh%otfy, qh%otfz, &
                   qh%ntfx, qh%ntfy, qh%ntfz   ) > 0) qh%err = .true.
    if (plane_elev(qh%cbx4, qh%cby4, qh%cbz4, &
                   qh%otfx, qh%otfy, qh%otfz, &
                   qh%ntfx, qh%ntfy, qh%ntfz   ) > 0) qh%err = .true.

    ! Poloidal front and back planes
    if (plane_elev(qh%ctx1, qh%cty1, qh%ctz1, &
                   qh%opfx, qh%opfy, qh%opfz, &
                   qh%npfx, qh%npfy, qh%npfz   ) > 0) qh%err = .true.
    if (plane_elev(qh%ctx2, qh%cty2, qh%ctz2, &
                   qh%opfx, qh%opfy, qh%opfz, &
                   qh%npfx, qh%npfy, qh%npfz   ) > 0) qh%err = .true.
    if (plane_elev(qh%cbx1, qh%cby1, qh%cbz1, &
                   qh%opfx, qh%opfy, qh%opfz, &
                   qh%npfx, qh%npfy, qh%npfz   ) > 0) qh%err = .true.
    if (plane_elev(qh%cbx2, qh%cby2, qh%cbz2, &
                   qh%opfx, qh%opfy, qh%opfz, &
                   qh%npfx, qh%npfy, qh%npfz   ) > 0) qh%err = .true.

end subroutine check_qhex

end module qhex_properties


!-------------------------------------------------------------------------------
! qhex_adjust.f90
!
! Contains subroutines that adjust the geometry of a quadrilaterally-faced
! hexahedron and ensure consistency of all impacted parameters.
!-------------------------------------------------------------------------------
module qhex_adjust

use magpie_globals, only: dp

implicit none

integer, parameter :: stell_transf_mode_len = 9

contains

!-------------------------------------------------------------------------------
! check_and_adjust(qh1, qh2, interval)
!
! Checks for overlaps between two quad hexahedra and retracts overlapping faces
! until no overlaps exist or until one or both become erroneous
!
! Updating parameters:
!     type(quad_hexahedron) :: qh1, qh2   -> the two qhex objects to be checked
!
! Input parameter:
!     real(dp) :: interval  -> Distance by which intersecting faces are to be
!                              retracted on each iteration when an overlap is
!                              found
!-------------------------------------------------------------------------------
subroutine check_and_adjust(qh1, qh2, interval)

    use qhex_properties, only: quad_hexahedron, nHexFaces
    use qhex_overlap,    only: qhex_overlapping

    implicit none

    type(quad_hexahedron), intent(INOUT) :: qh1, qh2
    real(dp), intent(IN) :: interval
    logical, dimension(nHexFaces) :: ovl1, ovl2

    do while (qhex_overlapping(qh1, qh2, ovl1, ovl2))
        if (qh1%err .or. qh2%err) exit
        call adjust_lateral_faces(qh1, qh2, ovl1, ovl2, interval)
    end do

end subroutine check_and_adjust

!-------------------------------------------------------------------------------
! check_and_adjust_pair(pair1, pair2, interval)
!
! Wrapper for check_and_adjust that carries out the procedure for two pairs
! of quadrilateral hexahedra. Only does comparisons with upper qhex if 
! applicable. Upper and lower qhex in the same pair are assumed not to overlap.
!
! Updating parameters:
!     type(quad_hexahedron) :: pair1, pair2 -> pairs of qhex to be checked
!
! Input parameter:
!     real(dp) :: interval  -> Distance by which intersecting faces are to be
!                              retracted on each iteration when an overlap is
!                              found
!-------------------------------------------------------------------------------
subroutine check_and_adjust_pair(pair1, pair2, interval)

    use qhex_properties, only: quad_hexahedron, qhex_pair

    implicit none

    type(qhex_pair), intent(INOUT) :: pair1, pair2
    real(dp) :: interval

    call check_and_adjust(pair1%lo, pair2%lo, interval)

    if (pair1%incl_upper) call check_and_adjust(pair1%up, pair2%lo, interval)
    if (pair2%incl_upper) call check_and_adjust(pair1%lo, pair2%up, interval)
    if (pair1%incl_upper .and. pair2%incl_upper) &
        call check_and_adjust(pair1%up, pair2%up, interval)

end subroutine check_and_adjust_pair

!-------------------------------------------------------------------------------
! reset_qhex_height(surf, qh)
!
! Determines the allowable height for a hexahedral box to avoid collision with
! the vessel, based on its edge lines. Will also ensure that a general maximum
! height is not exceeded, if applicable.
!
! This must be called when the face normal vectors are initially declared or
! if a change is made that might increase the size of the qhex.
!
! Input parameters:
!     type(surface) :: surf -> The toroidal surface representing the vacuum 
!                              vessel or other barrier that the hexahedra must
!                              not intersect
!
! In/out parameters:
!     type(quad_hexahedron) :: qh -> The structure representing the qhex whose
!                                    height is to be adjusted.
!-------------------------------------------------------------------------------
subroutine reset_qhex_height(surf, qh, max_height)

    use magpie_globals,  only: gap_lim
    use geometry,        only: plane_elev, unit_vector
    use surface_calc,    only: surface
    use surface_solve,   only: surf_dist_3d
    use qhex_properties, only: quad_hexahedron, qhex_corners

    implicit none

    type(surface), intent(IN) :: surf
    type(quad_hexahedron), intent(INOUT) :: qh
    real(dp), intent(IN) :: max_height
    real(dp) :: ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3, ax4, ay4, az4
    real(dp) :: length_1, length_2, length_3, length_4
    real(dp) :: l0=0.0, theta0, phi0, theta, phi, chi2
    real(dp) :: xv1, yv1, zv1, xv2, yv2, zv2, xv3, yv3, zv3, xv4, yv4, zv4
    real(dp) :: h1, h2, h3, h4, h5, elev_top

    ! Directional unit vector for edge 1
    call unit_vector(qh%cbx1 - qh%ctx1, qh%cby1 - qh%cty1, qh%cbz1 - qh%ctz1, &
                     ax1, ay1, az1)
    
    ! Directional unit vector for edge 2
    call unit_vector(qh%cbx2 - qh%ctx2, qh%cby2 - qh%cty2, qh%cbz2 - qh%ctz2, &
                     ax2, ay2, az2)
    
    ! Directional unit vector for edge 3
    call unit_vector(qh%cbx3 - qh%ctx3, qh%cby3 - qh%cty3, qh%cbz3 - qh%ctz3, &
                     ax3, ay3, az3)
    
    ! Directional unit vector for edge 4
    call unit_vector(qh%cbx4 - qh%ctx4, qh%cby4 - qh%cty4, qh%cbz4 - qh%ctz4, &
                     ax4, ay4, az4)
    
    ! Determine the intersection points on the vessel for each edge
    call surf_dist_3d(surf, qh%ctx1, qh%cty1, qh%ctz1, ax1, ay1, az1, &
                      l0, qh%ves_theta, qh%ves_phi,                   &
                      length_1, theta, phi, xv1, yv1, zv1, chi2        )
    call surf_dist_3d(surf, qh%ctx2, qh%cty2, qh%ctz2, ax2, ay2, az2, &
                      l0, qh%ves_theta, qh%ves_phi,                   &
                      length_2, theta, phi, xv2, yv2, zv2, chi2        )
    call surf_dist_3d(surf, qh%ctx3, qh%cty3, qh%ctz3, ax3, ay3, az3, &
                      l0, qh%ves_theta, qh%ves_phi,                   &
                      length_3, theta, phi, xv3, yv3, zv3, chi2        )
    call surf_dist_3d(surf, qh%ctx4, qh%cty4, qh%ctz4, ax4, ay4, az4, &
                      l0, qh%ves_theta, qh%ves_phi,                   &
                      length_4, theta, phi, xv4, yv4, zv4, chi2        )

    ! Evaluate the elevations of each vessel point above the vessel
    h1 = plane_elev(xv1, yv1, zv1, qh%vx, qh%vy, qh%vz,  &
                    -qh%snx, -qh%sny, -qh%snz)
    h2 = plane_elev(xv2, yv2, zv2, qh%vx, qh%vy, qh%vz,  &
                    -qh%snx, -qh%sny, -qh%snz)
    h3 = plane_elev(xv3, yv3, zv3, qh%vx, qh%vy, qh%vz,  &
                    -qh%snx, -qh%sny, -qh%snz)
    h4 = plane_elev(xv4, yv4, zv4, qh%vx, qh%vy, qh%vz,  &
                    -qh%snx, -qh%sny, -qh%snz)
    h5 = 0.0   ! The elevation of the reference point on the vessel
    elev_top = max(h1, h2, h3, h4, h5) + gap_lim

    ! Reset the reference point for the top of the qhex
    qh%otx = qh%vx - qh%snx * elev_top
    qh%oty = qh%vy - qh%sny * elev_top
    qh%otz = qh%vz - qh%snz * elev_top

    ! Find the height of the qhex
    qh%height = plane_elev(qh%otx, qh%oty, qh%otz, qh%obx, qh%oby, qh%obz, &
                           qh%snx, qh%sny, qh%snz)
                           
    ! Recalculate the top plane's reference point if height exceeds limit
    if (max_height > 0.0 .and. qh%height > max_height) then
        qh%otx = qh%obx + qh%snx * max_height
        qh%oty = qh%oby + qh%sny * max_height
        qh%otz = qh%obz + qh%snz * max_height
        qh%height = max_height
    end if

    ! Re-calculate the corners according to the adjusted top plane height
    call qhex_corners(qh)

end subroutine reset_qhex_height

!-------------------------------------------------------------------------------
! adjust_lateral_faces(qhex1, qhex2, ovl1, ovl2)
! 
! Wrapper for adjust_face: adjusts the faces of a pair of hexahedra as 
! indicated by corresponding logical input arrays.
!
! Input parameters:
!     type(quad_hexahedron) :: qhex1, qhex2  -> pair of qhex to be adjusted
!     logical, dimension(nHexFaces) :: ovl1, ovl2 
!                                      -> logical arrays whose values in the
!                                         positions corresponding to lateral
!                                         face IDs determine whether the 
!                                         respective lateral face should be 
!                                         adjusted
!-------------------------------------------------------------------------------
subroutine adjust_lateral_faces(qhex1, qhex2, ovl1, ovl2, interval)

    use qhex_properties, only: quad_hexahedron, nHexFaces, iTB, iTF, iPB, iPF

    implicit none

    type(quad_hexahedron), intent(INOUT) :: qhex1, qhex2
    logical, dimension(nHexFaces), intent(IN) :: ovl1, ovl2
    real(dp), intent(IN) :: interval

    ! Adjust the faces of hexahedron 1
    if (ovl1(iTB)) call adjust_face(qhex1, iTB, -abs(interval))
    if (ovl1(iTF)) call adjust_face(qhex1, iTF, -abs(interval))
    if (ovl1(iPB)) call adjust_face(qhex1, iPB, -abs(interval))
    if (ovl1(iPF)) call adjust_face(qhex1, iPF, -abs(interval))

    ! Adjust the faces of hexahedron 2
    if (ovl2(iTB)) call adjust_face(qhex2, iTB, -abs(interval))
    if (ovl2(iTF)) call adjust_face(qhex2, iTF, -abs(interval))
    if (ovl2(iPB)) call adjust_face(qhex2, iPB, -abs(interval))
    if (ovl2(iPF)) call adjust_face(qhex2, iPF, -abs(interval))

end subroutine adjust_lateral_faces

!-------------------------------------------------------------------------------
! adjust_face(qhex, face, dist)
!
! Modifies the geometry of a quad-hexahedron by moving one of the lateral faces
! (toroidal front/back or poloidal front/back) by a specified distance along
! the respective face's normal vector. The geometry of the hexahedron is 
! then updated to reflect this displacement.
!
! Input parameters:
!     type(quad_hexahedron) :: qhex  -> Trapezoid to be modified
!     integer :: face          -> Index of face to modify (iTB, iTF, iPB, iPF)
!     REAL(8) :: dist          -> (signed) displacement dist. along face normal
!------------------------------------------------------------------------------
subroutine adjust_face(qhex, face, dist)

    use qhex_properties, only: quad_hexahedron, qhex_corners, qhex_volume_ctr, &
                               check_qhex, iTOP, iBASE, iTB, iTF, iPB, iPF

    implicit none

    type(quad_hexahedron), intent(INOUT) :: qhex
    integer, intent(IN) :: face
    REAL(8), intent(IN) :: dist

    ! Change reference point of specified face to attain desired displacement
    select case (face)
        case (iTB)
            qhex%otbx = qhex%otbx + dist*qhex%ntbx
            qhex%otby = qhex%otby + dist*qhex%ntby
            qhex%otbz = qhex%otbz + dist*qhex%ntbz
        case (iTF)
            qhex%otfx = qhex%otfx + dist*qhex%ntfx
            qhex%otfy = qhex%otfy + dist*qhex%ntfy
            qhex%otfz = qhex%otfz + dist*qhex%ntfz
        case (iPB)
            qhex%opbx = qhex%opbx + dist*qhex%npbx
            qhex%opby = qhex%opby + dist*qhex%npby
            qhex%opbz = qhex%opbz + dist*qhex%npbz
        case (iPF)
            qhex%opfx = qhex%opfx + dist*qhex%npfx
            qhex%opfy = qhex%opfy + dist*qhex%npfy
            qhex%opfz = qhex%opfz + dist*qhex%npfz
        case (iBASE)
            qhex%obx  = qhex%obx  - dist*qhex%snx
            qhex%oby  = qhex%oby  - dist*qhex%sny
            qhex%obz  = qhex%obz  - dist*qhex%snz
        case (iTOP)
            qhex%otx  = qhex%otx  + dist*qhex%snx
            qhex%oty  = qhex%oty  + dist*qhex%sny
            qhex%otz  = qhex%otz  + dist*qhex%snz
        case default
            stop 'adjust_face: unrecognized face ID'
    end select
        
    ! Update the hexahedron geometry based on this modification
    call qhex_corners(qhex)
    call qhex_volume_ctr(qhex)
    call check_qhex(qhex)

end subroutine adjust_face

!-------------------------------------------------------------------------------
! stell_qhex_transform(mode, phi, qhIn, qhOut)
!
! Transforms a quadrilaterally-faced hexahedron according to the mode selected:
!     reflect:   Reflects a quadrilaterally-faced hexahedron in a given poloidal
!                plane presumed to form the boundary between two stellarator-
!                symmetryc half-periods. The output qhex will have the 
!                equivalent position and orientation associated with the 
!                adjacent half-period. Size and shape are preserved.
!     translate: Translates a quadrilaterally-faced hexahedron in the toroidal 
!                direction; i.e., along a circle of constant radius about the 
!                z-axis. The radial and toroidal orientation of the qhex is 
!                is preserved. Shape and size are preserved as well.
!
! Input parameters:
!     character(len=*) :: mode     -> 'reflect' or 'translate' (see above)
!     real(dp) :: phi              -> 'reflect' mode: toroidal angle (radians)
!                                         of the symmetry plane
!                                     'translate' mode: toroidal interval
!                                         (radians) along which to translate
!     type(quad_hexahedron) qhIn -> qhex to be reflected
!
! Output parameters:
!     type(quad_hexahedron) qhOut -> transformed qhex
!-------------------------------------------------------------------------------
subroutine stell_qhex_transform(mode, phi, qhIn, qhOut)

    use geometry,        only: stell_point_transform, stell_vector_transform
    use qhex_properties, only: quad_hexahedron, qhex_corners

    implicit none

    character(len=*), intent(IN) :: mode
    real(dp), intent(IN) :: phi
    type(quad_hexahedron), intent(IN) :: qhIn
    type(quad_hexahedron), intent(OUT) :: qhOut

    ! Transform the bounding planes
    call stell_vector_transform(mode, phi, qhIn%otfx,  qhIn%otfy,  qhIn%otfz, &
                                           qhIn%ntfx,  qhIn%ntfy,  qhIn%ntfz, &
                                          qhOut%otfx, qhOut%otfy, qhOut%otfz, &
                                          qhOut%ntfx, qhOut%ntfy, qhOut%ntfz  )
    call stell_vector_transform(mode, phi,  qhIn%otbx,  qhIn%otby,  qhIn%otbz, &
                                            qhIn%ntbx,  qhIn%ntby,  qhIn%ntbz, &
                                           qhOut%otbx, qhOut%otby, qhOut%otbz, &
                                           qhOut%ntbx, qhOut%ntby, qhOut%ntbz  )
    call stell_vector_transform(mode, phi,  qhIn%opfx,  qhIn%opfy,  qhIn%opfz, &
                                            qhIn%npfx,  qhIn%npfy,  qhIn%npfz, &
                                           qhOut%opfx, qhOut%opfy, qhOut%opfz, &
                                           qhOut%npfx, qhOut%npfy, qhOut%npfz  )
    call stell_vector_transform(mode, phi,  qhIn%opbx,  qhIn%opby,  qhIn%opbz, &
                                            qhIn%npbx,  qhIn%npby,  qhIn%npbz, &
                                           qhOut%opbx, qhOut%opby, qhOut%opbz, &
                                           qhOut%npbx, qhOut%npby, qhOut%npbz  )
    call stell_vector_transform(mode, phi,  qhIn%otx,  qhIn%oty,  qhIn%otz, &
                                            qhIn%snx,  qhIn%sny,  qhIn%snz, &
                                           qhOut%otx, qhOut%oty, qhOut%otz, &
                                           qhOut%snx, qhOut%sny, qhOut%snz   )
    call stell_vector_transform(mode, phi,  qhIn%obx,  qhIn%oby,  qhIn%obz, &
                                            qhIn%snx,  qhIn%sny,  qhIn%snz, &
                                           qhOut%obx, qhOut%oby, qhOut%obz, &
                                           qhOut%snx, qhOut%sny, qhOut%snz   )

    ! Calculate the corner points based on the bounding planes
    call qhex_corners(qhOut)

    ! Transform points on the vessel and the reference surface
    call stell_point_transform(mode, phi,  qhIn%vx,  qhIn%vy,  qhIn%vz, &
                                        qhOut%vx, qhOut%vy, qhOut%vz   )
    call stell_point_transform(mode, phi,  qhIn%rx,  qhIn%ry,  qhIn%rz, &
                                        qhOut%rx, qhOut%ry, qhOut%rz   )
    call stell_point_transform(mode, phi,  qhIn%ox,  qhIn%oy,  qhIn%oz, &
                                        qhOut%ox, qhOut%oy, qhOut%oz   )


    ! Transform the reference angles (note: mode is checked for validity in 
    ! calls to stell_vector_transform and stell_point_reflect)
    if (trim(mode) == 'reflect') then
        qhOut%ref_phi = phi - (phi + qhIn%ref_phi)
        qhOut%ves_phi = phi - (phi + qhIn%ves_phi)
        qhOut%ref_theta = -qhIn%ref_theta
        qhOut%ves_theta = -qhIn%ves_theta   
    else
        qhOut%ref_phi = qhIn%ref_phi + phi
        qhOut%ves_phi = qhIn%ves_phi + phi
        qhOut%ref_theta = qhIn%ref_theta
        qhOut%ves_theta = qhIn%ves_theta
    end if

    ! Values that remain the same
    qhOut%vol    = qhIn%vol
    qhOut%height = qhIn%height
    qhOut%err    = qhIn%err

end subroutine stell_qhex_transform

!-------------------------------------------------------------------------------
! stell_qhPair_transform(mode, phi, qhIn, qhOut)
!
! Wrapper routine that calls stell_qhex_transform for each of the qhexs within
! a qhex_pair structure. See documentation for stell_qhex_transform for more
! information. The 'incl_upper' field is preserved.
!
! Input parameters:
!     character(len=*) :: mode     -> 'reflect' or 'translate' 
!     real(dp) :: phi              -> 'reflect' mode: toroidal angle (radians)
!                                         of the symmetry plane
!                                     'translate' mode: toroidal interval
!                                         (radians) along which to translate
!     type(quad_hexahedron) qhPairIn -> qhex pair to be transformed
!
! Output parameters:
!     type(quad_hexahedron) qhPairOut -> transformed qhex pair
!-------------------------------------------------------------------------------
subroutine stell_qhPair_transform(mode, phi, qhPairIn, qhPairOut)

    use qhex_properties, only: qhex_pair

    implicit none

    character(len=*), intent(IN) :: mode
    real(dp), intent(IN) :: phi
    type(qhex_pair), intent(IN) :: qhPairIn
    type(qhex_pair), intent(OUT) :: qhPairOut

    call stell_qhex_transform(mode, phi, qhPairIn%lo, qhPairOut%lo)
    call stell_qhex_transform(mode, phi, qhPairIn%up, qhPairOut%up)
    qhPairOut%incl_upper = qhPairIn%incl_upper

end subroutine stell_qhPair_transform

end module qhex_adjust


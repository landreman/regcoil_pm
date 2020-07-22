!-------------------------------------------------------------------------------
! qhex_overlap.f90
!
! Module with subroutines for detecting overlapping quadrateral hexahedra.
!
! Author: K. C. Hammond
! Contact: khammond@pppl.gov
! Updated: 2020-01-02
!-------------------------------------------------------------------------------
module qhex_overlap

use magpie_globals,  only: dp

implicit none

contains

!-------------------------------------------------------------------------------
! qhex_overlapping(qh1, qh2, ovl1, ovl2)
!
! Checks whether two quadrilaterally-faced hexahedra overlap one another by
! iterating through each distinct pair of faces (one from each qhex) and 
! calling the subroutine qhex_face_intersect to detect face intersection.
!
! Input parameters:
!     type(quad_hexahedron) :: qh1, qh2 -> qhex structures for the 2 hexahedra
!                                          to be checked
!
! Output parameters:
!     logical, dimension(6) :: ovl1, ovl2 -> logical arrays indicating which
!                                            faces on each qhex intersect at
!                                            least one face on the other qhex
!-------------------------------------------------------------------------------
logical function qhex_overlapping(qh1, qh2, ovl1, ovl2)

    use qhex_properties, only: quad_hexahedron, nHexFaces, nQhexFaceProps, &
                               qhex_face_parameters 

    implicit none

    type(quad_hexahedron), intent(IN) :: qh1, qh2
    logical, dimension(nHexFaces), intent(OUT) :: ovl1, ovl2
    integer :: i, j
    real(dp), dimension(nHexFaces, nQhexFaceProps) :: fp1, fp2

    ! Collect the face parameters of each qhex
    call qhex_face_parameters(qh1, fp1)
    call qhex_face_parameters(qh2, fp2)

    ! Initialize overlap arrays
    ovl1 = .false.
    ovl2 = .false.

    ! Check all pairs of faces, returning true if at least one intersecting
    ! pair is identified
    qhex_overlapping = .false.
    do i = 1, nHexFaces
        do j = 1, nHexFaces
            if (qhex_face_intersect( &
                fp1(i,1),  fp1(i,2),  fp1(i,3),  fp1(i,4),  fp1(i,5),  &
                fp1(i,6),  fp1(i,7),  fp1(i,8),  fp1(i,9),  fp1(i,10), &
                fp1(i,11), fp1(i,12), fp1(i,13), fp1(i,14), fp1(i,15), &
                fp2(j,1),  fp2(j,2),  fp2(j,3),  fp2(j,4),  fp2(j,5),  &
                fp2(j,6),  fp2(j,7),  fp2(j,8),  fp2(j,9),  fp2(j,10), &
                fp2(j,11), fp2(j,12), fp2(j,13), fp2(j,14), fp2(j,15)   )) then
                    qhex_overlapping = .true.
                    ovl1(i) = .true.
                    ovl2(j) = .true.
            end if
        end do
    end do

end function qhex_overlapping

!-------------------------------------------------------------------------------
! qhex_face_intersect(f1nx, f1ny, f1nz, f1x1, f1y1, f1z1, f1x2, f1y2, f1z2, 
!                     f1x3, f1y3, f1z3, f1x4, f1y4, f1z4, 
!                     f2nx, f2ny, f2nz, f2x1, f2y1, f2z1, f2x2, f2y2, f2z2,
!                     f2x3, f2y3, f2z3, f2x4, f2y4, f2z4)
!
! Determines whether two faces of a quadrilateral hexahedron intersect/collide 
! with one another.
!
! Input parameters:
!     real(dp) :: f1nx, f1ny, f1nz -> x, y, z components, unit normal, face 1
!     real(dp) :: f1x1, f1y1, f1z1 -> x, y, z coords, vertex 1, face 1
!     real(dp) :: f1x2, f1y2, f1z2 -> x, y, z coords, vertex 2, face 1
!     real(dp) :: f1x3, f1y3, f1z3 -> x, y, z coords, vertex 3, face 1
!     real(dp) :: f1x4, f1y4, f1z4 -> x, y, z coords, vertex 4, face 1
!     real(dp) :: f2nx, f2ny, f2nz -> x, y, z components, unit normal, face 2
!     real(dp) :: f2x1, f2y1, f2z1 -> x, y, z coords, vertex 1, face 2
!     real(dp) :: f2x2, f2y2, f2z2 -> x, y, z coords, vertex 2, face 2
!     real(dp) :: f2x3, f2y3, f2z3 -> x, y, z coords, vertex 3, face 2
!     real(dp) :: f2x4, f2y4, f2z4 -> x, y, z coords, vertex 4, face 2
!-------------------------------------------------------------------------------
logical function qhex_face_intersect( &
                     f1nx, f1ny, f1nz, f1x1, f1y1, f1z1, f1x2, f1y2, f1z2, &
                                       f1x3, f1y3, f1z3, f1x4, f1y4, f1z4, &
                     f2nx, f2ny, f2nz, f2x1, f2y1, f2z1, f2x2, f2y2, f2z2, &
                                       f2x3, f2y3, f2z3, f2x4, f2y4, f2z4)

    use geometry, only: cross_prod, two_plane_intersect, segment_crossing_2d

    implicit none

    real(dp), intent(IN) :: f1nx, f1ny, f1nz, f1x1, f1y1, f1z1, f1x2, f1y2, f1z2
    real(dp), intent(IN) :: f1x3, f1y3, f1z3, f1x4, f1y4, f1z4
    real(dp), intent(IN) :: f2nx, f2ny, f2nz, f2x1, f2y1, f2z1, f2x2, f2y2, f2z2
    real(dp), intent(IN) :: f2x3, f2y3, f2z3, f2x4, f2y4, f2z4
    integer :: i, next, ind1, ind2
    real(dp) :: oxi, oyi, ozi, axi, ayi, azi, bx1, by1, bz1, bx2, by2, bz2
    real(dp) :: ais, bis, aie, bie
    real(dp), dimension(4) :: a1, b1, a2, b2
    real(dp), dimension(4) :: s1, s2, s1i, s2i
    logical, dimension(4) :: crossed_qhex1, crossed_qhex2
    integer :: nCrossed_qhex1, nCrossed_qhex2
    real(dp), allocatable :: cross_s1i(:), cross_s2i(:)
    logical :: ignore

    ! Determine the line of intersection between the planes of each face
    call two_plane_intersect(f1x1, f1y1, f1z1, f1nx, f1ny, f1nz, &
                             f2x1, f2y1, f2z1, f2nx, f2ny, f2nz, &
                              oxi,  oyi,  ozi,  axi,  ayi,  azi   )

    ! Determine 2d coordinates for vertices of first face, with the axis of
    ! the intersect line corresponding to the "a" coordinate and the reference
    ! point of the intersect line corresponding to the origin
    call cross_prod(f1nx, f1ny, f1nz, axi, ayi, azi, bx1, by1, bz1)
    a1(1) = (f1x1-oxi)*axi + (f1y1-oyi)*ayi + (f1z1-ozi)*azi
    b1(1) = (f1x1-oxi)*bx1 + (f1y1-oyi)*by1 + (f1z1-ozi)*bz1
    a1(2) = (f1x2-oxi)*axi + (f1y2-oyi)*ayi + (f1z2-ozi)*azi
    b1(2) = (f1x2-oxi)*bx1 + (f1y2-oyi)*by1 + (f1z2-ozi)*bz1
    a1(3) = (f1x3-oxi)*axi + (f1y3-oyi)*ayi + (f1z3-ozi)*azi
    b1(3) = (f1x3-oxi)*bx1 + (f1y3-oyi)*by1 + (f1z3-ozi)*bz1
    a1(4) = (f1x4-oxi)*axi + (f1y4-oyi)*ayi + (f1z4-ozi)*azi
    b1(4) = (f1x4-oxi)*bx1 + (f1y4-oyi)*by1 + (f1z4-ozi)*bz1

    ! Determine 2d coordinates for vertices of the second face, with the axis of
    ! the intersect line corresponding to the "a" coordinate and the reference
    ! point of the intersect line corresponding to the origin
    call cross_prod(f2nx, f2ny, f2nz, axi, ayi, azi, bx2, by2, bz2)
    a2(1) = (f2x1-oxi)*axi + (f2y1-oyi)*ayi + (f2z1-ozi)*azi
    b2(1) = (f2x1-oxi)*bx2 + (f2y1-oyi)*by2 + (f2z1-ozi)*bz2
    a2(2) = (f2x2-oxi)*axi + (f2y2-oyi)*ayi + (f2z2-ozi)*azi
    b2(2) = (f2x2-oxi)*bx2 + (f2y2-oyi)*by2 + (f2z2-ozi)*bz2
    a2(3) = (f2x3-oxi)*axi + (f2y3-oyi)*ayi + (f2z3-ozi)*azi
    b2(3) = (f2x3-oxi)*bx2 + (f2y3-oyi)*by2 + (f2z3-ozi)*bz2
    a2(4) = (f2x4-oxi)*axi + (f2y4-oyi)*ayi + (f2z4-ozi)*azi
    b2(4) = (f2x4-oxi)*bx2 + (f2y4-oyi)*by2 + (f2z4-ozi)*bz2

    ! Start (s) and end (e) points for a segment of unit (arbitrary) length on 
    ! the intersection axis, which are the same in the 2d coordinate systems for
    ! both face 1 and face 2
    ais = 0.0
    bis = 0.0
    aie = 1.0
    bie = 0.0

    ! Determine line intersections between edges of the quadrilateral and the 
    ! plane intersection line, noting whether the edges are crossed by the line
    nCrossed_qhex1 = 0
    nCrossed_qhex2 = 0
    do i = 1, 4
        if (i < 4) then
            next = i + 1
        else
            next = 1
        end if
        call segment_crossing_2d(a1(i), b1(i), ais, bis, &
                                 a1(next), b1(next), aie, bie, &
                                 s1(i), s1i(i), crossed_qhex1(i), ignore)
        call segment_crossing_2d(a2(i), b2(i), ais, bis, &
                                 a2(next), b2(next), aie, bie, &
                                 s2(i), s2i(i), crossed_qhex2(i), ignore)
        if (crossed_qhex1(i)) nCrossed_qhex1 = nCrossed_qhex1 + 1
        if (crossed_qhex2(i)) nCrossed_qhex2 = nCrossed_qhex2 + 1
    end do

    ! Populate arrays containing the positions (s{1,2}i) of the segment 
    ! crossing points along the plane intersect line 
    allocate(cross_s1i(nCrossed_qhex1), cross_s2i(nCrossed_qhex2))
    ind1 = 1
    ind2 = 1
    do i = 1, 4
        if (crossed_qhex1(i)) then
            cross_s1i(ind1) = s1i(i)
            ind1 = ind1 + 1
        end if
        if (crossed_qhex2(i)) then
            cross_s2i(ind2) = s2i(i)
            ind2 = ind2 + 1
        end if
    end do

    ! Determine whether the quadrilateral faces intersect one another; 
    ! i.e., whether the positions of edge crossings with the line of plane 
    ! intersection (if any crossings exist) overlap one another
    if (nCrossed_qhex1 == 0 .or. nCrossed_qhex2 == 0) then
        qhex_face_intersect = .false.
    else if (maxval(cross_s1i) < minval(cross_s2i) .or. &
             maxval(cross_s2i) < minval(cross_s1i)       ) then
        qhex_face_intersect = .false.
    else
        qhex_face_intersect = .true.
    end if
    
    deallocate(cross_s1i, cross_s2i)

end function qhex_face_intersect

!-------------------------------------------------------------------------------
! object_overlap_0(qh, overlaps)
!
! Checks whether any of the corners of a qhex intersect with any applicable
! object (port, coil, etc.).
!
! NOTE: for correct functionality, calling programs must call the 
! intersections_to_check subroutine once prior to calls to object_overlap.
!
! NOTE: deprecated in favor of object_overlap (see below).
!
! Input parameters:
!     type(quad_hexahedron) :: qh  -> Qhex whose corners are to be checked
!
! Output parameter:
!     logical :: overlap -> True if any of the corners intersect an applicable
!                           object.
!-------------------------------------------------------------------------------
subroutine object_overlap_0(qh, overlaps)

    use qhex_properties, only: quad_hexahedron, nQhexCorners
    use intersections,   only: check_intersections

    implicit none

    type(quad_hexahedron), intent(IN) :: qh
    logical, intent(OUT) :: overlaps
    logical, dimension(nQhexCorners) :: corner_overlaps
    integer :: i

    overlaps = .false.

    call check_intersections(qh%ctx1, qh%cty1, qh%ctz1, corner_overlaps(1))
    call check_intersections(qh%ctx2, qh%cty2, qh%ctz2, corner_overlaps(2))
    call check_intersections(qh%ctx3, qh%cty3, qh%ctz3, corner_overlaps(3))
    call check_intersections(qh%ctx4, qh%cty4, qh%ctz4, corner_overlaps(4))
    call check_intersections(qh%cbx1, qh%cby1, qh%cbz1, corner_overlaps(5))
    call check_intersections(qh%cbx2, qh%cby2, qh%cbz2, corner_overlaps(6))
    call check_intersections(qh%cbx3, qh%cby3, qh%cbz3, corner_overlaps(7))
    call check_intersections(qh%cbx4, qh%cby4, qh%cbz4, corner_overlaps(8))

    do i = 1, nQhexCorners
        if (corner_overlaps(i)) overlaps = .true.
    end do

end subroutine object_overlap_0

!-------------------------------------------------------------------------------
! object_overlap(qh, maxCheckGap, overlaps)
!
! Checks whether the faces of a qhex intersect with any applicable object
! (port, coil, etc.).
!
! NOTE: for correct functionality, calling procedures must call the 
! intersections_to_check subroutine once prior to calls to object_overlap.
!
! Input parameters:
!     type(quad_hexahedron) :: qh -> Qhex whose faces are to be checked
!     real(dp) :: maxCheckGap     -> Max spacing between adjacent points on 
!                                    each face for overlap checking
!
! Output parameters:
!     logical :: overlaps  -> True if any faces are found to intersect an 
!                             applicable object
!-------------------------------------------------------------------------------
subroutine object_overlap(qh, maxCheckGap, overlaps)

    use qhex_properties, only: quad_hexahedron, nHexFaces, &
                               nQhexFaceProps, qhex_face_parameters 

    implicit none

    type(quad_hexahedron), intent(IN) :: qh
    real(dp), intent(IN) :: maxCheckGap
    logical, intent(OUT) :: overlaps
    real(dp), dimension(nHexFaces,nQhexFaceProps) :: fparams
    integer :: i

    overlaps = .false.

    call qhex_face_parameters(qh, fparams)

    do i = 1, nHexFaces
        overlaps = overlaps_face(fparams(i,:), maxCheckGap)
        if (overlaps) return
    end do

end subroutine object_overlap

!-------------------------------------------------------------------------------
! overlaps_face(fpRow, maxCheckGap)
!
! Checks whether any of a set of regularly-spaced points on the face of a qhex 
! intersects an external object.
!
! Input parameters:
!     real(dp), dimension(nQhexFaceProps) :: fpRow
!                             -> Row of an array of face parameters as generated
!                                by the qhex_face_parameters subroutine 
!                                (qhex_properties module)
!     real(dp) :: maxCheckGap -> Maximum linear distance between adjacent points
!                                in the grid to check for intersections
!
! Returns:
!     .true. if any of the query points intersect an external object;
!     .false. otherwise
!-------------------------------------------------------------------------------
logical function overlaps_face(fpRow, maxCheckGap)

    use geometry,        only: linspace3d
    use qhex_properties, only: iCX1, iCY1, iCZ1, iCX2, iCY2, iCZ2, &
                               iCX3, iCY3, iCZ3, iCX4, iCY4, iCZ4, &
                               nQhexFaceProps
    use intersections,   only: check_intersections

    implicit none

    real(dp), dimension(nQhexFaceProps), intent(IN) :: fpRow
    real(dp), intent(IN) :: maxCheckGap
    real(dp) :: edgeLg12sq, edgeLg23sq, edgeLg34sq, edgeLg41sq
    real(dp) :: evenEdgeLgMaxSq, oddEdgeLgMaxSq
    real(dp), dimension(:), allocatable :: ox0, oy0, oz0, ox1, oy1, oz1
    real(dp), dimension(:), allocatable :: xq, yq, zq
    integer :: i, j, nEvenEdgePts, nOddEdgePts

    ! Initialize the return value
    overlaps_face = .false.

    ! Determine number of cuts on "even" and "odd" edges based on edge length
    edgeLg12sq = (fpRow(iCX2)-fpRow(iCX1))**2 + (fpRow(iCY2)-fpRow(iCY1))**2 + &
                 (fpRow(iCZ2)-fpRow(iCZ1))**2
    edgeLg23sq = (fpRow(iCX3)-fpRow(iCX2))**2 + (fpRow(iCY3)-fpRow(iCY2))**2 + &
                 (fpRow(iCZ3)-fpRow(iCZ2))**2
    edgeLg34sq = (fpRow(iCX4)-fpRow(iCX3))**2 + (fpRow(iCY4)-fpRow(iCY3))**2 + &
                 (fpRow(iCZ4)-fpRow(iCZ3))**2
    edgeLg41sq = (fpRow(iCX1)-fpRow(iCX4))**2 + (fpRow(iCY1)-fpRow(iCY4))**2 + &
                 (fpRow(iCZ1)-fpRow(iCZ4))**2
    oddEdgeLgMaxSq  = max( edgeLg12sq, edgeLg34sq )
    evenEdgeLgMaxSq = max( edgeLg23sq, edgeLg41sq )
    nOddEdgePts  = ceiling( sqrt(oddEdgeLgMaxSq)  / maxCheckGap ) + 1
    nEvenEdgePts = ceiling( sqrt(evenEdgeLgMaxSq) / maxCheckGap ) + 1

    ! Allocate arrays of start (0) and end (1) points on the two "odd" edges
    allocate(ox0(nOddEdgePts),  oy0(nOddEdgePts),  oz0(nOddEdgePts), &
             ox1(nOddEdgePts),  oy1(nOddEdgePts),  oz1(nOddEdgePts), &
             xq(nEvenEdgePts),  yq(nEvenEdgePts),  zq(nEvenEdgePts)   )
    
    ! Distribute start and end points evenly along each odd edge
    call linspace3d(nOddEdgePts,  fpRow(iCX1), fpRow(iCY1), fpRow(iCZ1), &
                    fpRow(iCX2), fpRow(iCY2), fpRow(iCZ2), ox0, oy0, oz0)
    call linspace3d(nOddEdgePts,  fpRow(iCX4), fpRow(iCY4), fpRow(iCZ4), &
                    fpRow(iCX3), fpRow(iCY3), fpRow(iCZ3), ox1, oy1, oz1)

    do i = 1, nOddEdgePts

        ! Space query points evenly along lines connecting odd start/end points
        call linspace3d(nEvenEdgePts, ox0(i), oy0(i), oz0(i), &
                        ox1(i), oy1(i), oz1(i), xq, yq, zq)


        ! Check each query point for intersections; exit immediately if found
        do j = 1, nEvenEdgePts
            call check_intersections(xq(j), yq(j), zq(j), overlaps_face)
            if (overlaps_face) exit 
        end do

        if (overlaps_face) exit

    end do

    deallocate(ox0, oy0, oz0, ox1, oy1, oz1, xq, yq, zq)

end function overlaps_face

end module qhex_overlap


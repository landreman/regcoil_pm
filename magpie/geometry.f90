!-------------------------------------------------------------------------------
! geometry.f90
!
! Module containing procedures for geometric calculations used in MAGPIE
!
! Author: K. C. Hammond
! Contact: khammond@pppl.gov
! Updated: 2020-01-01
!-------------------------------------------------------------------------------

module geometry

use magpie_globals, only: dp

implicit none

contains

!-------------------------------------------------------------------------------
! unit_vector(vx, vy, vz, ux, uy, uz)
!
! Scales an input vector to a unit-length vector of the same direction.
!
! Input parameters:
!     real(dp) :: vx, vy, vz   -> x, y, z components of the input vector
!
! Output parameters:
!     real(dp) :: ux, uy, uz   -> x, y, z components of the rescaled unit vector
!-------------------------------------------------------------------------------
subroutine unit_vector(vx, vy, vz, ux, uy, uz)

    implicit none

    real(dp), intent(IN)  :: vx, vy, vz
    real(dp), intent(OUT) :: ux, uy, uz
    real(dp) :: norm

    norm = sqrt(vx**2 + vy**2 + vz**2)
    ux = vx / norm
    uy = vy / norm
    uz = vz / norm

end subroutine unit_vector

!-------------------------------------------------------------------------------
! cross_prod(ax, ay, az, bx, by, bz, cx, cy, cz)
! 
! Computes the cross product of two Cartesian 3-dimensional vectors
!
! Input parameters
!     real(dp) :: ax, ay, az  -> x, y, and z components of the first vector
!     real(dp) :: bx, by, bz  -> x, y, and z components of the second vector
!
! Output parameters
!     real(dp) :: cx, cy, cz  -> x, y, and z components of a cross b
!-------------------------------------------------------------------------------
subroutine cross_prod(ax, ay, az, bx, by, bz, cx, cy, cz)

    implicit none

    real(dp), intent(IN) :: ax, ay, az, bx, by, bz
    real(dp), intent(OUT) :: cx, cy, cz

    cx =  ay*bz - az*by
    cy = -ax*bz + az*bx
    cz =  ax*by - ay*bx

end subroutine cross_prod

!-------------------------------------------------------------------------------
! segment_crossing_2d(x1s, y1s, x2s, y2s, x1e, y1e, x2e, y2e, 
!                     s1, s2, seg1_crossed, seg2_crossed)
!
! Determines the crossing points of lines defined by two segments
! in the x-y plane.
!
! Input parameters:
!     real(dp) :: x1s, y1s -> x, y coordinates, start point, segment 1
!     real(dp) :: x2s, y2s -> x, y coordinates, start point, segment 2
!     real(dp) :: x1e, y1e -> x, y coordinates, end point, segment 1
!     real(dp) :: x2e, y2e -> x, y coordinates, end point, segment 2
!
! Output parameters:
!     real(dp) :: s1, s2   -> Distance from the start point of line 1, 2 to the
!                            intersection, in the direction toward the end point
!     logical :: seg1_crossed  -> True if the first segment is crossed by the
!                                 line parametrized by the second segment
!     logical :: seg2_crossed  -> True if the second segment is crossed by the
!                                 line parametrized by the first segment
!-------------------------------------------------------------------------------
subroutine segment_crossing_2d(x1s, y1s, x2s, y2s, x1e, y1e, x2e, y2e, &
                               s1, s2, seg1_crossed, seg2_crossed)

    use algebra, only: inverse_2x2

    implicit none

    real(dp), intent(IN)  :: x1s, y1s, x2s, y2s, x1e, y1e, x2e, y2e
    real(dp), intent(OUT) :: s1, s2
    logical, intent(OUT) :: seg1_crossed, seg2_crossed
    real(dp) :: ax1, ay1, ax2, ay2, length1, length2
    real(dp) :: inv11, inv12, inv21, inv22

    ! Lengths of each segment
    length1 = sqrt( (x1e-x1s)**2 + (y1e-y1s)**2 )
    length2 = sqrt( (x2e-x2s)**2 + (y2e-y2s)**2 )

    ! Unit axis vectors for each segment
    ax1 = (x1e-x1s) / length1
    ay1 = (y1e-y1s) / length1
    ax2 = (x2e-x2s) / length2
    ay2 = (y2e-y2s) / length2

    ! Inverse of the matrix defined by [ ax1  -ax2;  ay1  -ay2  ]
    call inverse_2x2(ax1, -ax2, ay1, -ay2, inv11, inv12, inv21, inv22)

    ! Distances of the intersection from start points of each line segment
    s1 = inv11*(x2s-x1s) + inv12*(y2s-y1s)
    s2 = inv21*(x2s-x1s) + inv22*(y2s-y1s)

    ! Determine whether the line intersection point is within both segments
    if (s1 >= 0.0 .and. s1 <= length1) then
        seg1_crossed = .true. 
    else
        seg1_crossed = .false.
    end if
    if (s2 >= 0.0 .and. s2 <= length2) then
        seg2_crossed = .true.
    else
        seg2_crossed = .false.
    end if

end subroutine segment_crossing_2d

!-------------------------------------------------------------------------------
! segments_crossed(x1s, y1s, x2s, y2s, x1e, y1e, x2e, y2e)
!
! Wrapper for segment_crossings_2d that returns a logical if two line segments
! in the x-y plane cross one another.
!
! Input parameters:
!     real(dp) x1s, y1s   -> x, y coordinates, start point, segment 1
!     real(dp) x2s, y2s   -> x, y coordinates, start point, segment 2
!     real(dp) x1e, y1e   -> x, y coordinates, end point, segment 1
!     real(dp) x2e, y2e   -> x, y coordinates, end point, segment 2
!-------------------------------------------------------------------------------
logical function segments_crossed(x1s, y1s, x2s, y2s, x1e, y1e, x2e, y2e)

    implicit none

    real(dp), intent(IN) :: x1s, y1s, x2s, y2s, x1e, y1e, x2e, y2e
    real(dp) :: s1, s2
    logical :: seg1_crossed, seg2_crossed

    ! Find the crossings points for each segment
    call segment_crossing_2d(x1s, y1s, x2s, y2s, x1e, y1e, x2e, y2e, &
                             s1, s2, seg1_crossed, seg2_crossed)

    ! Set return value according to the outcome of the subroutine call
    segments_crossed = (seg1_crossed .and. seg2_crossed)

end function segments_crossed

!-------------------------------------------------------------------------------
! plane_elev(ox, oy, oz, nx, ny, nz)
!
! Calculates the elevation of a point above a plane in 3D space.
!
! Input parameters:
!     real(dp) :: x, y, z    -> x, y, z coords of the query/input point
!     real(dp) :: ox, oy, oz -> x, y, z coords of a reference point on the plane
!     real(dp) :: nx, ny, nz -> x, y, z components of the plane normal vector
!
! Returns:
!     real(dp) :: plane_elev -> elevation of the input point above the plane
!-------------------------------------------------------------------------------
real(dp) function plane_elev(x, y, z, ox, oy, oz, nx, ny, nz)

    implicit none

    real(dp), intent(IN)  :: x, y, z, ox, oy, oz, nx, ny, nz
    
    plane_elev = (x - ox)*nx + (y - oy)*ny + (z - oz)*nz
    
end function plane_elev

!-------------------------------------------------------------------------------
! three_plane_intersect(ox1, oy1, oz1, nx1, ny1, nz1, &
!                       ox2, oy2, oz2, nx2, ny2, nz2, &
!                       ox3, oy3, oz3, nx3, ny3, nz3, &
!                       xi, yi, zi)
!
! Calculates the point of intersection between three planes. The approach is
! to solve a matrix equation of the form Ax = b, where:
!     A is a matrix in which each row is a normal vector
!     x contains the components of the intersection point
!     b vector in which each component is the dot product of a normal vector
!       with its corresponding intersection point.
!
! Equivalently, solves the system of 3 equations of the form 
!     (I - On) dot Nn = 0,
! where I is the intersection point, On is the origin point of plane n, and 
! Nn is the normal vector of plane n.
!
! Note: normal vectors are assumed to be linearly independent from one another;
! i.e., a solution is assumed to exist.
!
! Input parameters:
!     real(dp) :: ox1, oy1, oz1 -> x, y, and z coords, ref point on plane 1
!     real(dp) :: nx1, ny1, nz1 -> x, y, and z components, norm vector to plane 1
!     real(dp) :: ox2, oy2, oz2 -> x, y, and z coords, ref point on plane 2
!     real(dp) :: nx2, ny2, nz2 -> x, y, and z components, norm vector to plane 2
!     real(dp) :: ox3, oy3, oz3 -> x, y, and z coords, ref point on plane 3
!     real(dp) :: nx3, ny3, nz3 -> x, y, and z components, norm vector to plane 3
!
! Output parameters:
!     real(dp) :: xi, yi, zi      -> x, y, and z coords, intersection point
!-------------------------------------------------------------------------------
subroutine three_plane_intersect(ox1, oy1, oz1, nx1, ny1, nz1, &
                                 ox2, oy2, oz2, nx2, ny2, nz2, &
                                 ox3, oy3, oz3, nx3, ny3, nz3, &
                                 xi, yi, zi)


    use algebra, only: inverse_3x3

    implicit none

    real(dp), intent(IN)  :: ox1, oy1, oz1, nx1, ny1, nz1
    real(dp), intent(IN)  :: ox2, oy2, oz2, nx2, ny2, nz2
    real(dp), intent(IN)  :: ox3, oy3, oz3, nx3, ny3, nz3
    real(dp), intent(OUT) :: xi, yi, zi
    real(dp) :: ni11, ni12, ni13, ni21, ni22, ni23, ni31, ni32, ni33
    real(dp) :: b1, b2, b3

    ! The right-hand side of the matrix equation (see doc string)
    b1 = nx1*ox1 + ny1*oy1 + nz1*oz1
    b2 = nx2*ox2 + ny2*oy2 + nz2*oz2
    b3 = nx3*ox3 + ny3*oy3 + nz3*oz3

    ! The inverse of the matrix of normal vectors
    call inverse_3x3(nx1,  ny1,  nz1,  nx2,  ny2,  nz2,  nx3,  ny3,  nz3, &
                     ni11, ni12, ni13, ni21, ni22, ni23, ni31, ni32, ni33)

    ! Solution
    xi = ni11*b1 + ni12*b2 + ni13*b3
    yi = ni21*b1 + ni22*b2 + ni23*b3
    zi = ni31*b1 + ni32*b2 + ni33*b3

end subroutine three_plane_intersect

!-------------------------------------------------------------------------------
! two_plane_intersect(ox1, oy1, oz1, nx1, ny1, nz1, &
!                     ox2, oy2, oz2, nx2, ny2, nz2, &
!                     oxi, oyi, ozi, axi, ayi, azi   )
!
! Calculates the line of intersection between two planes. 
!
! Note: if the two planes are parallel, all outputs will be NaN. 
!
! Input parameters:
!     real(dp) :: ox1, oy1, oz1 -> x, y, and z coords, ref point on plane 1
!     real(dp) :: nx1, ny1, nz1 -> x, y, and z components, norm vector to plane 1
!     real(dp) :: ox2, oy2, oz2 -> x, y, and z coords, ref point on plane 2
!     real(dp) :: nx2, ny2, nz2 -> x, y, and z components, norm vector to plane 2
!
! Output parameters:
!     real(dp) :: oxi, oyi, ozi -> x, y, z coords, ref. point, intersect line
!     real(dp) :: axi, ayi, azi -> x, y, z comps, unit vector, intersect axis
!-------------------------------------------------------------------------------
subroutine two_plane_intersect(ox1, oy1, oz1, nx1, ny1, nz1, &
                               ox2, oy2, oz2, nx2, ny2, nz2, &
                               oxi, oyi, ozi, axi, ayi, azi   )

    implicit none

    real(dp), intent(IN)  :: ox1, oy1, oz1, nx1, ny1, nz1
    real(dp), intent(IN)  :: ox2, oy2, oz2, nx2, ny2, nz2
    real(dp), intent(OUT) :: oxi, oyi, ozi, axi, ayi, azi
    real(dp) :: axi_vec, ayi_vec, azi_vec
    real(dp) :: l, bx1, by1, bz1, diff_dot_n2, b_dot_n2

    ! The direction of the line of intersection
    call cross_prod(nx1, ny1, nz1, nx2, ny2, nz2, axi_vec, ayi_vec, azi_vec)
    call unit_vector(axi_vec, ayi_vec, azi_vec, axi, ayi, azi)

    ! Vector perpendicular to plane 1's normal vector and the intersect axis
    call cross_prod(nx1, ny1, nz1, axi, ayi, azi, bx1, by1, bz1)

    ! Distance l along b vector from plane 1 origin to reach intersection line:
    !     (o1 + l*b - o2).n2 = 0
    !           (o1 - o2).n2 = -l*b.n2
    !                      l = (o2 - o1).n2 / b.n2
    diff_dot_n2 = (ox2 - ox1)*nx2 + (oy2 - oy1)*ny2 + (oz2 - oz1)*nz2
    b_dot_n2 = bx1*nx2 + by1*ny2 + bz1*nz2
    l = diff_dot_n2 / b_dot_n2

    oxi = ox1 + l*bx1
    oyi = oy1 + l*by1
    ozi = oz1 + l*bz1

end subroutine two_plane_intersect

!-------------------------------------------------------------------------------
! in_polygon(np, xp, yp, xq, yq), result(inside)
!
! Determines whether a query point lies within a closed, two-dimensional 
! polygon. Points on bounding segments are included. The input points need
! not be closed (i.e., it is not necessary for the last polygon point to be
! equal to the first), although the function will work either way.
!
! Uses the winding number method described here by D. Sunday:
! http://geomalgorithms.com/a03-_inclusion.html  (accessed 2020-04-09)
! but modified to consider all points and bounding segments as "inside" the
! polygon.
!
! Input parameters:
!     integer :: np       -> Number of points in the polygon
!     real(dp) :: xp, yp  -> x, y coordinates of the points in the polygon
!     real(dp) :: xq, yq  -> x, y coordinates of the query point
! 
! Return value:
!     logical :: inside   -> True if the query point is inside the polygon; 
!                            false if it is outside (or on a bounding segment)
!-------------------------------------------------------------------------------
function in_polygon(np, xp, yp, xq, yq) result(inside)

    implicit none

    integer, intent(IN) :: np
    real(dp), dimension(np), intent(IN) :: xp, yp
    real(dp), intent(IN) :: xq, yq
    logical :: inside
    integer :: windNo, i, next
    real(dp) :: ccw

    ! Set winding number counter to zero
    windNo = 0

    do i = 1, np

        next = i + 1
        if (i==np) next = 1

        ! Consider a horizontal ray extending in the positive x direction from
        ! the query point.

        ! For upward segments that cross the ray, increment the winding number
        if (yp(i) <= yq) then
            if (yp(next) > yq) then
                ccw = counterclockwise(xp(i), yp(i), &
                                       xp(next), yp(next), xq, yq)
                if (ccw > 0) then
                    windNo = windNo + 1
                elseif (ccw == 0) then
                    inside = .true.
                    return
                endif
            endif

        ! For downward segments that cross the ray, decrement the winding number
        else
            if (yp(next) <= yq) then
                ccw = counterclockwise(xp(i), yp(i), &
                                       xp(next), yp(next), xq, yq)
                if (ccw < 0) then
                    windNo = windNo - 1
                elseif (ccw == 0) then
                    inside = .true.
                    return
                endif
            endif
        endif
 
    end do

    if (windNo /= 0) then
        inside = .true.
    else
        inside = .false.
    endif

end function in_polygon

!-------------------------------------------------------------------------------
! counterclockwise(x0, y0, x1, y1, xq, yq), result(ccw)
!
! Determines whether a query point is counterclockwise relative to a ray in
! 2d space defined by a pair of input points.
!
! Input parameters:
!     real(dp) :: x0, y0    -> x, y coordinates of the beginning of the ray
!     real(dp) :: x1, y1    -> x, y coordinates of a second point along the ray
!     real(dp) :: xq, yq    -> x, y coordinates of the query point
!
! Return value:
!     real(dp) :: ccw -> >0 if the query point is counterclockwise from the ray;
!                        =0 if the query point is collinear with the ray;
!                        <0 if the query point is clockwise from the ray
!-------------------------------------------------------------------------------
function counterclockwise(x0, y0, x1, y1, xq, yq) result(ccw)

    implicit none

    real(dp), intent(IN) :: x0, y0, x1, y1, xq, yq
    real(dp) :: ccw

    ccw = (x1-x0)*(yq-y0) - (y1-y0)*(xq-x0)

end function counterclockwise

!-------------------------------------------------------------------------------
! point_segment_lengths(nPoints, s_ideal, r, z, s, meanSqResid)
!
! For a given set of input points in two dimensions, calculates a cumulative
! length array for the segments between successive vertices.
!
! Input parameters:
!     integer :: nPoints     -> Number of points in the poloidal array
!     real(dp) :: s_ideal(:) -> Desired normalized cumulative length btwn points
!     real(dp) :: r(:), z(:) -> r and z coordinates of the input points
!
! Output parameters:
!     real(dp) :: s           -> Actual normalized cumulative length btwn points
!     real(dp) :: meanSqResid -> Mean squared residual, scaled to total length
!-------------------------------------------------------------------------------
subroutine point_segment_lengths(nPoints, s_ideal, r, z, s, meanSqResid)
    
    implicit none

    integer, intent(IN) :: nPoints
    real(dp), dimension(nPoints), intent(IN)  :: s_ideal, r, z
    real(dp), dimension(nPoints), intent(OUT) :: s
    real(dp), intent(OUT) :: meanSqResid
    integer :: i
    real(dp) :: dl, l_total
    real(dp), dimension(nPoints) :: l

    ! Compute the individual (dl) and cumulative (l) segment ("arc") lengths
    l(1) = 0.0
    do i = 1, nPoints-1
        dl = sqrt((r(i+1) - r(i))**2 + (z(i+1) - z(i))**2)
        l(i+1) = l(i) + dl
    end do

    ! Normalize the arc lengths to the interval [0,1]
    l_total = l(nPoints)
    s = l / l_total

    ! Calculate the mean squared residual
    meanSqResid = 0.0
    do i = 2, nPoints-1
        meanSqResid = meanSqResid + (s(i) - s_ideal(i))**2
    end do

    ! Calculate the mean and rescale to the squared total length
    meanSqResid = l_total**2 * meanSqResid/(nPoints - 2)

end subroutine point_segment_lengths

!-------------------------------------------------------------------------------
! linspace3d(n, x0, y0, z0, x1, y1, z1, xArr, yArr, zArr
!
! Generates arrays of Cartesian coordinates of points that are evenly spaced
! along a given line segment.
!
! Behavior is analogous to the linspace() functions in Matlab and Python Numpy.
!
! Input parameters:
!     integer :: n
!     real(dp) :: x0, y0, z0  -> x, y, z coordinates, beginning of the segment
!     real(dp) :: x1, y1, z1  -> x, y, z coordinates, end of the segment
!
! Output parameters:
!     real(dp), dimension(:) :: xArr, yArr, zArr
!                             -> Arrays of the coordinates of the set of points
!                                beginning with x0, y0, and z0, respectively, 
!                                and ending with x1, y1, and z1, respectively.
!-------------------------------------------------------------------------------
subroutine linspace3d(n, x0, y0, z0, x1, y1, z1, xArr, yArr, zArr)

    implicit none

    integer, intent(IN) :: n
    real(dp), intent(IN) :: x0, y0, z0, x1, y1, z1
    real(dp), dimension(:), intent(OUT) :: xArr, yArr, zArr
    real(dp) :: dx, dy, dz
    integer :: i

    dx = (x1 - x0)/(n - 1)
    dy = (y1 - y0)/(n - 1)
    dz = (z1 - z0)/(n - 1)

    do i = 1, n
        xArr(i) = x0 + (i-1)*dx
        yArr(i) = y0 + (i-1)*dy
        zArr(i) = z0 + (i-1)*dz
    end do

end subroutine linspace3d

!-------------------------------------------------------------------------------
! stell_point_xform(nfp, xform_hp, xIn, yIn, zIn, xOut, yOut, zOut)
!
! Transforms a point along a given number of consecutive half-periods in a 
! stellarator-symmetric configuration.
!
! Input parameters:
!     integer  :: nfp              -> # field periods in the configuration
!     integer  :: xform_hp         -> # half-periods along which to translate
!     real(dp) :: xIn, yIn, zIn    -> x, y, z coordinates of the origin point
!
! Output parameters:
!     real(dp) :: xOut, yOut, zOut -> x, y, z coords of the transformed point
!-------------------------------------------------------------------------------
subroutine stell_point_xform(nfp, xform_hp, xIn, yIn, zIn, xOut, yOut, zOut)

    implicit none

    integer, intent(IN) :: nfp, xform_hp
    real(dp), intent(IN) :: xIn, yIn, zIn
    real(dp), intent(OUT) :: xOut, yOut, zOut
    real(dp) :: rIn, phiIn, rOut, phiOut, trSign

    ! Convert to cylindrical coordinates
    rIn   = sqrt(xIn**2 + yIn**2)
    phiIn = atan2(yIn, xIn)

    ! r coordinate is always preserved
    rOut = rIn

    ! Transform the phi and z coordinate
    call stell_phi_z_xform(nfp, xform_hp, phiIn, zIn, phiOut, zOut, trSign)

    ! Transform back to Cartesian coordinates
    xOut = rOut * cos(phiOut)
    yOut = rOut * sin(phiOut)

end subroutine stell_point_xform

!-------------------------------------------------------------------------------
! stell_vector_xform(nfp, xform_hp, xIn, yIn, zIn, vxIn, vyIn, vzIn, 
!                    xOut, yOut, zOut, vxOut, vyOut, vzOut)
!
! Transforms a vector with a defined origin point along a given number of 
! consecutive half-periods in a stellarator-symmetric configuration.
!
! Input parameters:
!     integer  :: nfp              -> # field periods in the configuration
!     integer  :: xform_hp         -> # half-periods along which to translate
!     real(dp) :: xIn, yIn, zIn    -> x, y, z coordinates of the origin point
!     real(dp) :: vxIn, vyIn, vzIn -> x, y, z components of the vector
!
! Output parameters:
!     real(dp) :: xOut, yOut, zOut    -> x, y, z coordinates, origin of 
!                                        the transformed vector
!     real(dp) :: vxOut, vyOut, vzOut -> x, y, z components of the transformed
!                                        vector
!-------------------------------------------------------------------------------
subroutine stell_vector_xform(nfp, xform_hp, xIn, yIn, zIn, vxIn, vyIn, vzIn, &
                              xOut, yOut, zOut, vxOut, vyOut, vzOut)

    implicit none

    integer, intent(IN) :: nfp, xform_hp
    real(dp), intent(IN) :: xIn, yIn, zIn, vxIn, vyIn, vzIn
    real(dp), intent(OUT) :: xOut, yOut, zOut, vxOut, vyOut, vzOut
    real(dp) ::   rhat_x_in,   rhat_y_in,   rhat_x_out,   rhat_y_out 
    real(dp) :: phihat_x_in, phihat_y_in, phihat_x_out, phihat_y_out 
    real(dp) :: rIn, phiIn, vrIn, vphiIn, rOut, phiOut, vrOut, vphiOut, trSign

    ! Radial coordinate and toroidal angle of input point
    rIn   = sqrt(xIn**2 + yIn**2)
    phiIn = atan2(yIn, xIn)

    ! Cylindrical unit vectors at the input toroidal angle
    rhat_x_in = cos(phiIn)
    rhat_y_in = sin(phiIn)
    phihat_x_in = -sin(phiIn)
    phihat_y_in = cos(phiIn)

    ! Radial and toroidal components of the input vector
    vrIn   = vxIn*rhat_x_in   + vyIn*rhat_y_in
    vphiIn = vxIn*phihat_x_in + vyIn*phihat_y_in

    ! r coordinate is always preserved
    rOut = rIn

    ! Transform the toroidal angle and z coordinate of the origin point
    call stell_phi_z_xform(nfp, xform_hp, phiIn, zIn, phiOut, zOut, trSign)

    ! Transformed x and y components of the origin point
    xOut = rOut * cos(phiOut)
    yOut = rOut * sin(phiOut)

    ! Cylindrical unit vectors at the output toroidal angle
    rhat_x_out = cos(phiOut)
    rhat_y_out = sin(phiOut)
    phihat_x_out = -sin(phiOut)
    phihat_y_out = cos(phiOut)

    ! Radial and toroidal components of the output vector
    vrOut = vrIn
    vphiOut = trSign*vphiIn

    ! Cartesian components of the output vector
    vxOut = vrOut*rhat_x_out + vphiOut*phihat_x_out
    vyOut = vrOut*rhat_y_out + vphiOut*phihat_y_out
    vzOut = trSign*vzIn

end subroutine stell_vector_xform

!-------------------------------------------------------------------------------
! stell_phi_z_xform(nfp, xform_hp, phiIn, zIn, phiOut, zOut, trSign)
!
! Transforms the toroidal-angle (phi) and vertical (z) coordinates of a point
! along a given number of consecutive half-periods in a stellarator-symmetric
! configuration.
!
! Note: the poloidal angle may used in place of the z coordinate, as the 
! the same transformation applies.
!
! Input parameters:
!     integer  :: nfp      -> number of field periods in the configuration
!     integer  :: xform_hp -> number of half-periods along which to translate
!     real(dp) :: phiIn    -> toroidal (phi) coordinate of input point (radians)
!     real(dp) :: zIn      -> vertical (z) coordinate of the input point
!
! Output parameters:
!     real(dp) :: phiOut   -> phi coordinate of the transformed point (radians)
!     real(dp) :: zOut     -> z coordinate of the transformed point
!     real(dp) :: trSign   -> -1 if the transformed coordinates reversed sign;
!                             otherwise +1
!-------------------------------------------------------------------------------
subroutine stell_phi_z_xform(nfp, xform_hp, phiIn, zIn, phiOut, zOut, trSign)

    use magpie_globals, only: pi

    implicit none

    integer, intent(IN) :: nfp, xform_hp
    real(dp), intent(IN) :: phiIn, zIn
    real(dp), intent(OUT) :: phiOut, zOut, trSign
    real(dp) :: dPhi_hp, phi_plane, phi_refl, dPhi_refl

    ! Toroidal angle within one half-period
    dPhi_hp = pi / real(nfp, dp)

    ! Determine whether the xform is a whole number of field periods (even # hp)
    if (mod(xform_hp, 2) == 0) then

        ! For whole-period transforms, simply push everything along in phi
        phiOut = phiIn + xform_hp*dPhi_hp
        zOut   = zIn
        trSign  = 1.0

    else

        ! For half-period transforms, reflect in half-period plane
        phi_plane = floor(phiIn/dPhi_hp) * dPhi_hp
        dPhi_refl = phiIn - phi_plane
        phi_refl = phi_plane - dPhi_refl
        phiOut = phi_refl + (xform_hp + sign(real(1,dp), dPhi_refl))*dPhi_hp

        zOut   = -zIn
        trSign  = -1.0

    end if

end subroutine stell_phi_z_xform

!-------------------------------------------------------------------------------
! stell_point_transform(mode, phi, xIn, yIn, zIn, xOut, yOut, zOut)
!
! Transforms a point in one of two ways, depending on the mode selected:
!     reflect:   Reflects a point in a given poloidal plane presumed to form 
!                the boundary between two stellarator-symmetryc half-periods. 
!                The output point will have the equivalent location associated 
!                with the adjacent half-period. 
!     translate: Translates a point in the toroidal direction (i.e., along
!                a circle with a fixed radius about the z axis) by a given
!                angle.
!
! Input parameters:
!     character(len=*) :: mode     -> 'reflect' or 'translate' (see above)
!     real(dp) :: phi              -> 'reflect' mode: toroidal angle (radians)
!                                         of the symmetry plane
!                                     'translate' mode: toroidal interval
!                                         (radians) along which to translate
!     real(dp) :: xIn, yIn, zIn    -> x, y, z coordinates of the origin point
!     real(dp) :: vxIn, vyIn, vzIn -> x, y, z components of the vector
!
! Output parameters:
!     real(dp) :: xOut, yOut, zOut    -> x, y, z coordinates, origin of 
!                                        the reflected vector
!     real(dp) :: vxOut, vyOut, vzOut -> x, y, z components of the reflected
!                                        vector
!-------------------------------------------------------------------------------
subroutine stell_point_transform(mode, phi, xIn, yIn, zIn, xOut, yOut, zOut)

    implicit none

    character(len=*), intent(IN) :: mode
    real(dp), intent(IN) :: phi, xIn, yIn, zIn
    real(dp), intent(OUT) :: xOut, yOut, zOut
    real(dp) :: rIn, phiIn, rOut, phiOut
    logical :: refl

    if (trim(mode) == 'reflect') then
        refl = .true.
    else if (trim(mode) == 'translate') then
        refl = .false.
    else
        write(*,*) 'stell_point_transform: unrecognized mode ' // &
                   '(must be ''translate'' or ''reflect'')'
        stop
    end if

    ! Radial coordinate and toroidal angle of the input point
    rIn = sqrt(xIn**2 + yIn**2)
    phiIn = atan2(yIn, xIn)

    ! Radial coordinate and toroidal angle of the output point
    rOut = rIn
    if (refl) then
        phiOut = phi + (phi - phiIn)
        zOut   = -zIn
    else
        phiOut = phiIn + phi
        zOut   = zIn
    end if

    ! Cartesian coordinates of the output point
    xOut = rOut * cos(phiOut)
    yOut = rOut * sin(phiOut)

end subroutine stell_point_transform

!-------------------------------------------------------------------------------
! stell_vector_transform(mode, phi, xIn, yIn, zIn, vxIn, vyIn, vzIn, 
!                        xOut, yOut, zOut, vxOut, vyOut, vzOut)
!
! Transforms a vector in one of two ways, depending on the mode selected:
!     reflect:   Reflects a vector with a defined origin point in a given 
!                poloidal plane presumed to form the boundary between two 
!                stellarator-symmetryc half-periods. The output vector will 
!                have the equivalent origin point and direction associated with
!                the adjacent half-period. Vector magnitude is preserved.
!     translate: Translates a vector with a defined origin in the toroidal 
!                direction. The origin thus moves along a circle of fixed 
!                radius about the z axis by a given angle. The azimuthal 
!                components of the vector change in order to preserve the
!                radial and toroidal components. The vertical component of the
!                vector, as well as its length, are held fixed.
!
! Input parameters:
!     character(len=*) :: mode     -> 'reflect' or 'translate' (see above)
!     real(dp) :: phi              -> 'reflect' mode: toroidal angle (radians)
!                                         of the symmetry plane
!                                     'translate' mode: toroidal interval
!                                         (radians) along which to translate
!     real(dp) :: xIn, yIn, zIn    -> x, y, z coordinates of the origin point
!     real(dp) :: vxIn, vyIn, vzIn -> x, y, z components of the vector
!
! Output parameters:
!     real(dp) :: xOut, yOut, zOut    -> x, y, z coordinates, origin of 
!                                        the reflected vector
!     real(dp) :: vxOut, vyOut, vzOut -> x, y, z components of the reflected
!                                        vector
!-------------------------------------------------------------------------------
subroutine stell_vector_transform(mode, phi, xIn, yIn, zIn, vxIn, vyIn, vzIn, &
                                  xOut, yOut, zOut, vxOut, vyOut, vzOut)

    implicit none

    character(len=*), intent(IN) :: mode
    real(dp), intent(IN) :: phi, xIn, yIn, zIn, vxIn, vyIn, vzIn
    real(dp), intent(OUT) :: xOut, yOut, zOut, vxOut, vyOut, vzOut
    real(dp) :: rIn, vrIn, phiIn, vphiIn, rOut, vrOut, phiOut, vphiOut
    real(dp) ::  rhat_x_in,  rhat_y_in,  phihat_x_in,  phihat_y_in
    real(dp) :: rhat_x_out, rhat_y_out, phihat_x_out, phihat_y_out
    logical :: refl

    if (trim(mode) == 'reflect') then
        refl = .true.
    else if (trim(mode) == 'translate') then
        refl = .false.
    else
        write(*,*) 'stell_vector_transform: unrecognized mode ' // &
                   '(must be ''translate'' or ''reflect'')'
        stop
    end if

    ! Radial coordinate and toroidal angle of the input origin point
    rIn = sqrt(xIn**2 + yIn**2)
    phiIn = atan2(yIn, xIn)

    ! Radial and toroidal unit vectors at the origin point
    rhat_x_in   =  cos(phiIn)
    rhat_y_in   =  sin(phiIn)
    phihat_x_in = -sin(phiIn)
    phihat_y_in =  cos(phiIn)

    ! Radial and toroidal components of the input vector
    vrIn   = vxIn*rhat_x_in   + vyIn*rhat_y_in
    vphiIn = vxIn*phihat_x_in + vyIn*phihat_y_in

    ! Output origin coordinates & vector components; radial, toroidal, vertical
    rOut    = rIn
    vrOut   = vrIn
    if (refl) then
        phiOut  = phi + (phi - phiIn)
        vphiOut = -vphiIn
        zOut    = -zIn
        vzOut   = -vzIn
    else
        phiOut  = phi + phiIn
        vphiOut = vphiIn
        zOut    = zIn
        vzOut   = vzIn
    end if

    ! Radial and toroidal unit vectors at the output origin point
    rhat_x_out   =  cos(phiOut)
    rhat_y_out   =  sin(phiOut)
    phihat_x_out = -sin(phiOut)
    phihat_y_out =  cos(phiOut)

    ! Cartesian x/y coordinates of the output origin point
    xOut = rOut * cos(phiOut)
    yOut = rOut * sin(phiOut)   

    ! Cartesian x/y components of the output vector
    vxOut = vrOut*rhat_x_out + vphiOut*phihat_x_out
    vyOut = vrOut*rhat_y_out + vphiOut*phihat_y_out

end subroutine stell_vector_transform

end module geometry


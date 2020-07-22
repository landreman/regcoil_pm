module output

use magpie_globals, only: dp

implicit none

contains

!-------------------------------------------------------------------------------
! write_focus_file(magnets, nMagnets, base_name)
!
! Prints magnet information to a file suitable for input to the FOCUS code for
! optimization.
!
! Input paramters:
!     type(magnet), dimension(:) :: magnets -> array of magnet structures
!     integer :: nMagnets                   -> number of magnets in the array
!     character(len=len_fname) :: base_name -> output filename w/o extension
!-------------------------------------------------------------------------------
subroutine write_focus_file(magnets, nMagnets, base_name)

    use magpie_globals,    only: unit_focus, len_fname, focus_type, suffix_focus
    use magnet_properties, only: magnet

    implicit none

    type(magnet), dimension(:), intent(IN) :: magnets
    integer, intent(IN) :: nMagnets
    character(len=len_fname), intent(IN) :: base_name
    character(len=len_fname) :: filename
    integer :: i, open_status, write_status, rho_free, axis_free

    ! Open the file
    filename = trim(base_name) // trim(suffix_focus)
    open(unit=unit_focus, file=trim(filename), action='write', &
         iostat=open_status)
    if (open_status /= 0) then
        write(*, *) 'write_focus_file: unable to open file ' // trim(filename)
        stop
    end if

    write(*,fmt='(A)') '    Writing magnet data to file ' // trim(filename)

    ! Header lines
    write(unit=unit_focus, fmt='(A)') '# Total number of dipoles '
    write(unit=unit_focus, fmt='(I0)', iostat=write_status) nMagnets
    if (write_status /= 0) then
        write(*, *) 'write_focus_file: unable to write ' &
                    // 'first line of ' // filename
        stop
    end if

    write(unit=unit_focus, fmt='(A)') &
         '# coiltype, symmetry,  coilname,  ox,  oy,  oz,  Ic,  M_0,  ' &
         // 'rho,  Lc,  mp,  mt '
    
    ! Lines of magnet data
    do i = 1, nMagnets

        rho_free  = merge(1, 0, magnets(i)%moment_free)
        axis_free = merge(1, 0, magnets(i)%axis_free)

        write(unit=unit_focus, &
              fmt='(I0, A, I0, 3A, ES15.8, A, ES15.8, A, ES15.8, A, '      // &
                  ' I0, A, ES15.8, A, ES15.8, A, I0, A, ES15.8, A, ES15.8)', &
              iostat=write_status) &
                  focus_type,                ', ', &
                  magnets(i)%symm,           ', ', &
                  trim(magnets(i)%label),    ', ', &
                  magnets(i)%ox,             ', ', &
                  magnets(i)%oy,             ', ', &
                  magnets(i)%oz,             ', ', &
                  rho_free,                  ', ', &
                  magnets(i)%moment,         ', ', &
                  magnets(i)%rho_initial,    ', ', &
                  axis_free,                 ', ', &
                  magnets(i)%phi,            ', ', &
                  magnets(i)%theta
        if (write_status /= 0) then
            write(*, fmt='(A, I0, A)') &
                'write_focus_file: unable to write line ', i, &
                ' to ' // trim(filename)
            stop
        end if
    end do

    close(unit_focus)

end subroutine write_focus_file

!-------------------------------------------------------------------------------
! write_corners_file(qhexes, nQhexes, incl_err, incl_ovl, ovl, base_name)
!
! Writes a file containing the Cartesian coordinates of the corners of an
! array of quadrilaterally-faced hexahedra. Each line contains the x, y, and
! z cordinates of each of the eight corners of a single hexahedron in the 
! following order (24 numbers per line altogether):
!     ctx1, cty1, ctz1, ctx2, cty2, ctz2, ctx3, cty3, ctz3, ctx4, cty4, ctz4,
!     cbx1, cby1, cbz1, cbx2, cby2, cbz2, cbx3, cby3, cbz3, cbx4, cby4, cbz4
!
! For further information on the orientation of the corners as defined above,
! see the documentation in qhex_properties.f90.
!
! Input parameters:
!     type(quad_hexahedron), dimension(:) :: qhexes
!                            -> Array of qhex structures to be printed to file
!     integer :: nQhexes     -> Number of qhexes (incl. erroneous) in array
!     logical :: incl_err    -> If false, erroneous qhexes won't be printed
!     logical :: incl_ovl    -> If false, qhexes with overlaps won't be printed
!     logical, dimension(nQhexes) :: ovl 
!                            -> indicates which qhexes have overlaps
!     character :: base_name -> base name of the file (not including suffix)
!-------------------------------------------------------------------------------
subroutine write_corners_file(qhexes, nQhexes, incl_err, incl_ovl, ovl, &
                              base_name)

    use magpie_globals,    only: unit_corners, len_fname, suffix_corners
    use qhex_properties,   only: quad_hexahedron

    implicit none

    type(quad_hexahedron), dimension(:), intent(IN) :: qhexes
    integer, intent(IN) :: nQhexes
    logical, intent(IN) :: incl_err, incl_ovl
    logical, dimension(nQhexes), intent(IN) :: ovl
    character(len=len_fname), intent(IN) :: base_name
    character(len=len_fname) :: filename
    integer :: i, open_status, write_status

    ! Open the file
    filename = trim(base_name) // trim(suffix_corners)
    open(unit=unit_corners, file=trim(filename), action='write', &
         iostat=open_status)
    if (open_status /= 0) then
        write(*, *) 'write_corners_file: unable to open file ' // trim(filename)
        stop
    end if

    write(*,fmt='(A)') '    Writing corners data to file ' // trim(filename)

    ! Write data on the lower qhexes
    do i = 1, nQhexes

        if ((incl_err .or. (.not. qhexes(i)%err)) .and. &
            (incl_ovl .or. (.not. ovl(i))       )        ) then
            call write_corners_file_line(unit_corners, qhexes(i), write_status)
            if (write_status /= 0) then
                write(*, fmt='(A, I0, A)') &
                    'write_corners_file: unable to write line ', i, &
                    ' to file ' // trim(filename)
                stop
            end if
        end if

    end do

    close(unit_corners)

end subroutine write_corners_file

!-------------------------------------------------------------------------------
! write_corners_file_line(file_unit, qhex, write_status)
!
! Writes one line of a corners file, containing the Cartesian coordinates of the
! corners of one quadrilaterally-faced hexahedron.
!
! Input paramters:
!     integer :: file_unit          -> Reference for the file to be written to
!     type(quad_hexahedron) :: qhex -> qhex whose information is to be written
!
! Output parameters:
!     integer :: write_status  -> Status upon completion of the write statement
!-------------------------------------------------------------------------------
subroutine write_corners_file_line(file_unit, qhex, write_status)

    use qhex_properties, only: quad_hexahedron

    implicit none

    integer, intent(IN) :: file_unit
    type(quad_hexahedron), intent(IN) :: qhex
    integer, intent(OUT) :: write_status

    write(unit=file_unit, &
        fmt='(ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
            ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
            ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
            ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
            ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X )',           &
        iostat=write_status) &
        qhex%ctx1, qhex%cty1, qhex%ctz1, qhex%ctx2, qhex%cty2, qhex%ctz2, &
        qhex%ctx3, qhex%cty3, qhex%ctz3, qhex%ctx4, qhex%cty4, qhex%ctz4, &
        qhex%cbx1, qhex%cby1, qhex%cbz1, qhex%cbx2, qhex%cby2, qhex%cbz2, &
        qhex%cbx3, qhex%cby3, qhex%cbz3, qhex%cbx4, qhex%cby4, qhex%cbz4

end subroutine write_corners_file_line

!-------------------------------------------------------------------------------
! write_base_grid_file(bg, base_name)
!
! Prints base grid data to a text file.
!
! Each line stores the following data:
!     xVert, yVert, zVert, xSurf, ySurf, zSurf, nx, ny, nz
! where [x,y,z]Vert are the coordinates of the vertices of the base grid,
! [x,y,z]Surf are the points on the reference surface corresponding to the 
! base grid vertices, and n[x,y,z] are unit normal vectors pointing from
! the surface points to the grid points.
!
! Input paramters:
!     type(base_grid_struct)   :: bg -> base grid whose data is to be printed
!     character(len=len_fname) :: base_name -> output filename w/o extension
!-------------------------------------------------------------------------------
subroutine write_base_grid_file(bg, base_name)

    use magpie_globals,    only: unit_bg, len_fname, suffix_bg
    use base_grid,         only: base_grid_struct

    implicit none

    type(base_grid_struct), intent(IN) :: bg
    character(len=len_fname), intent(IN) :: base_name
    character(len=len_fname) :: filename
    integer :: i, j, open_status, write_status

    ! Open the file
    filename = trim(base_name) // trim(suffix_bg)
    open(unit=unit_bg, file=trim(filename), action='write', &
         iostat=open_status)
    if (open_status /= 0) then
        write(*, *) 'write_base_grid_file: unable to open file ' &
                    // trim(filename)
        stop
    end if

    write(*,fmt='(A)') '    Writing base grid data to file ' // trim(filename)

    do i = 1, bg%nTheta
        do j = 1, bg%nPhi
            write(unit=unit_bg,                                &
                  fmt='(ES15.8, X, ES15.8, X, ES15.8, X, ' //  &
                       'ES15.8, X, ES15.8, X, ES15.8, X, ' //  &
                       'ES15.8, X, ES15.8, X, ES15.8)',        &
                  iostat=write_status                        ) &
                bg%x(i,j),     bg%y(i,j),     bg%z(i,j),       &
                bg%xSurf(i,j), bg%ySurf(i,j), bg%zSurf(i,j),   &
                bg%nx(i,j),    bg%ny(i,j),    bg%nz(i,j)
            if (write_status /= 0) then
                write(*, fmt='(A, I0, A)') &
                    'write_base_grid_file: unable to write line ', i, &
                    ' to ' // trim(filename)
            end if
        end do
    end do

    close(unit_bg)

end subroutine write_base_grid_file

!-------------------------------------------------------------------------------
! write_cbrick_file(cbricks, nCbricks, incl_ovl, ovl, base_name)
!
! Writes a file containing the bounding values of the radial, vertical, and 
! toroidal coordinates of an input set of curved bricks represented by cbrick
! structures. Each line contains the following data:
!     ri, ro, zl, zu, phib, phif
! ri and ro are the radial coordinate (m) of the inboard and outboard faces;
! zl and zu are the vertical coordinates (m) of the lower and upper faces;
! phib and phif are the toroidal angles (rad) of the back and front faces.
!
! Input parameters:
!     type(cbrick), dimension(:) :: cbricks
!                            -> Array of cbrick structures to be printed to file
!     integer :: nCbricks    -> Number of cbricks in the array
!     logical :: incl_ovl    -> If false, cbricks with overlaps won't be printed
!     logical, dimension(nQhexes) :: ovl 
!                            -> indicates which cbricks have overlaps
!     character :: base_name -> base name of the file (not including suffix)
!-------------------------------------------------------------------------------
subroutine write_cbrick_file(cbricks, nCbricks, incl_overlap, ovl, base_name)

    use magpie_globals,    only: unit_cbrick, len_fname, suffix_cbrick
    use cbrick_properties, only: cbrick

    implicit none

    type(cbrick), dimension(:), intent(IN) :: cbricks
    integer, intent(IN) :: nCbricks
    logical, intent(IN) :: incl_overlap
    logical, dimension(:), intent(IN) :: ovl
    character(len=len_fname), intent(IN) :: base_name
    character(len=len_fname) :: filename
    integer :: i, open_status, write_status

    ! Open the file
    filename = trim(base_name) // trim(suffix_cbrick)
    open(unit=unit_cbrick, file=trim(filename), action='write', &
         iostat=open_status)
    if (open_status /= 0) then
        write(*, *) 'write_cbrick_file: unable to open file ' // trim(filename)
        stop
    end if

    write(*,fmt='(A)') &
        '    Writing cbrick geometric data to file ' // trim(filename)

    ! Write data on the lower qhexes
    do i = 1, nCbricks
        if (incl_overlap .or. (.not. ovl(i))) then
            write(unit_cbrick, fmt='(ES15.8, X, ES15.8, X, ES15.8, X, ' //    &
                                   ' ES15.8, X, ES15.8, X, ES15.8)',          &
                  iostat=write_status)                                        &
                cbricks(i)%ri,   cbricks(i)%ro, cbricks(i)%zl, cbricks(i)%zu, &
                cbricks(i)%phib, cbricks(i)%phif
            if (write_status /= 0) then
                write(*, fmt='(A, I0, A)') &
                    'write_corners_file: unable to write line ', i, &
                    ' to file ' // trim(filename)
                stop
            end if
        end if

    end do

    close(unit_cbrick)

end subroutine write_cbrick_file

!-------------------------------------------------------------------------------
! write_cbrick_vtk(cbricks, nCbricks, incl_ovl, ovl, base_name)
!
! Writes a VTK file with geometric data for a three-dimensional solid 
! representation of each curved-brick magnet in an input set. Also contains
! two scalar datasets containing the r and z coordinates, respectively, of
! the center-of-mass of each brick.
!
! Input parameters:
!     type(cbrick), dimension(:) :: cbricks
!                            -> Array of cbrick structures to be printed to file
!     integer :: nCbricks    -> Number of cbricks in the array
!     logical :: incl_ovl    -> If false, cbricks with overlaps won't be printed
!     logical, dimension(nQhexes) :: ovl 
!                            -> indicates which cbricks have overlaps
!     character :: base_name -> base name of the file (not including suffix)
!-------------------------------------------------------------------------------
subroutine write_cbrick_vtk(cbricks, nCbricks, incl_overlap, ovl, base_name)

    use magpie_globals,    only: unit_vtk, len_fname, suffix_mag_vtk, &
                                 arc_min_pts_per_deg, arcs_per_cbrick, pi
    use cbrick_properties, only: cbrick

    implicit none

    type(cbrick), dimension(:), intent(IN) :: cbricks
    integer, intent(IN) :: nCbricks
    logical, intent(IN) :: incl_overlap
    logical, dimension(:), intent(IN) :: ovl
    character(len=len_fname), intent(IN) :: base_name
    character(len=len_fname) :: filename
    integer :: i, j, open_status, write_status
    integer :: nPts_arc, nPts_total, nPoly_total
    integer, dimension(nCbricks) :: nPtsArc, nPts, nPoly
    integer, dimension(nCbricks) :: indInLo, indOutLo, indInUp, indOutUp
    real(dp) :: deltaPhi_deg, phiStep, phi

    ! Open the file
    filename = trim(base_name) // trim(suffix_mag_vtk)
    open(unit=unit_vtk, file=trim(filename), action='write', &
         iostat=open_status)
    if (open_status /= 0) then
        write(*, *) 'write_cbrick_vtk: unable to open file ' // trim(filename)
        stop
    end if

    write(*,fmt='(A)') &
        '    Writing cbrick geometric data to VTK file ' // trim(filename)

    ! Header lines
    write(unit=unit_vtk, fmt='(A)', iostat=write_status) &
        '# vtk DataFile Version 2.0' // new_line('a') // &
        trim(filename)               // new_line('a') // &
        'ASCII'                      // new_line('a') // &
        'DATASET POLYDATA'
    if (write_status /= 0) then
        write(*, *) 'write_cbrick_vtk: unable to write header lines of ' &
                    // filename
        stop
    end if

    ! Determine the number of points, polygons, and lines for each brick
    nPts_total = 0
    nPoly_total = 0
    indInLo(1) = 0
    do i = 1, nCbricks
        if (.not. ovl(i) .or. incl_overlap) then
            deltaPhi_deg = (cbricks(i)%phif - cbricks(i)%phib)*180./pi
            nPtsArc(i) = ceiling(deltaPhi_deg/real(arc_min_pts_per_deg,dp))
            nPts(i) = arcs_per_cbrick*nPtsArc(i)
            indOutLo(i) = indInLo(i)  + nPtsArc(i)
            indInUp(i)  = indOutLo(i) + nPtsArc(i)
            indOutUp(i) = indInUp(i)  + nPtsArc(i)
            if (i < nCbricks) indInLo(i+1) = indOutUp(i) + nPtsArc(i)
            nPoly(i) = arcs_per_cbrick*(nPtsArc(i)-1) + 2
            nPts_total = nPts_total + nPts(i)
            nPoly_total = nPoly_total + nPoly(i)
        else
            nPts(i) = 0
            if (i < nCbricks) indInLo(i+1) = indInLo(i)
            indOutLo(i) = 0
            indOutUp(i) = 0
            indInLo(i)  = 0
            indInUp(i)  = 0
        end if
    end do

    ! Print point data to file
    write(unit=unit_vtk, fmt='(A,I0,A)', iostat=write_status) &
        'POINTS ', nPts_total, ' double'
    if (write_status /= 0) then
        write(*, *) 'write_cbrick_vtk: unable to write point header of ' &
                    // filename
        stop
    end if
    do i = 1, nCbricks

        if (ovl(i) .and. .not. incl_overlap) cycle

        phiStep = (cbricks(i)%phif - cbricks(i)%phib)/real(nPtsArc(i)-1,dp)

        ! Inner, lower arc
        do j = 1, nPtsArc(i)
            phi = cbricks(i)%phib + (j-1)*phiStep
            write(unit=unit_vtk, fmt='(ES16.8,X,ES16.8,X,ES16.8)') &
                cbricks(i)%ri*cos(phi), cbricks(i)%ri*sin(phi), cbricks(i)%zl
        end do

        ! Outer, lower arc
        do j = 1, nPtsArc(i)
            phi = cbricks(i)%phib + (j-1)*phiStep
            write(unit=unit_vtk, fmt='(ES16.8,X,ES16.8,X,ES16.8)') &
                cbricks(i)%ro*cos(phi), cbricks(i)%ro*sin(phi), cbricks(i)%zl
        end do

        ! Inner, upper arc
        do j = 1, nPtsArc(i)
            phi = cbricks(i)%phib + (j-1)*phiStep
            write(unit=unit_vtk, fmt='(ES16.8,X,ES16.8,X,ES16.8)') &
                cbricks(i)%ri*cos(phi), cbricks(i)%ri*sin(phi), cbricks(i)%zu
        end do

        ! Outer, upper arc
        do j = 1, nPtsArc(i)
            phi = cbricks(i)%phib + (j-1)*phiStep
            write(unit=unit_vtk, fmt='(ES16.8,X,ES16.8,X,ES16.8)') &
                cbricks(i)%ro*cos(phi), cbricks(i)%ro*sin(phi), cbricks(i)%zu
        end do

    end do

    ! Print polygon data to file
    write(unit=unit_vtk, fmt='(A,I0,X,I0)', iostat=write_status) & 
        'POLYGONS ', nPoly_total, 5*nPoly_total
    if (write_status /= 0) then
        write(*, *) 'write_cbrick_vtk: unable to write polygon header of ' &
                    // filename
        stop
    end if
    do i = 1, nCbricks

        if (ovl(i) .and. .not. incl_overlap) cycle
        if (ovl(i)) stop

        do j = 0, nPtsArc(i)-2
            write(unit=unit_vtk, fmt='(I0,X,I0,X,I0,X,I0,X,I0,X)') 4, &
                indOutLo(i)+j, indOutLo(i)+j+1, indOutUp(i)+j+1, indOutUp(i)+j
            write(unit=unit_vtk, fmt='(I0,X,I0,X,I0,X,I0,X,I0,X)') 4, &
                indOutUp(i)+j, indOutUp(i)+j+1, indInUp(i)+j+1,  indInUp(i)+j
            write(unit=unit_vtk, fmt='(I0,X,I0,X,I0,X,I0,X,I0,X)') 4, &
                indInUp(i)+j,  indInUp(i)+j+1,  indInLo(i)+j+1,  indInLo(i)+j
            write(unit=unit_vtk, fmt='(I0,X,I0,X,I0,X,I0,X,I0,X)') 4, &
                indInLo(i)+j,  indInLo(i)+j+1,  indOutLo(i)+j+1, indOutLo(i)+j
        end do
        write(unit=unit_vtk, fmt='(I0,X,I0,X,I0,X,I0,X,I0,X)') 4, &
            indOutLo(i), indOutUp(i), indInUp(i), indInLo(i)
        write(unit=unit_vtk, fmt='(I0,X,I0,X,I0,X,I0,X,I0,X)') 4, &
            indOutLo(i)+nPtsArc(i)-1, indOutUp(i)+nPtsArc(i)-1,     &
            indInUp(i) +nPtsArc(i)-1, indInLo(i) +nPtsArc(i)-1
    end do

    ! Print test data
    write(unit=unit_vtk, fmt='(A,I0)') 'CELL_DATA ', nPoly_total
    write(unit=unit_vtk, fmt='(A)') 'SCALARS r double 1'
    write(unit=unit_vtk, fmt='(A)') 'LOOKUP_TABLE default'
    do i = 1, nCbricks
        if (ovl(i) .and. .not. incl_overlap) cycle
        do j = 1, nPoly(i)
            write(unit=unit_vtk, fmt='(ES16.8)') cbricks(i)%cmr
        end do
    end do
    write(unit=unit_vtk, fmt='(A)') 'SCALARS z double 1'
    write(unit=unit_vtk, fmt='(A)') 'LOOKUP_TABLE default'
    do i = 1, nCbricks
        if (ovl(i) .and. .not. incl_overlap) cycle
        do j = 1, nPoly(i)
            write(unit=unit_vtk, fmt='(ES16.8)') cbricks(i)%cmz
        end do
    end do

     close(unit_vtk)

end subroutine write_cbrick_vtk

!-------------------------------------------------------------------------------
! write_cbrick_wf_vtk(cbricks, nCbricks, incl_ovl, ovl, base_name)
!
! Writes a VTK file with geometric data for a three-dimensional wireframe
! representation of each curved-brick magnet in an input set. 
!
! Input parameters:
!     type(cbrick), dimension(:) :: cbricks
!                            -> Array of cbrick structures to be printed to file
!     integer :: nCbricks    -> Number of cbricks in the array
!     logical :: incl_ovl    -> If false, cbricks with overlaps won't be printed
!     logical, dimension(nQhexes) :: ovl 
!                            -> indicates which cbricks have overlaps
!     character :: base_name -> base name of the file (not including suffix)
!-------------------------------------------------------------------------------
subroutine write_cbrick_wf_vtk(cbricks, nCbricks, incl_overlap, ovl, base_name)

    use magpie_globals,    only: unit_wf_vtk, len_fname, suffix_wf_vtk, &
                                 arc_min_pts_per_deg, arcs_per_cbrick, pi
    use cbrick_properties, only: cbrick

    implicit none

    type(cbrick), dimension(:), intent(IN) :: cbricks
    integer, intent(IN) :: nCbricks
    logical, intent(IN) :: incl_overlap
    logical, dimension(:), intent(IN) :: ovl
    character(len=len_fname), intent(IN) :: base_name
    character(len=len_fname) :: filename
    integer :: i, j, open_status, write_status
    integer :: nPts_arc, nPts_total, nLine_total, lineSize
    integer, dimension(nCbricks) :: nPtsArc, nPts, nPoly, nLine
    integer, dimension(nCbricks) :: indInLo, indOutLo, indInUp, indOutUp
    real(dp) :: deltaPhi_deg, phiStep, phi

    ! Open the file
    filename = trim(base_name) // trim(suffix_wf_vtk)
    open(unit=unit_wf_vtk, file=trim(filename), action='write', &
         iostat=open_status)
    if (open_status /= 0) then
        write(*, *) 'write_cbrick_vtk: unable to open file ' // trim(filename)
        stop
    end if

    write(*,fmt='(A)') &
        '    Writing cbrick wireframe data to VTK file ' // trim(filename)

    ! Header lines
    write(unit=unit_wf_vtk, fmt='(A)', iostat=write_status) &
        '# vtk DataFile Version 2.0' // new_line('a') // &
        trim(filename)               // new_line('a') // &
        'ASCII'                      // new_line('a') // &
        'DATASET POLYDATA'
    if (write_status /= 0) then
        write(*, *) 'write_cbrick_wf_vtk: unable to write header lines of ' &
                    // filename
        stop
    end if

    ! Determine the number of points, polygons, and lines for each brick
    nPts_total = 0
    nLine_total = 0
    lineSize = 0
    indInLo(1) = 0
    do i = 1, nCbricks
        if (.not. ovl(i) .or. incl_overlap) then
            deltaPhi_deg = (cbricks(i)%phif - cbricks(i)%phib)*180./pi
            nPtsArc(i) = ceiling(deltaPhi_deg/real(arc_min_pts_per_deg,dp))
            nPts(i) = arcs_per_cbrick*nPtsArc(i)
            indOutLo(i) = indInLo(i)  + nPtsArc(i)
            indInUp(i)  = indOutLo(i) + nPtsArc(i)
            indOutUp(i) = indInUp(i)  + nPtsArc(i)
            if (i < nCbricks) indInLo(i+1) = indOutUp(i) + nPtsArc(i)
            nLine(i) = arcs_per_cbrick + 2
            nPts_total = nPts_total + nPts(i)
            nLine_total = nLine_total + nLine(i)
            lineSize = lineSize + arcs_per_cbrick*nPtsArc(i) + 10
        else
            nPts(i) = 0
            if (i < nCbricks) indInLo(i+1) = indInLo(i)
            indOutLo(i) = 0
            indOutUp(i) = 0
            indInLo(i)  = 0
            indInUp(i)  = 0
        end if
    end do

    ! Print point data to file
    write(unit=unit_wf_vtk, fmt='(A,I0,A)', iostat=write_status) &
        'POINTS ', nPts_total, ' double'
    if (write_status /= 0) then
        write(*, *) 'write_cbrick_vtk: unable to write point header of ' &
                    // filename
        stop
    end if
    do i = 1, nCbricks

        if (ovl(i) .and. .not. incl_overlap) cycle

        phiStep = (cbricks(i)%phif - cbricks(i)%phib)/real(nPtsArc(i)-1,dp)

        ! Inner, lower arc
        do j = 1, nPtsArc(i)
            phi = cbricks(i)%phib + (j-1)*phiStep
            write(unit=unit_wf_vtk, fmt='(ES16.8,X,ES16.8,X,ES16.8)') &
                cbricks(i)%ri*cos(phi), cbricks(i)%ri*sin(phi), cbricks(i)%zl
        end do

        ! Outer, lower arc
        do j = 1, nPtsArc(i)
            phi = cbricks(i)%phib + (j-1)*phiStep
            write(unit=unit_wf_vtk, fmt='(ES16.8,X,ES16.8,X,ES16.8)') &
                cbricks(i)%ro*cos(phi), cbricks(i)%ro*sin(phi), cbricks(i)%zl
        end do

        ! Inner, upper arc
        do j = 1, nPtsArc(i)
            phi = cbricks(i)%phib + (j-1)*phiStep
            write(unit=unit_wf_vtk, fmt='(ES16.8,X,ES16.8,X,ES16.8)') &
                cbricks(i)%ri*cos(phi), cbricks(i)%ri*sin(phi), cbricks(i)%zu
        end do

        ! Outer, upper arc
        do j = 1, nPtsArc(i)
            phi = cbricks(i)%phib + (j-1)*phiStep
            write(unit=unit_wf_vtk, fmt='(ES16.8,X,ES16.8,X,ES16.8)') &
                cbricks(i)%ro*cos(phi), cbricks(i)%ro*sin(phi), cbricks(i)%zu
        end do

    end do

    ! Print line data to file
    write(unit=unit_wf_vtk, fmt='(A,I0,X,I0)', iostat=write_status) & 
        'LINES ', nLine_total, lineSize+nLine_total
    if (write_status /= 0) then
        write(*, *) 'write_cbrick_vtk: unable to write line header of ' &
                    // filename
        stop
    end if
    do i = 1, nCbricks

        if (ovl(i) .and. .not. incl_overlap) cycle

        ! Outer, lower arc
        write(unit=unit_wf_vtk, fmt='(I0,X)', advance='no') nPtsArc(i)
        do j = 1, nPtsArc(i)
            write(unit=unit_wf_vtk, fmt='(I0,X)', advance='no') indOutLo(i)+j-1
        end do
        write(unit=unit_wf_vtk, fmt='(X)')

        ! Outer, upper arc
        write(unit=unit_wf_vtk, fmt='(I0,X)', advance='no') nPtsArc(i)
        do j = 1, nPtsArc(i)
            write(unit=unit_wf_vtk, fmt='(I0,X)', advance='no') indOutUp(i)+j-1
        end do
        write(unit=unit_wf_vtk, fmt='(X)')

        ! Inner, upper arc
        write(unit=unit_wf_vtk, fmt='(I0,X)', advance='no') nPtsArc(i)
        do j = 1, nPtsArc(i)
            write(unit=unit_wf_vtk, fmt='(I0,X)', advance='no') indInUp(i)+j-1
        end do
        write(unit=unit_wf_vtk, fmt='(X)')

        ! Inner, lower arc
        write(unit=unit_wf_vtk, fmt='(I0,X)', advance='no') nPtsArc(i)
        do j = 1, nPtsArc(i)
            write(unit=unit_wf_vtk, fmt='(I0,X)', advance='no') indInLo(i)+j-1
        end do
        write(unit=unit_wf_vtk, fmt='(X)')

        ! Back and front rectangles
        write(unit=unit_wf_vtk, fmt='(I0,X,I0,X,I0,X,I0,X,I0,X,I0,X)') &
            5, indOutLo(i), indOutUp(i), indInUp(i), indInLo(i), indOutLo(i)
        write(unit=unit_wf_vtk, fmt='(I0,X,I0,X,I0,X,I0,X,I0,X,I0,X)') &
            5, indOutLo(i) + nPtsArc(i)-1, indOutUp(i) + nPtsArc(i)-1,  &
               indInUp(i)  + nPtsArc(i)-1, indInLo(i)  + nPtsArc(i)-1,  &
               indOutLo(i) + nPtsArc(i)-1
     end do

     close(unit_wf_vtk)

end subroutine write_cbrick_wf_vtk

!-------------------------------------------------------------------------------
! write_trec_corners_file(trecs, nTrecs, incl_ovl, ovl, base_name)
!
! Writes a file containing the Cartesian coordinates of the corners of an
! array of rectangular prisms. The format is the same as for corners files for
! the more general case of quadrilaterally-faced hexahedra. Each line contains 
! the x, y, and z cordinates of each of the eight corners of a single prism in 
! the following order (24 numbers per line altogether):
!     ctx1, cty1, ctz1, ctx2, cty2, ctz2, ctx3, cty3, ctz3, ctx4, cty4, ctz4,
!     cbx1, cby1, cbz1, cbx2, cby2, cbz2, cbx3, cby3, cbz3, cbx4, cby4, cbz4
!
! For further information on the orientation of the corners as defined above,
! see the documentation in qhex_properties.f90.
!
! Input parameters:
!     type(trec), dimension(:) :: trecs
!                            -> Array of trec structures to be printed to file
!     integer :: nTrecs      -> Number of trecs in the array
!     logical :: incl_ovl    -> If false, trecs with overlaps won't be printed
!     logical, dimension(nTrecs) :: ovl 
!                            -> indicates which trecs have overlaps
!     character :: base_name -> base name of the file (not including suffix)
!-------------------------------------------------------------------------------
subroutine write_trec_corners_file(trecs, nTrecs, incl_ovl, ovl, base_name)

    use magpie_globals,    only: unit_corners, len_fname, suffix_corners
    use trec_properties,   only: trec

    implicit none

    type(trec), dimension(:), intent(IN) :: trecs
    integer, intent(IN) :: nTrecs
    logical, intent(IN) :: incl_ovl
    logical, dimension(nTrecs), intent(IN) :: ovl
    character(len=len_fname), intent(IN) :: base_name
    character(len=len_fname) :: filename
    integer :: i, open_status, write_status

    ! Open the file
    filename = trim(base_name) // trim(suffix_corners)
    open(unit=unit_corners, file=trim(filename), action='write', &
         iostat=open_status)
    if (open_status /= 0) then
        write(*, *) 'write_trec_corners_file: unable to open file ' &
                    // trim(filename)
        stop
    end if

    write(*,fmt='(A)') '    Writing corners data to file ' // trim(filename)

    ! Write data on the lower qhexes
    do i = 1, nTrecs

        if ( incl_ovl .or. (.not. ovl(i)) ) then
            call write_trec_corners_file_line(unit_corners, trecs(i), &
                                              write_status)
            if (write_status /= 0) then
                write(*, fmt='(A, I0, A)') &
                    'write_trec_corners_file: unable to write line ', i, &
                    ' to file ' // trim(filename)
                stop
            end if
        end if

    end do

    close(unit_corners)

end subroutine write_trec_corners_file

!-------------------------------------------------------------------------------
! write_trec_corners_file_line(file_unit, tr, write_status)
!
! Writes one line of a corners file, containing the Cartesian coordinates of the
! corners of one rectangular prism.
!
! Input paramters:
!     integer :: file_unit          -> Reference for the file to be written to
!     type(quad_hexahedron) :: tr   -> trec whose information is to be written
!
! Output parameters:
!     integer :: write_status  -> Status upon completion of the write statement
!-------------------------------------------------------------------------------
subroutine write_trec_corners_file_line(file_unit, tr, write_status)

    use trec_properties, only: trec

    implicit none

    integer, intent(IN) :: file_unit
    type(trec), intent(IN) :: tr
    integer, intent(OUT) :: write_status

    write(unit=file_unit, &
        fmt='(ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
            ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
            ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
            ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
            ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X )',           &
        iostat=write_status) &
        tr%xif, tr%yif, tr%zl, tr%xib, tr%yib, tr%zl, &
        tr%xib, tr%yib, tr%zu, tr%xif, tr%yif, tr%zu, &
        tr%xof, tr%yof, tr%zl, tr%xob, tr%yob, tr%zl, &
        tr%xob, tr%yob, tr%zu, tr%xof, tr%yof, tr%zu 

end subroutine write_trec_corners_file_line

end module output


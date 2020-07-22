module input

use magpie_globals, only: dp

implicit none

contains

!-------------------------------------------------------------------------------
! read_namelists(fname)
!
! Parses through a file and reads all namelists currently accepted by MAGPIE.
!
! Input parameter:
!     character(len=len_fname) :: fname -> name of the file to read
!-------------------------------------------------------------------------------
subroutine read_namelists(fname)

    use magpie_globals, only: magpie, len_fname, len_line, unit_in, magnet_type

    implicit none

    character(len=len_fname) :: fname
    integer :: open_stat, read_stat
    character(len=len_line) :: line
    logical :: magpie_read=.false.

    open(unit=unit_in, file=fname, status='old', action='read', &
         iostat=open_stat)
    if (open_stat /= 0) then
        write(*,*) 'Unable to open namelist file ' // trim(fname)
        stop
    end if

    ! Scan through the file and look for namelist declarations
    do while (.true.)

        read(unit=unit_in, fmt='(A)', iostat=read_stat) line

        ! Exit if end of file is reached
        if (read_stat < 0) then 
            if (.not. magpie_read) &
                write(*,fmt='(A)') &
                    '    Warning: no MAGPIE list found in file ', &
                    trim(fname), '; using only default values.'
            exit

        ! MAGPIE namelist
        else if (to_lowercase(trim(line)) == '&magpie') then
            write(*,fmt='(A)') '    Reading MAGPIE namelist from file'
            backspace(unit_in)
            read(unit=unit_in, nml=magpie, iostat=read_stat)
            if (read_stat /= 0) stop 'Problem reading magpie namelist'
            magpie_read = .true.

        end if

    end do

    ! Sanitize some string inputs
    magnet_type = to_lowercase(magnet_type)

    close(unit_in)

end subroutine read_namelists

!-------------------------------------------------------------------------------
! read_surface_file(fname, helicity, symm, surf)
!
! Reads in Fourier coefficients parametrizing the plasma lcfs from 
! a file. The file must have one row containing the integer number of modes N,
! followed by N rows with 4 or columns:
!     1. poloidal (m) mode number
!     2. toroidal (n) mode number
!     3. cosine coefficient, radial coordinate
!     4. sine coefficient, z coordinate
!     5. sine coefficient, r coordinate    (optional if stellarator symmetric)
!     6. cosine coefficient, z coordinate  (optional if stellarator symmetric)
! All other rows must begin with the commenting character (!)
!
! Input parameters
!     character :: fname     -> Path to the file with the table of coefficients.
!     integer   :: nfp       -> Number of field periods in the surface
!     integer   :: helicity  -> 1 if phase of modes is m*theta + nfp*n*phi; 
!                               -1 if phase is m*theta - nfp*n*phi
!     logical   :: symm      -> True if the surface is stellarator symmetric
!
! Output paramters:
!     type(surface) :: surface -> structure storing the parameters from the file
!-------------------------------------------------------------------------------
subroutine read_surface_file(fname, nfp, helicity, symm, surf)

    use magpie_globals, only: unit_surf, len_line
    use surface_calc,   only: surface

    implicit none

    character(len=*), intent(IN) :: fname
    integer, intent(IN) :: nfp, helicity
    logical, intent(IN) :: symm
    type(surface), intent(OUT) :: surf
    integer :: stat, nModes, i, m, n
    character(len=len_line) :: line, phase
    real(dp) :: rs, zc, rc, zs

    ! Check input
    if (.not. (helicity == 1 .or. helicity == -1)) &
        stop 'read_surface_file: helicity input must be -1 or 1'

    surf%symm = symm
    surf%nfp = nfp

    open(unit=unit_surf, file=trim(fname), iostat=stat, action='read')
    if (stat /= 0) stop 'read_surface_file: unable to open file ' // trim(fname)

    phase = merge('m*theta + nfp*n*phi', 'm*theta - nfp*n*phi', helicity > 0)
    write(*,fmt='(A)') '    Reading surface data from file ' // trim(fname) 
    write(*,fmt='(A)') '        assuming modes with phase ' // trim(phase)

    ! Read lines until the number of modes is obtained
    do while (.true.)
        read(unit=unit_surf, fmt=*, iostat=stat) nModes
        if (stat > 0)  cycle
        if (stat == 0) exit
        if (stat < 0) &
            stop 'read_surface_file: problem reading from file ' // trim(fname)
    end do

    ! Allocate arrays in the surface structure according to the number of modes
    surf%nModes = nModes
    allocate( surf%m(nModes), surf%n(nModes), surf%rc(nModes), surf%zs(nModes) )
    if (.not. symm) allocate( surf%rs(nModes), surf%zc(nModes) )

    ! Read the number of lines and allocate arrays accordingly
    i = 0
    do while (i < nModes)

        if (symm) then
            read(unit=unit_surf, fmt=*, iostat=stat) m, n, rc, zs
        else
            read(unit=unit_surf, fmt=*, iostat=stat) m, n, rc, zs, rs, zc
        end if

        if (stat > 0) cycle
        if (stat < 0) &
           stop 'read_surface_file: file ended before expected number of ' // &
                'modes read (' // trim(fname) // ')'

        i = i + 1
        surf%m(i)  = m
        surf%n(i)  = helicity*n
        surf%rc(i) = rc
        surf%zs(i) = zs
        if (.not. symm) then
            surf%rs(i) = rs
            surf%zc(i) = zc
        end if

    end do

    close(unit=unit_surf)

end subroutine read_surface_file

function to_lowercase(str_in) result(str_out)

    implicit none

    character(*), intent(IN) :: str_in
    character(len(str_in)) :: str_out
    character(len=26), parameter :: uppers = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=26), parameter :: lowers = 'abcdefghijklmnopqrstuvwxyz'
    integer :: i, n

    do i = 1, len(str_in)
        n = index(uppers, str_in(i:i))
        if (n > 0) then
            str_out(i:i) = lowers(n:n)
        else
            str_out(i:i) = str_in(i:i)
        end if
    end do

end function to_lowercase

end module input

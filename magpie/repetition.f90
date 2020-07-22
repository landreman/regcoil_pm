module repetition

use magpie_globals, only: dp

implicit none

contains

!-------------------------------------------------------------------------------
! geometry_repetition(repeat_to_fill, nfp, stell_symm, tor_symm, 
!                     nQhexes,     all_qhexes,     overlaps,
!                     nQhexes_rep, all_qhexes_rep, overlaps_rep)
!
! Expands arrays of magnet geometric properties to fill a selected extent
! around a torus. Presently supports only arrays of quadrilaterally-faced
! hexahedra.
!
! Input parameters:
!     character(len=*) :: repeat_to_fill 
!         -> portion of the torus to be filled through repetition:
!                'none':         no repetition; outputs will match inputs
!                'half-period':  only available for stellarator-symmetric 
!                                configurations; behavior is the same as 'none'
!                'field-period': only available for stellarator- or toroidally
!                                symmetric configurations
!                'half-torus':   only available for stellarator- or toroidally
!                                symmetric configurations. For toroidally 
!                                symmetric configs, will output the greatest
!                                number of field periods that forms less than or
!                                equal to half the torus.
!                'torus':        available for all symmetry cases
!     integer :: nfp -> number of field periods
!     logical :: stell_symm   -> true if the configuration is stellarator
!                                symmetric; input geometry is assumed to 
!                                characterize one half-period
!     logical :: tor_symm     -> true if the configuration is toroidally 
!                                symmetric. If stell_symm == .false., the
!                                input geometry is assumed to characterize one
!                                field period.
!     integer :: nQhexes      -> number of qhex units in the input arrays
!     type(quad_hexahedron), dimension(:) all_qhexes 
!                             -> array of input qhex structures 
!     logical, dimension(:) overlaps
!                             -> indicates which of the qhexes in all_qhexes
!                                overlaps with ports or other objects
!
! Output parameters:
!     integer :: nQhexes_rep  -> total number of qhex units in the expanded set
!     type(quad_hexahedron), allocatable, dimension(:) :: all_qhexes_rep
!                             -> expanded set of qhexes
!     logical, allocatable, dimension(:) :: overlaps_rep
!                             -> overlap array for the expanded set
!-------------------------------------------------------------------------------
subroutine geometry_repetition_qhex_ovl( &
               repeat_to_fill, nfp, stell_symm, tor_symm,  &
               nQhexes,     all_qhexes,     overlaps, &
               nQhexes_rep, all_qhexes_rep, overlaps_rep)

    use magpie_globals,  only: pi
    use qhex_properties, only: quad_hexahedron
    use qhex_adjust,     only: stell_qhex_transform
    use input,           only: to_lowercase

    implicit none

    character(len=*), intent(IN) :: repeat_to_fill
    integer, intent(IN) :: nfp
    logical, intent(IN) :: stell_symm, tor_symm
    integer, intent(IN) :: nQhexes
    type(quad_hexahedron), dimension(:), intent(IN) :: all_qhexes
    logical, dimension(:), intent(IN) :: overlaps
    integer, intent(OUT) :: nQhexes_rep
    type(quad_hexahedron), dimension(:), allocatable, &
                                        intent(OUT) :: all_qhexes_rep
    logical, dimension(:), allocatable, intent(OUT) :: overlaps_rep
    integer :: i, j

    ! Handle requests for magnet repetition in the corners file
    select case (to_lowercase(repeat_to_fill))

        case ('half-period')
            if (.not. stell_symm) then
                write(*,*) 'half-period geometric repetition is only ' // &
                           'available if stellarator symmetry is assumed'
                stop
            end if
            nQhexes_rep = nQhexes
            allocate(all_qhexes_rep(nQhexes_rep), &
                     overlaps_rep(nQhexes_rep)     )
            all_qhexes_rep   = all_qhexes
            overlaps_rep     = overlaps

        case ('field-period')
            if (.not. (stell_symm .or. tor_symm)) then
                write(*,*) 'full-period geometric repetition is only ' // &
                           'available if toroidal or stellarator symmetry ' // &
                           'is assumed'
                stop
            end if
            if (stell_symm) then
                nQhexes_rep = 2*nQhexes
                allocate(all_qhexes_rep(nQhexes_rep), &
                         overlaps_rep(nQhexes_rep)     )
                all_qhexes_rep(1:nQhexes) = all_qhexes
                overlaps_rep(1:nQhexes)    = overlaps
                do i = 1, nQhexes
                    call stell_qhex_transform('reflect', pi/nfp, &
                             all_qhexes(i), all_qhexes_rep(i+nQhexes))
                end do
                overlaps_rep(nQhexes+1:nQhexes_rep) = overlaps
            else 
                nQhexes_rep = nQhexes
                allocate(all_qhexes_rep(nQhexes_rep), &
                         overlaps_rep(nQhexes_rep)     )
                all_qhexes_rep   = all_qhexes
                overlaps_rep     = overlaps
            end if

        case ('half-torus') 

            if (stell_symm) then
                nQhexes_rep = nfp*nQhexes
                allocate(all_qhexes_rep(nQhexes_rep), &
                         overlaps_rep(nQhexes_rep)     )
                ! Iterate over half periods up to i = nfp (half the torus)
                do i = 1, nfp
                    if (mod(i,2) == 1) then
                        do j = 1, nQhexes
                            call stell_qhex_transform('translate', &
                                     (2*pi*(i/2))/nfp,             &
                                     all_qhexes(j),                &
                                     all_qhexes_rep((i-1)*nQhexes+j))
                        end do
                    else
                        do j = 1, nQhexes
                            call stell_qhex_transform('reflect',      &
                                     (i-1)*pi/nfp,                    &
                                     all_qhexes_rep((i-2)*nQhexes+j), &
                                     all_qhexes_rep((i-1)*nQhexes+j)) 
                        end do
                    end if
                    overlaps_rep((i-1)*nQhexes+1:i*nQhexes) = overlaps
                end do

            ! Note: will generate magnets to fill the largest number of 
            ! field periods that covers less than or equal to half the torus
            else if (tor_symm) then
                nQhexes_rep = (nfp/2)*nQhexes
                allocate(all_qhexes_rep(nQhexes_rep), &
                         overlaps_rep(nQhexes_rep)     )
                do i = 1, nfp/2
                    do j = 1, nQhexes
                        call stell_qhex_transform('translate', &
                                 (2*pi*(i-1))/nfp,            &
                                 all_qhexes(i),               &
                                 all_qhexes_rep((i-1)*nQhexes+j))
                    end do
                    overlaps_rep((i-1)*nQhexes+1:i*nQhexes) = overlaps
                end do

            else
                write(*,*) 'half-torus geometric repetition is only ' // &
                           'available if toroidal or stellarator symmetry ' // &
                           'is assumed'
                stop
            end if

        case ('torus') 

            if (stell_symm) then
                nQhexes_rep = 2*nfp*nQhexes
                allocate(all_qhexes_rep(nQhexes_rep), &
                         overlaps_rep(nQhexes_rep)     )
                do i = 1, nfp
                    do j = 1, nQhexes
                        call stell_qhex_transform('translate', &
                                 (2*pi*(i-1))/nfp,            &
                                 all_qhexes(j),               &
                                 all_qhexes_rep((2*i-2)*nQhexes+j))
                        call stell_qhex_transform('reflect',   &
                                 (2*(i-1)+1)*pi/nfp,          &
                                 all_qhexes_rep((2*i-2)*nQhexes+j), &
                                 all_qhexes_rep((2*i-1)*nQhexes+j)) 
                    end do
                    overlaps_rep((2*i-2)*nQhexes+1:(2*i-1)*nQhexes) &
                        = overlaps
                    overlaps_rep((2*i-1)*nQhexes+1:2*i*nQhexes) &
                        = overlaps
                end do

            else if (tor_symm) then
                nQhexes_rep = nfp*nQhexes
                allocate(all_qhexes_rep(nQhexes_rep), &
                         overlaps_rep(nQhexes_rep)     )
                do i = 1, nfp
                    do j = 1, nQhexes
                        call stell_qhex_transform('translate', &
                                 (2*pi*(i-1))/nfp,            &
                                 all_qhexes(i),               &
                                 all_qhexes_rep((i-1)*nQhexes+j))
                    end do
                    overlaps_rep((i-1)*nQhexes+1:i*nQhexes) = overlaps
                end do

            else
                nQhexes_rep = nQhexes
                allocate(all_qhexes_rep(nQhexes_rep), &
                         overlaps_rep(nQhexes_rep)     )
                all_qhexes_rep = all_qhexes
                overlaps_rep = overlaps
            end if

        case ('none')
            nQhexes_rep = nQhexes
            allocate(all_qhexes_rep(nQhexes_rep), &
                     overlaps_rep(nQhexes_rep)     )
            all_qhexes_rep = all_qhexes
            overlaps_rep = overlaps

        case default
            write(*,*) 'Unrecognized value for repeat_to_fill'
            stop

    end select

end subroutine geometry_repetition_qhex_ovl

subroutine geometry_repetition_qhex( &
               repeat_to_fill, nfp, stell_symm, tor_symm,  &
               nQhexes, all_qhexes, nQhexes_rep, all_qhexes_rep )

    use qhex_properties, only: quad_hexahedron
    use qhex_adjust,     only: stell_qhex_transform
    use input,           only: to_lowercase

    implicit none

    character(len=*), intent(IN) :: repeat_to_fill
    integer, intent(IN) :: nfp
    logical, intent(IN) :: stell_symm, tor_symm
    integer, intent(IN) :: nQhexes
    type(quad_hexahedron), dimension(:), intent(IN) :: all_qhexes
    type(quad_hexahedron), dimension(:), allocatable, &
                                        intent(OUT) :: all_qhexes_rep
    logical, dimension(:), allocatable :: overlaps, overlaps_rep
    integer, intent(OUT) :: nQhexes_rep

    allocate(overlaps(nQhexes))

    call geometry_repetition_qhex_ovl( &
               repeat_to_fill, nfp, stell_symm, tor_symm,  &
               nQhexes,     all_qhexes,     overlaps, &
               nQhexes_rep, all_qhexes_rep, overlaps_rep)

    deallocate(overlaps, overlaps_rep)

end subroutine geometry_repetition_qhex

end module repetition


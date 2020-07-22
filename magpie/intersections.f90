!-------------------------------------------------------------------------------
! intersections.f90
!
! Module containing top-level subroutines for managing which intersections to
! check for query points in the magnet array
!-------------------------------------------------------------------------------
module intersections

use magpie_globals, only: dp

implicit none

contains

!-------------------------------------------------------------------------------
! intersections_to_check
!
! For each family of possible intersections (i.e., coil sets for a device, port
! sets for a device), checks all of the input namelist variables corresponding
! to individual items within each family and sets the values of higher-level
! logical variables to indicate which modules to invoke when checking for
! intersection
!-------------------------------------------------------------------------------
subroutine intersections_to_check

    use magpie_globals, &
        only: incl_port_ncsx_nb, incl_port_ncsx_02, incl_port_ncsx_04, &
              incl_port_ncsx_05, incl_port_ncsx_06, incl_port_ncsx_07, &
              incl_port_ncsx_08, incl_port_ncsx_09, incl_port_ncsx_10, &
              incl_port_ncsx_11, incl_port_ncsx_12, incl_port_ncsx_15, &
              incl_port_ncsx_17, incl_port_ncsx_18, incl_port_ncsx_dome, &
              check_ncsx_ports

    implicit none

    if (incl_port_ncsx_nb .or. incl_port_ncsx_02 .or. incl_port_ncsx_04 .or. &
        incl_port_ncsx_05 .or. incl_port_ncsx_06 .or. incl_port_ncsx_07 .or. &
        incl_port_ncsx_08 .or. incl_port_ncsx_09 .or. incl_port_ncsx_10 .or. &
        incl_port_ncsx_11 .or. incl_port_ncsx_12 .or. incl_port_ncsx_15 .or. &
        incl_port_ncsx_17 .or. incl_port_ncsx_18 .or. incl_port_ncsx_dome    ) &
            check_ncsx_ports = .true.

end subroutine intersections_to_check

!-------------------------------------------------------------------------------
! check_intersections(x, y, z, intersects)
!
! Checks for intersections with any applicable set of objects.
!
! Currently supports just NCSX ports but can in principle be expanded to
! include other object models.
!
! Input parameters:
!     real(dp) :: x, y, z    -> x, y, and z coordinate of the query point
!
! Output parameter:
!     logical :: intersects  -> True if the query point intersects any 
!                               applicable object
!-------------------------------------------------------------------------------
subroutine check_intersections(x, y, z, intersects)

    use magpie_globals, only: check_ncsx_ports
    use ports_ncsx_eval, only: in_ncsx_port

    implicit none

    real(dp), intent(IN) :: x, y, z
    logical, intent(OUT) :: intersects
    logical :: intersects_subset

    intersects = .false.

    if (check_ncsx_ports) then
        call in_ncsx_port(x, y, z, intersects_subset)
        if (intersects_subset) intersects = .true.
    end if

end subroutine check_intersections

end module intersections


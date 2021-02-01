
! ============================================================================
!  Module for topography data
! ============================================================================

module topo_module

    implicit none

    logical, private :: module_setup = .false.

    real(kind=8), allocatable :: xtopo(:), ztopo(:)
    integer :: mx_topo, ntopofiles, test_topography
    real(kind=8) :: topo_missing


contains

!   ============================================================
    subroutine read_topo_settings(file_name)

    implicit none

    character(len=*), intent(in), optional :: file_name

    ! Locals
    integer, parameter :: iunit = 7
    character (len=150) :: topofname

    if (present(file_name)) then
        call opendatafile(iunit, file_name)
    else
        call opendatafile(iunit, 'topo.data')
    endif

    ! Read in value to use in place of no_data_value in topofile
    read(iunit,*) topo_missing

    ! not used now:
    read(iunit,*) test_topography

    read(iunit,*) ntopofiles
    if (ntopofiles /= 1) then
        write(6,*) '*** Currently require exactly 1 topofile'
        stop
    endif

    read(iunit,*) topofname
    close(iunit)

    call read_topo_file(topofname)

    end subroutine read_topo_settings


!   ============================================================
    subroutine read_topo_file(fname)

    implicit none

    character(len=150), intent(in) :: fname
    integer :: i
    integer, parameter :: iunit = 7

    open(unit=iunit, file=trim(fname), status='unknown',form='formatted')

    read(iunit,*) mx_topo
    allocate(xtopo(mx_topo), ztopo(mx_topo))

    write(6,"('Reading ',i6,' topo values from ')") mx_topo
    write(6,*) '    ',trim(fname)
    
    do i=1,mx_topo
        read(iunit,*) xtopo(i), ztopo(i)
    enddo
    close(iunit)

    end subroutine read_topo_file

!   ============================================================
    subroutine topo_integrate(mx,mbc,maux,aux)

    use grid_module, only: xp_edge,z_edge,mx_edge
    implicit none

    integer, intent(in) :: mx, mbc, maux
    real(kind=8), intent(inout) ::  aux(maux,1-mbc:mx+mbc)

    integer :: i,itopo1,itopo2,ifull
    real(kind=8) :: x1,x2,topo_int,dx1,dx2,slope1

    itopo1 = 1
    do i=1,mx
        ! edges of computational cell:
        x1 = xp_edge(i)
        x2 = xp_edge(i+1)

        do while (xtopo(itopo1) <= x1)
            itopo1 = itopo1 + 1
        enddo
        ! now itopo1 points to first topo point in cell

        if (itopo1 == 1) then
            write(6,*) '*** topo does not cover domain'
            write(6,*) '*** x1, xtopo(1): ',x1,xtopo(1)
            stop
        endif

        if (xtopo(itopo1) >= x2) then
            ! xtopo(itopo1-1) to left of x1, so topo points bracket cell
            slope1 = (ztopo(itopo1) - ztopo(itopo1-1)) &
                   / (xtopo(itopo1) - xtopo(itopo1-1))
            dx1 = (x1 - xtopo(itopo1-1))
            dx2 = (x2 - xtopo(itopo1-1))
            aux(1,i) = ztopo(itopo1-1) + 0.5d0*(dx1+dx2)*slope1

        else
            itopo2 = itopo1
            do while (xtopo(itopo2) < x2)
                itopo2 = itopo2 + 1
            enddo
            ! now itopo2 points to first topo point to right of cell
    
    
            ! compute integral of pw linear topo function
            ! start with portion between x1 and xtopo(itopo1):
            dx1 = xtopo(itopo1) - x1
            slope1 = (ztopo(itopo1) - ztopo(itopo1-1)) &
                   / (xtopo(itopo1) - xtopo(itopo1-1))
            topo_int = dx1 * (ztopo(itopo1) - 0.5d0*dx1*slope1)
    
            ! add in contributions from full topo intervals between x1,x2:
            do ifull=itopo1+1,itopo2-1
                dx1 = xtopo(ifull) - xtopo(ifull-1)
                topo_int = topo_int + dx1*0.5d0*(ztopo(ifull-1)+ztopo(ifull))
            enddo
    
            ! finish with portion between xtopo(itopo2-1) and x2:
            dx1 = x2 - xtopo(itopo2-1)
            slope1 = (ztopo(itopo2) - ztopo(itopo2-1)) &
                   / (xtopo(itopo2) - xtopo(itopo2-1))
            topo_int = topo_int + dx1 * (ztopo(itopo2-1) + 0.5d0*dx1*slope1)
    
            ! divide by width of computational cell to get cell average topo value:
            aux(1,i) = topo_int / (x2 - x1)
        endif
    enddo

    end subroutine topo_integrate

end module topo_module


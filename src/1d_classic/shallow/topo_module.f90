
! ============================================================================
!  Module for topography data
! ============================================================================

module topo_module

    implicit none

    logical, private :: module_setup = .false.

    real(kind=8), allocatable :: xtopo(:), ztopo(:)
    integer :: mx_topo, ntopofiles, test_topography
    real(kind=8) :: topo_missing

    ! for dtopo:
    integer :: mt_dtopo, mx_dtopo
    real(kind=8), allocatable :: x_dtopo(:),t_dtopo(:),dz_dtopo(:,:)

    logical :: topo_finalized

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

    use grid_module, only: mbc,mx,zcell
    implicit none

    character(len=150), intent(in) :: fname
    integer :: i,ibc
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
    
    write(6,*) '+++ call cell_average...'
    ! compute zcell by cell-averaging pw linear function
    ! defined by xtopo,ztopo:
    call cell_average(mx_topo, xtopo, ztopo, zcell)
    
    ! extrapolate to ghost cells:
    do ibc=1,mbc
        zcell(1-ibc) = zcell(1)
        zcell(mx+ibc) = zcell(mx)
    enddo

    end subroutine read_topo_file


!   ============================================================
    subroutine cell_average(mx_pwlin, x_pwlin, z_pwlin, z_cell)

    ! take arrays x_pwlin and z_pwlin with mx_pwlin points that define
    ! a pw linear function, and return z_cell containing cell averages
    ! of this function, for the grid cells defined in grid_module.

    use grid_module, only: xp_edge,mx_edge
    implicit none

    integer, intent(in) :: mx_pwlin
    real(kind=8), dimension(1:mx_pwlin), intent(in) ::  x_pwlin, z_pwlin
    real(kind=8), intent(out) ::  z_cell(-1:mx_edge+1)

    integer :: i,i1,i2,ifull
    real(kind=8) :: x1,x2,cell_int,dx1,dx2,slope1

    i1 = 1
    do i=1,mx_edge-1
        ! edges of computational cell:
        x1 = xp_edge(i)
        x2 = xp_edge(i+1)

        do while (x_pwlin(i1) <= x1)
            i1 = i1 + 1
        enddo
        ! now i1 points to first topo point in cell

        if (i1 == 1) then
            write(6,*) '*** pwlin function does not cover domain'
            write(6,*) '*** x1, x_pwlin(1): ',x1,x_pwlin(1)
            stop
        endif

        if (x_pwlin(i1) >= x2) then
            ! x_pwlin(i1-1) to left of x1, so nodal points bracket cell
            slope1 = (z_pwlin(i1) - z_pwlin(i1-1)) &
                   / (x_pwlin(i1) - x_pwlin(i1-1))
            dx1 = (x1 - x_pwlin(i1-1))
            dx2 = (x2 - x_pwlin(i1-1))
            z_cell(i) = z_pwlin(i1-1) + 0.5d0*(dx1+dx2)*slope1
            !write(6,*) '+++ i, z_cell, z_pwlin: ',i,z_cell(i),z_pwlin(i)

        else
            i2 = i1
            do while (x_pwlin(i2) < x2)
                i2 = i2 + 1
            enddo
            ! now i2 points to first nodal point to right of cell
    
    
            ! compute integral of pw linear function
            ! start with portion between x1 and x_pwlin(i1):
            dx1 = x_pwlin(i1) - x1
            slope1 = (z_pwlin(i1) - z_pwlin(i1-1)) &
                   / (x_pwlin(i1) - x_pwlin(i1-1))
            cell_int = dx1 * (z_pwlin(i1) - 0.5d0*dx1*slope1)
    
            ! add in contributions from full intervals between x1,x2:
            do ifull=i1+1,i2-1
                dx1 = x_pwlin(ifull) - x_pwlin(ifull-1)
                cell_int = cell_int + dx1*0.5d0*(z_pwlin(ifull-1)+z_pwlin(ifull))
            enddo
    
            ! finish with portion between x_pwlin(i2-1) and x2:
            dx1 = x2 - x_pwlin(i2-1)
            slope1 = (z_pwlin(i2) - z_pwlin(i2-1)) &
                   / (x_pwlin(i2) - x_pwlin(i2-1))
            cell_int = cell_int + dx1 * (z_pwlin(i2-1) + 0.5d0*dx1*slope1)

            ! divide by width of computational cell to get cell average topo value:
            z_cell(i) = cell_int / (x2 - x1)
            !write(6,*) '+++ i, z_cell, z_pwlin: ',i,z_cell(i),z_pwlin(i)
        endif
    enddo

    end subroutine cell_average

!   ============================================================
    subroutine read_dtopo_settings(file_name)

    implicit none

    character(len=*), intent(in), optional :: file_name

    ! Locals
    integer, parameter :: iunit = 7
    character (len=150) :: dtopofname
    integer :: ndtopofiles, dtopotype

    if (present(file_name)) then
        call opendatafile(iunit, file_name)
    else
        call opendatafile(iunit, 'dtopo.data')
    endif


    read(iunit,*) ndtopofiles

    if (ndtopofiles == 0) then
        topo_finalized = .true.

    else if (ndtopofiles == 1) then
        read(iunit,*) dtopofname, dtopotype
        close(iunit)
        call read_dtopo_file(dtopofname,dtopotype)
        topo_finalized = .false.

    else if (ndtopofiles > 1) then
        write(6,*) '*** Currently require at most 1 dtopofile'
        stop

    endif

    end subroutine read_dtopo_settings



!   ============================================================
    subroutine read_dtopo_file(fname,dtopotype)

    use grid_module, only: dz_cell

    implicit none

    integer, intent(in) :: dtopotype
    character(len=150), intent(in) :: fname
    integer :: i, k, status, dtopo_size
    integer, parameter :: iunit = 7
    real(kind=8) :: t0,tf,t,x,x1,x2,dt_dtopo


    if (dtopotype == 1) then
        ! determine how many x values at each time:
        open(unit=iunit, file=trim(fname), status='unknown',form='formatted')
        read(iunit, *) t0,x1
        dtopo_size = 1
        mx_dtopo = 1
        t = t0
        status = 0
        do while (status == 0)
            read(iunit,fmt=*,iostat=status) t,x
            if (t == t0 ) then
                mx_dtopo = mx_dtopo + 1
            endif
            dtopo_size = dtopo_size + 1
        end do
        write(6,*) '+++ mx_dtopo = ',mx_dtopo
        write(6,*) '+++ dtopo_size = ',dtopo_size
        mt_dtopo = (dtopo_size-1)/mx_dtopo
        write(6,*) '+++ mt_dtopo = ',mt_dtopo
        x2 = x
        tf = t
        dt_dtopo = (tf-t0)/(mt_dtopo-1)
        write(6,*) '+++ dt_dtopo = ',dt_dtopo
        close(iunit)
        

        allocate(x_dtopo(mx_dtopo), dz_dtopo(mt_dtopo,mx_dtopo), &
                 t_dtopo(mt_dtopo))

        write(6,"('Reading ',i6,' dtopo values from ')") mx_dtopo
        write(6,*) '    ',trim(fname)
        
        do k=1,mt_dtopo
            do i=1,mx_dtopo
                read(iunit,*) t_dtopo(k), x_dtopo(i), dz_dtopo(k,i)
            enddo
            ! x_dtopo, dz_dtopo defines a pw linear function,
            ! compute cell averages based on computational cells:
            call cell_average(mx_dtopo, x_dtopo, dz_dtopo, dz_cell)
        enddo
        close(iunit)

    else
        write(6,*) 'Unrecognized dtopotype = ',dtopotype
        stop
    endif


    end subroutine read_dtopo_file

!   ============================================================
    subroutine topo_update(t)

    !integer :: mt_dtopo, mx_dtopo
    !real(kind=8), allocatable :: x_dtopo(:),t_dtopo(:),dz_dtopo(:,:)

    ! assumes dz_cell(1:mt_dtopo, 1-mbc:mx_grid+mbc) has been set
    ! to cell average of dz_dtopo on each grid cell at each dtopo time.

    ! updates aux(1,:) to be (original zcell) + dz_cell at 
    ! time interpolated from dtopo times.

    use grid_module, only: zcell

    implicit none
    
    real(kind=8), intent(in) :: t

    integer :: i,k,kd1,kd2
    real(kind=8) :: tau,td1,td2

    if ((mt_dtopo == 0) .or. topo_finalized) return

    if (mt_dtopo == 1) then
        kd1 = 1
        kd2 = 1
        td1 = t
        td2 = t
        tau = 0.d0
    else
        do kd1=1,mt_dtopo
            if (t >= t_dtopo(kd1)) exit
        enddo
        if (kd1 == mt_dtopo) then
            topo_finalized = .true.
            return
        endif

        kd2 = kd1+1
        tau = (t - t_dtopo(kd1)) / (t_dtopo(kd2) - t_dtopo(kd1))
    endif

    !do i=1,mx_  FIX!!

    end subroutine topo_update

end module topo_module

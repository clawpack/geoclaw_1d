module grid_module

    implicit none
    save

    integer :: mx, mbc, grid_type
    real(kind=8) :: xlower, xupper
    logical :: mapped_grid

    real(kind=8), allocatable, dimension(:) :: xcell, xc_edge, xp_edge
    real(kind=8), allocatable, dimension(:) :: zcell, z_edge

    integer :: mx_edge

    ! to keep track of max depth, speed over all time:
    real(kind=8), allocatable, dimension(:) ::  hmax, smax

    real(kind=8) :: total_mass

contains

subroutine set_grid(mx,dx)

    use geoclaw_module, only: earth_radius, coordinate_system
    implicit none

    integer, intent(in) :: mx
    real(kind=8), intent(in) :: dx

    character(len=9) :: fname_grid
    character(len=150) :: fname_celledges
    integer :: i,j,mb
    real(kind=8) :: rim,rip,ric,c0i,cmi,cpi,r

    fname_grid = 'grid.data'
    call opendatafile(58,fname_grid)
    !read(58,*) mapped_grid
    read(58,*) grid_type
        ! grid_type == 0: uniform grid
        ! grid_type == 1: mapc2p.f90 used for mapping
        ! grid_type == 2: celledges.txt file gives edges
    close(unit=58)

    write(6,*) '+++ grid_type = ',grid_type

    if (grid_type == 2) then
    
        fname_celledges = 'celledges.txt'
        write(6,*) 'Reading cell edges from ',trim(fname_celledges)
        open(unit=58, file=fname_celledges, status='old',form='formatted')

        !read(58,*) radial  ! OLD, now using coordinate_system == -1
        !read(58,*) uniform_grid   ! OLD
    
        read(58,*) mx_edge

        if (mx_edge /= mx+1) then
            write(6,*) '*** expect mx_edge in grid.txt to equal mx+1'
            write(6,*) '*** mx_edge, mx:  ', mx_edge,mx
            stop
        endif 

    
        allocate(xc_edge(1-mbc:mx_edge+mbc))
        allocate(xp_edge(1-mbc:mx_edge+mbc))
        allocate(z_edge(1-mbc:mx_edge+mbc))
        allocate(xcell(1-mbc:mx+mbc))
        allocate(zcell(1-mbc:mx+mbc))
    
        write(6,*) 'Reading grid with mx = ',mx
        ! read in xc values and corresponding xp and z values
        do i=1,mx_edge
            read(58,*) xc_edge(i),xp_edge(i),z_edge(i)
        enddo
    else
        ! uniform grid with mx cells and mx+1 edges:
        if (dx /= (xupper - xlower) / mx) then
            write(6,*) '*** dx has unexpected value'
            stop
        endif

        mx_edge = mx+1
        do i=1,mx_edge
            xp_edge(i) = xlower + (i-1)*dx
            z_edge(i) = -4000.d0   ! NEED TO FIX FOR GENERAL TOPO
        enddo
    endif

    if ((coordinate_system == -1) .and. (xp_edge(1) < 0.d0)) then
        write(6,*) 'coordinate_system == -1, radial requires x >= 0'
        stop
    endif


    ! extend to ghost cells, to have same width as first interior cell:

    do mb=1,mbc
        xp_edge(1-mb) = 2.d0*xp_edge(2-mb) - xp_edge(3-mb)
        xp_edge(mx_edge+mb) = 2.d0*xp_edge(mx_edge+mb-1) - xp_edge(mx_edge+mb-2)
        z_edge(1-mb) = z_edge(1)
        z_edge(mx_edge+mb) = z_edge(mx_edge)
    enddo


    do i=1-mbc,mx+mbc
        ! cell centers based on grid edges:
        xcell(i) = 0.5d0*(xp_edge(i) + xp_edge(i+1))
        zcell(i) = 0.5d0*(z_edge(i) + z_edge(i+1))
    enddo


    ! for keeping track of max depth, speed over all time:
    ! initialize here and update in b4step1
    allocate(hmax(mx), smax(mx))
    hmax(:) = 0.d0
    smax(:) = 0.d0
    
end subroutine set_grid

end module grid_module

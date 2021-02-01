module grid_module

    implicit none
    save

    integer, parameter :: mx_grid_max=100000
    integer :: mx_grid
    logical :: radial, uniform_grid
    real(kind=8), dimension(0:mx_grid_max+1) :: dxm, dxc, cm, cp, c0

    ! assume mbc==2 ghost cells:
    real(kind=8), dimension(-1:mx_grid_max+3) ::  xgrid, zgrid, xcell, zcell

    ! to keep track of max depth, speed over all time:
    real(kind=8), allocatable, dimension(:) ::  hmax, smax


contains

subroutine set_grid(mx,dx)

    implicit none

    integer, intent(in) :: mx
    real(kind=8), intent(in) :: dx

    character(len=9) :: fname_grid
    integer :: i,j
    real(kind=8) :: rim,rip,ric,c0i,cmi,cpi,r

    
    fname_grid = 'grid.data'

    open(unit=58, file=fname_grid, status='old',form='formatted')
    read(58,*) radial
    read(58,*) uniform_grid

    if (uniform_grid) then
        write(6,*) '*** uniform grid not implemented yet'
        stop        
    else 
        read(58,*) mx_grid
        write(6,*) 'Reading grid with mx_grid = ',mx_grid
        if (mx_grid+1 > mx_grid_max) then
            write(6,*) '*** too many topo values'
            stop
        endif

        ! read in grid cell edges and topo value at each edge:
        do i=1,mx_grid+1
            read(58,*) xgrid(i),zgrid(i)
        enddo
    endif


    ! extend to ghost cells, to have same width as first interior cell:
    xgrid(0) = 2.d0*xgrid(1) - xgrid(2)
    xgrid(-1) = 2.d0*xgrid(0) - xgrid(1)
    xgrid(mx_grid+2) = 2.d0*xgrid(mx_grid+1) - xgrid(mx_grid)
    xgrid(mx_grid+3) = 2.d0*xgrid(mx_grid+2) - xgrid(mx_grid+1)
    zgrid(0) = zgrid(1)
    zgrid(-1) = zgrid(1)
    zgrid(mx_grid+2) = zgrid(mx_grid+1)
    zgrid(mx_grid+3) = zgrid(mx_grid+1)

    do i=-1,mx_grid+2
        ! cell centers based on grid edges:
        xcell(i) = 0.5d0*(xgrid(i) + xgrid(i+1))
        zcell(i) = 0.5d0*(zgrid(i) + zgrid(i+1))
    enddo

    do i=0,mx_grid+2
        dxm(i) = xcell(i) - xcell(i-1)
    enddo

    write(6,*) '+++ radial = ',radial

    do i=0,mx_grid+1
    
        ! define coefficients cm(i), cp(i), c0(i) so that 
        !   q_{xx} \approx cm(i)*Q(i-1) + c0(i)*Q(i) + cp(i)*Q(i+1)
        ! or in radial case,
        !   (1/r (r q_x))_x \approx cm(i)*Q(i-1) + c0(i)*Q(i) + cp(i)*Q(i+1)
        ! Note that approximation is not centered at xcell(i), so only 1st order
        ! in general, but should be 2nd order if grid is smoothly varying (?).

        if (uniform_grid) then
            dxc(i) = 2*dx
            cm(i) = 1.d0/dx**2
            c0(i) = -2.d0/dx**2
            cp(i) = 1.d0/dx**2
        else
            dxc(i) = (xcell(i+1) - xcell(i-1))  ! = 2*dx for uniform grid!
            cm(i) = 2.d0 / (dxm(i)*dxc(i))
            cp(i) = 2.d0 / (dxm(i+1)*dxc(i))
            c0(i) = -(cm(i) + cp(i))
        endif

        if (radial) then
            ! include factors for ((1/r)*(r*q)_r)_r 
            ! using form q_{rr} + (1/r)q_r - (1/r**2)q
            r = xcell(i) - xgrid(1)  ! assuming radial about left edge!
            cm(i) = cm(i) - 1.d0/(r*dxc(i))
            cp(i) = cp(i) + 1.d0/(r*dxc(i))
            c0(i) = c0(i) - 1.d0/(r**2)
        endif
    enddo

    ! for keeping track of max depth, speed over all time:
    ! initialize here and update in b4step1
    allocate(hmax(mx_grid), smax(mx_grid))
    hmax(:) = 0.d0
    smax(:) = 0.d0
    
end subroutine set_grid

end module grid_module

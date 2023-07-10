subroutine setaux(mbc,mx,xlower,dx,maux,aux)

    ! Called at start of computation before calling qinit, and
    ! when AMR is used, also called every time a new grid patch is created.
    ! Use to set auxiliary arrays aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc).
    ! Note that ghost cell values may need to be set if the aux arrays
    ! are used by the Riemann solver(s).
    !
    ! This version sets aux(1,:) to b(x) for shallow flow
    ! and aux(2,:) to the ratio of cell width to dxc (capacity function).

    !use geoclaw_module, only: dry_tolerance !uncomment if needed
    !use geoclaw_module, only: grav  !uncomment if needed
    use geoclaw_module, only: earth_radius, coordinate_system, DEG2RAD
    use grid_module, only: xp_edge,zcell,mx_edge
    use topo_module, only: topo_integrate

    implicit none
    integer, intent(in) :: mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(inout) ::  aux(maux,1-mbc:mx+mbc)

    !locals
    integer :: i,i0,i1,j
    real(kind=8) :: xcell,a

    if (mx+1 .ne. mx_edge) then
        write(6,*) 'mx_edge from grid.data must agree with mx+1'
        write(6,*) 'mx_edge = ',mx_edge
        write(6,*) 'mx = ',mx
        stop
        endif

    ! compute topo values and store in aux(1,:):
    !call topo_integrate(mx,mbc,maux,aux)
    

    do i=1,mx
        aux(1,i) = zcell(i)
        !write(6,*) '+++ i,zcell: ',i,zcell(i)
        aux(2,i) = (xp_edge(i+1) - xp_edge(i))/dx
        if (coordinate_system == 2) then
            ! convert degrees to meters:
            aux(2,i) = aux(2,i)* DEG2RAD * earth_radius
        endif
        if (aux(2,i) <= 0.d0) then
            write(6,*) '+++ i,xp_edge(i),xp_edge(i+1): ',i,xp_edge(i),xp_edge(i+1)
            endif
    enddo

    aux(:,0) = aux(:,1)
    aux(:,-1) = aux(:,1)
    aux(:,mx+1) = aux(:,mx)
    aux(:,mx+2) = aux(:,mx)

end subroutine setaux

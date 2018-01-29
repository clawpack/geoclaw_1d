subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.
    ! This default version simply sets eta = max(h + b,0)

    ! For more specific initial conditions
    !  copy this to an application directory and
    !  loop over all grid cells to set values of q(1:meqn, 1:mx).

    !use geoclaw_module, only: dry_tolerance !uncomment if needed
    use geoclaw_module, only: grav  !uncomment if needed
    use grid_module, only: xgrid,zgrid,mx_grid
    !use setprob_module, only: RC,RD,DC

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i
    real(kind=8) :: xcell,r,x0

    real(kind=8) :: eta, width
    real(kind=8) :: dz(mx)

    open(unit=33,file='qinit_okada.data',status='old',form='formatted')
    do i=1,mx
        read(33,*) dz(i)
        enddo

    do i=1,mx
      q(1,i) = max(0.d0, dz(i) - aux(1,i))
      q(2,i) = 0.d0

   enddo


end subroutine qinit

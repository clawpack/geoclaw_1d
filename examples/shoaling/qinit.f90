subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.
    ! This default version simply sets eta = max(h + b,0)

    ! For more specific initial conditions
    !  copy this to an application directory and
    !  loop over all grid cells to set values of q(1:meqn, 1:mx).

    !use geoclaw_module, only: dry_tolerance !uncomment if needed
    !use geoclaw_module, only: grav  !uncomment if needed
    use grid_module, only: xp_edge,z_edge,mx_edge

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i
    real(kind=8) :: xcell,r,x0

    real(kind=8) :: eta, width

    width = 2.d3    ! controls width of Gaussian
    x0 = 0.d3   ! initial location of Gaussian

    do i=1,mx
      xcell = 0.5*(xp_edge(i) + xp_edge(i+1))
      r = xcell  ! in meters, based on xlower=0
      eta = 2.d0 * exp(-((r-x0)/width)**2)
      q(1,i) = max(0.0, eta - aux(1,i))
      q(2,i) = 0.d0  !eta*sqrt(grav*q(1,i))  ! right-going

   enddo


end subroutine qinit

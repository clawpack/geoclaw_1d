subroutine b4step1(mbc,mx,meqn,q,xlower,dx,t,dt,maux,aux)

    ! Called before each call to step1.
    ! Use to set time-dependent aux arrays or perform other tasks.

    ! this version checks for negative depths 

    use geoclaw_module, only: dry_tolerance
    use grid_module, only: hmax, smax

    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc)

    !local variables
    integer :: i,ig,j,m,mvars
    real(kind=8) :: speed


      do i=1-mbc,mx+mbc
         if (q(1,i)<=dry_tolerance) then
            q(1,i) = max(q(1,i),0.0)
            do m=2,meqn
               q(m,i)=0.d0
            enddo
            speed = 0.d0
         else
            speed = abs(q(2,i)/q(1,i))
         endif
         
         if ((i>=1) .and. (i<=mx)) then
             ! keep track of max depth and speed:
             hmax(i) = dmax1(hmax(i), q(1,i))
             smax(i) = dmax1(smax(i), speed)
         endif

      enddo


end subroutine b4step1


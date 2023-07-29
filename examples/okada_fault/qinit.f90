subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.

    use geoclaw_module, only: sea_level
    use grid_module, only: xcell

    use geoclaw_module, only: dry_tolerance, grav

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i

    real(kind=8) :: eta
    real(kind=8) :: dz(mx)

    ! assume displaced surface agrees with sea floor deformation
    open(unit=33,file='dtopo_okada.data',status='old',form='formatted')

    write(6,*) '+++ mx = ',mx
    do i=1,mx
        read(33,*) dz(i)
        enddo

    do i=1,mx
        eta = dz(i)
        q(1,i) = max(0.d0, eta - aux(1,i))
        q(2,i) = 0.d0

   enddo


end subroutine qinit

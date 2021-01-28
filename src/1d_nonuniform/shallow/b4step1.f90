subroutine b4step1(mbc,mx,meqn,q,xlower,dx,t,dt,maux,aux)

    ! Called before each call to step1.
    ! Use to set time-dependent aux arrays or perform other tasks.

    ! this version checks for negative depths 

    use geoclaw_module, only: dry_tolerance
    use geoclaw_module, only: DEG2RAD, earth_radius, coordinate_system, pi
    use grid_module, only: hmax, smax, xp_edge, total_mass

    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc)

    !local variables
    integer :: i,ig,j,m,mvars
    real(kind=8) :: speed, cell_mass, x, capa_latitude, zeta

    total_mass = 0.d0 
    if (coordinate_system == 2) then
        capa_latitude = 2.d0*pi*earth_radius**2 ! to scale to surface area
    endif

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

             ! total mass
             !cell_mass = q(1,i) * (xp_edge(i+1) - xp_edge(i))

             ! total mass deviation based on zeta:
             if (aux(1,i) < 0.d0) then
                 zeta = q(1,i) + aux(1,i)  ! = eta = h+B
             else
                 zeta = q(1,i)             ! = h
             endif
             cell_mass = zeta * (xp_edge(i+1) - xp_edge(i))

             if (coordinate_system == 2) then
                 !cell_mass = cell_mass * DEG2RAD * earth_radius
                 !x = 0.5d0*(xp_edge(i) + xp_edge(i+1))
                 !cell_mass = cell_mass * cos(x*DEG2RAD) * earth_radius
                 ! for total mass on full sphere:
                 cell_mass = cell_mass * (sin(xp_edge(i+1)*DEG2RAD) - &
                             sin(xp_edge(i)*DEG2RAD))/dx * capa_latitude
             endif
             total_mass = total_mass + cell_mass
         endif

      enddo

    write(76,*) t,total_mass

end subroutine b4step1


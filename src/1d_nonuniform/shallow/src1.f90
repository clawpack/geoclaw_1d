
subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)

    ! Called to update q by solving source term equation
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version integrates manning friction or other friction terms if present
    
    ! Version from 1D branch of GeoClaw originally from Dave George

    ! Also handles radial source term, if desired.
    ! Note: assumes radial about lower boundary, wall bc should be imposed

    
    use geoclaw_module, only: dry_tolerance, grav, DEG2RAD
    use geoclaw_module, only: friction_forcing
    use geoclaw_module, only: frictioncoeff => friction_coefficient
    use geoclaw_module, only: earth_radius, coordinate_system
    use grid_module, only: mx_edge, xp_edge


    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)

    !Locals
    real(kind=8) :: gamma, rcell, u, xcell, tanxR
    integer :: i

      if (frictioncoeff.gt.0.d0 .and. friction_forcing) then
          ! integrate source term based on Manning formula
            do i=1,mx
               if (q(1,i)<=dry_tolerance) then
                  q(2,i) = 0.0
               else
                  gamma= dsqrt(q(2,i)**2)*(grav*frictioncoeff**2)/(q(1,i)**(7.0/3.0))
                  q(2,i)= q(2,i)/(1.d0 + dt*gamma)
              endif
            enddo
        endif

!      ----------------------------------------------------------------

    !if (radial) then  ! OLD
    if (coordinate_system == -1) then
        ! radial source term for SWE:
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                ! assume x is radial coordinate in meters, x>=0, 
                ! u = radial velocity
                rcell = 0.5*(xp_edge(i) + xp_edge(i+1))
                !rcell = 0.5*(xp_edge(i) + xp_edge(i+1)) - xp_edge(1) ! OLD
                q(1,i) = q(1,i) - dt/rcell * q(2,i)
                u = q(2,i)/q(1,i)
                q(2,i) = q(2,i) - dt/rcell * q(1,i)*u**2
            endif
         enddo
     endif

!      ----------------------------------------------------------------

    if (coordinate_system == 2) then
        ! source term for x = latitude in degrees -90 <= x <= 90,
        ! u = velocity in latitude direction (m/s) on sphere:
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                xcell = 0.5*(xp_edge(i) + xp_edge(i+1))
                tanxR = tan(xcell*DEG2RAD) / earth_radius
                q(1,i) = q(1,i) + dt * tanxR * q(2,i)
                u = q(2,i)/q(1,i)
                q(2,i) = q(2,i) + dt * tanxR * q(1,i)*u**2
            endif
         enddo
     endif



end subroutine src1



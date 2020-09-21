
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
    use grid_module, only: mx_grid, xgrid, radial, uniform_grid


    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)

    !Locals
    real(kind=8) :: gamma, rcell, u
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

    if (radial) then
        ! radial source term for SWE:
        do i=1,mx
            ! assume radial about left edge!
            rcell = 0.5*(xgrid(i) + xgrid(i+1)) - xgrid(1) 
            q(1,i) = q(1,i) - dt/rcell * q(2,i)
            if (q(1,i) .gt. dry_tolerance) then
                u = q(2,i)/q(1,i)
            else
                u = 0.d0
            endif
            q(2,i) = q(2,i) - dt/rcell * q(1,i)*u**2
         enddo
     endif



end subroutine src1



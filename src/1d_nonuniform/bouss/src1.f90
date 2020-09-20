
subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)

    use bouss_module, only: bouss,solve_tridiag
    use grid_module, only: mx_grid, xgrid, radial, uniform_grid
    use geoclaw_module, only: dry_tolerance, grav, DEG2RAD
    use geoclaw_module, only: friction_forcing
    use geoclaw_module, only: frictioncoeff => friction_coefficient

    implicit none
        
    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(inout) ::   q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(in) :: xlower, dx, t, dt

    real(kind=8) gamma, rcell, u

    real(kind=8)  q0(meqn,1-mbc:mx+mbc)
    real(kind=8) psi(mx+2)

    real(kind=8)   rk_stage(1:mx,4)
    real(kind=8) :: delt
    integer ::  i,k,ii,nstep,rk_order
    
    logical verbose


    verbose=.False.
       
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

!      ----------------------------------------------------------------

    if (bouss) then

        ! Boussinesq terms 

        nstep = 1  ! set > 1 to do subcycling with source term 
        delt = dt / nstep

        rk_order = 1  ! 1 for Forward Euler, 2 for second-order RK

        do ii=1,nstep

          rk_stage = 0.d0
          q0  = q

          !-----------------------
          if (rk_order == 1) then
              ! Forward Euler

              call solve_tridiag(mx,meqn,mbc,dx,q0,maux,aux,psi)
              
              ! Forward Euler update:
              ! Note component 1 of psi used for BC, so psi(2) updates q(:,1):
              q(2,1:mx) = q(2,1:mx) + delt*psi(2:mx+1)
              
          endif  ! rk_order == 1
          
          !-----------------------
          if (rk_order == 2) then
          
              ! Second-order explicit R-K
              ! First Stage

              call solve_tridiag(mx,meqn,mbc,dx,q0,maux,aux,psi)


              rk_stage(1:mx,1) = psi(2:mx+1)

              q0(2,1:mx)=q(2,1:mx) + delt/2.d0*rk_stage(1:mx,1)

              ! Second Stage

              call solve_tridiag(mx,meqn,mbc,dx,q0,maux,aux,psi)

              rk_stage(1:mx,2)=psi(2:mx+1)

              q(2,1:mx) = q(2,1:mx) + delt*rk_stage(1:mx,2)

          endif ! rk_order == 2

        enddo

    endif ! end of Bouss terms
  
end subroutine src1



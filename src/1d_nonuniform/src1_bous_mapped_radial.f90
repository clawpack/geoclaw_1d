
  subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)

      use bous_module
      use grid_module, only: mx_grid, xgrid, radial
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
      integer(kind=4) INFO
     
      INTEGER            LDB
      INTEGER            IPIV( mx+2 )
      DOUBLE PRECISION   D(mx+2), DL(mx+1), DU(mx+1), DU2(1:mx)             
      
      logical verbose
      real(kind=8) :: DLi


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
          ! radial source term:
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
!
!      Boussinesq terms 
    
     !nstep=ceiling(dt/dx*2*sqrt(2.))  ! ???
     !write(6,*) '+++ nstep = ',nstep
     nstep = 1

     k = 1

     delt=dt/nstep

     LDB = mx+2
     D  =0.d0
     DU =0.d0
     DL =0.d0
     DU2=0.d0

     rk_order = 1
    
     do ii=1,nstep

          call set_diag(mx,meqn,mbc,dx,q,maux,aux, DL, D, DU)
          
  661     format(3e16.6)
          if (.false.) then
              write(66,*) 'D matrix:'
              do i=1,mx+2
                  if (i>1) then
                      DLi = DL(i-1)
                  else
                      DLi = -9999.
                  endif
                  write(66,661) DLi, D(i), DU(i)
              enddo
          endif

          call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )
    
          psi = 0.d0
          rk_stage = 0.d0
          q0  = q

          !-----------------------
          if (rk_order == 1) then
          ! Forward Euler
      
          call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)

          if (.false.) then
              write(66,*) 'RHS psi:'
              do i=1,mx+2
                  write(66,661) psi(i)
              enddo
          endif
              

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )
                      
                      
          if (.false.) then
              write(66,*) 'Solution psi, dt = ', delt
              do i=1,mx+2
                  write(66,661) psi(i)
              enddo
          endif
          
          q(2,1:mx) = q(2,1:mx) - delt*psi(2:mx+1)
          endif
          
          !-----------------------
          if (rk_order == 2) then
          ! RK2   

          ! First Stage
      
          call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)


          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )


          rk_stage(1:mx,1) = psi(2:mx+1)
        
          q0(2,1:mx)=q(2,1:mx)-delt/2.d0*rk_stage(1:mx,1)
        
          ! Second Stage

          call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          rk_stage(1:mx,2)=psi(2:mx+1)
        
          q(2,1:mx) = q(2,1:mx)- delt*rk_stage(1:mx,2)

          endif ! rk_order==2

          !-----------------------
          if (rk_order == 4) then
          ! RK4   

          ! First Stage
      
          call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)


          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )


          rk_stage(1:mx,1) = psi(2:mx+1)
        
          q0(2,1:mx)=q(2,1:mx)-delt/2.d0*rk_stage(1:mx,1)
        
          ! Second Stage

          call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          rk_stage(1:mx,2)=psi(2:mx+1)
        
          q0(2,1:mx)=q(2,1:mx)-delt/2.d0*rk_stage(1:mx,2)
        
          ! Third Stage
        
          call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )
                
          rk_stage(1:mx,3)=psi(2:mx+1)
        
          q0(2,1:mx)=q(2,1:mx)-delt*rk_stage(1:mx,3)
        
          ! Fourth Stage
        
          call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          rk_stage(1:mx,4)=psi(2:mx+1)

                              
          q(2,1:mx) = q(2,1:mx)- delt/6.d0*(rk_stage(1:mx,1) &
               + 2.d0*rk_stage(1:mx,2) &
               + 2.d0*rk_stage(1:mx,3) + rk_stage(1:mx,4))

          endif ! rk_order==4

      enddo
      
    endif ! end of Bouss terms
      
  end subroutine src1

  
! =========================================================
      subroutine set_diag(mx,meqn,mbc,dx,q,maux,aux, DL, D, DU)
! =========================================================
      use bous_module, only: B_param, sw_depth0, bc_xlo, bc_xhi
      use geoclaw_module, only: sea_level
      use grid_module, only: xcell, dxm, dxc, cm, cp, c0

      implicit none
     
      integer(kind=4), intent(in) :: mx,meqn,mbc,maux
      real(kind=8), intent(in) ::  dx
      real(kind=8), intent(in) ::    q(meqn,1-mbc:mx+mbc)
      real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
      real(kind=8), intent(out) ::   D(mx+2), DL(mx+1), DU(mx+1)                
      
      real(kind=8), dimension(1-mbc:mx+mbc) :: h0, hh, hhh
      real(kind=8) :: tol
      integer(kind=4) :: i,k,mxy

      h0 = sea_level - aux(1,:)
      tol = 1d-4   ! for testing h0 and q(1).  Good value?
     
        do i=1-mbc,mx+mbc
           hh(i)=(max(0.,h0(i)))**2
           hhh(i)=(max(0.,h0(i)))**3
        enddo
      
        ! initialize to identity matrix:
        D = 1.d0
        DU= 0.d0
        DL= 0.d0
   
        do i=1,mx
            
            if ((h0(i) > sw_depth0) .and. &
                (minval(q(1,i-2:i+2)) >= tol) .and. &
                (minval(h0(i-2:i+2)) >= tol)) then
                
              ! replace this row of identity matrix with dispersion terms:

              D(i+1) = 1.d0 + c0(i)*((B_param+.5d0)*hh(i) &
                    - 1.d0/6.d0*hhh(i)/h0(i)) !+ 1.d0/xcell(i) ! new term
     
              DU(i+1)= -cp(i)*((B_param+.5d0)*hh(i) &
                    - 1.d0/6.d0*hhh(i)/h0(i+1))
     
              DL(i)= -cm(i)*((B_param+.5d0)*hh(i) &
                    - 1.d0/6.d0*hhh(i)/h0(i-1))
     
            endif
            
        enddo
        
      ! left boundary
      ! No change ==> Dirichlet value 0 in ghost cell
      

      if (bc_xlo==3) then
         ! for wall-reflecting BC at left: impose Q_0 = -Q_1
         DU(1) = 1.d0     
         endif
      if (bc_xlo==1) then
         ! for Neumann BC at left: impose Q_0 = Q_1
         DU(1) = -1.d0     
         endif
         
      ! other things tested:
      if (.false.) then
         ! for wall-reflecting BC at left: impose Q_0 = -Q_1
         D(2) = D(2) - DL(1)   
         DL(1) = 0.d0     
         endif
      if (.false.) then
         ! for zero Dirichlet at left: impose Q_0 = 0
         DL(1) = 0.d0     
         endif
      if (.false.) then
         ! for Neumann BC at left: impose Q_0 = Q_1
         D(2) = D(2) + DL(1)   
         DL(1) = 0.d0     
         endif

    
      ! right boundary
      ! No change ==> Dirichlet value 0 in ghost cell
            
      if (bc_xhi==1) then
         ! for Neumann at right:
         DL(mx+1) = -1.d0
         endif
      if (bc_xhi==3) then
         ! wall at right:
         DL(mx+1) = 1.d0
         endif

      if (.false.) then
         ! for Neumann at right:
         D(mx+1) = D(mx+1) + DU(mx+1)
         DU(mx+1) = 0.d0
         endif
                  
      return
      end subroutine set_diag
   
!======================================================================
      subroutine set_psi(mx,meqn,mbc,dx,q,maux,aux,psi)

      use geoclaw_module, only: g => grav, sea_level, dry_tolerance
      use bous_module, only: B_param, sw_depth0, sw_depth1
      use grid_module, only: xcell, dxm, dxc, cm, cp, c0
     
      implicit none

      integer, intent(in) :: mx,meqn,mbc,maux
      real(kind=8), intent(in) :: dx
     
      real(kind=8), intent(in) ::  q(meqn,1-mbc:mx+mbc)
      real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
      real(kind=8), intent(out) :: psi(mx+2)

      real(kind=8)  depth
      real(kind=8), dimension(1-mbc:mx+mbc) :: h0, hh, hhh, eta, hu2
      real(kind=8), dimension(1-mbc:mx+mbc) :: hetax,detax,s1,s1_h
      real(kind=8)  detaxxx,s1xx,s1_hxx
      integer :: i,j,k,kk,iL,iR
    
      h0 = sea_level - aux(1,:)
     
         do i=1-mbc,mx+mbc
           if (q(1,i).gt.dry_tolerance) then
              hu2(i)= q(2,i)**2/q(1,i)
             else
              hu2(i)= 0.d0
             endif
           !hu2(i) = 0.d0  ! Debug: turning off nonlinear terms
           eta(i)= q(1,i) - h0(i)
           hh(i)=(max(0.d0,h0(i)))**2
           hhh(i)=(max(0.d0,h0(i)))**3
        enddo
     
     hetax = 0.d0
     detax = 0.d0
     !write(6,*) '+++ eta(0:1): ',eta(0),eta(1)

         ! extend to one set of ghost cells:
         do i= 0,mx+1
         if (minval(h0(i-1:i+1)) > 0.d0) then
            if (i<0) then
                !write(6,*) '+++ dxm = ', (dxm(j), j=1,4)
                hetax(i)= q(1,i)*(eta(i+1)-eta(i))/dxm(i+1)
                detax(i)= h0(i)*(eta(i+1)-eta(i))/dxm(i+1)
            elseif (i>mx+1) then
                hetax(i)= q(1,i)*(eta(i)-eta(i-1))/dxm(i)
                detax(i)= h0(i)*(eta(i)-eta(i-1))/dxm(i)
            else
                hetax(i)= q(1,i)*(eta(i+1)-eta(i-1))/dxc(i)
                detax(i)= h0(i)*(eta(i+1)-eta(i-1))/dxc(i)
            endif
         endif
        enddo
     
     s1=0.d0
     s1_h=0.d0

     !write(6,*) '+++ hu2(0:1): ',hu2(0),hu2(1)
     !write(66,*) 'detax, h0, eta, s1, s1_h:'
     do i=0,mx+1
          if (i<0) then
             s1(i)= (hu2(i+1)-hu2(i))/dxm(i+1)
          elseif (i>mx+1) then
             s1(i)= (hu2(i)-hu2(i-1))/dxm(i)
          else
             !s1(i)= (hu2(i+1)-hu2(i-1))/dxc(i)  ! original, missing term
             s1(i)= (hu2(i+1)-hu2(i-1))/dxc(i) + hu2(i)/xcell(i)
             !s1(i)= (xcell(i+1)*hu2(i+1)-xcell(i-1)*hu2(i-1)) &
             !        / (xcell(i)*dxc(i))
          endif
         
          s1(i)= s1(i)+g*hetax(i)
      
           if (h0(i) > 0.1d0) then 
              s1_h(i)=s1(i)/h0(i)
           else
              s1_h(i)=0.d0
           endif
          !write(66,662) detax(i),h0(i),eta(i),s1(i),s1_h(i)
 662      format(5e16.6)
 663      format(i3,4e16.6)
     enddo
     
     !write(6,*) '+++ s1(0:1): ',s1(0),s1(1)
     !write(6,*) '+++ s1_h(0:1): ',s1_h(0),s1_h(1)
     !write(6,*) '+++ detax(0:1): ',detax(0),detax(1)
     !s1(0) = s1(1)
     !s1_h(0) = s1_h(1)
     !detax(0) = detax(1)
     
     !s1(mx+1) = s1(mx)
     !s1_h(mx+1) = s1_h(mx)
     !detax(mx+1) = detax(mx)

      psi = 0.d0

        do i=1,mx

           k = i+1

           if (h0(i) <= sw_depth0) then
              ! no dispersive term:
              psi(k) = 0.d0

           else
            
              ! check if we need to use one-sided differences:
              if (i<1) then
                 j = i+1
              elseif (i>mx) then
                 j = i-1
              else
                 j = i
              endif

              s1xx = cp(j)*s1(j+1) - c0(j)*s1(j) + cm(j)*s1(j-1)
              detaxxx = cp(j)*detax(j+1) - c0(j)*detax(j) + &
                        cm(j)*detax(j-1)
              s1_hxx = cp(j)*s1_h(j+1) - c0(j)*s1_h(j) + cm(j)*s1_h(j-1)
              
              !if ((i==1) .or. (i==mx)) then
              if (.false.) then
                 s1xx =  0.d0
                 detaxxx = 0.d0
                 s1_hxx = 0.d0
              endif

              psi(k) = (B_param+.5d0)*hh(i)*s1xx - B_param*g*hh(i)*detaxxx &
                  -hhh(i)/6.d0 * s1_hxx
                                   
           endif


        if ((h0(i) > sw_depth0) .and. (h0(i) < sw_depth1)) then
            ! reduce psi linearly between sw_depth0 and sw_depth1:
            psi(k) = psi(k) * (h0(i)-sw_depth0)/(sw_depth1-sw_depth0)
            endif
                
        enddo

      return
      end subroutine set_psi

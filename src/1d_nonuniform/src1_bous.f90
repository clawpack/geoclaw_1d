
  subroutine src1_bous(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)

      use bous_module

      implicit none
            
      integer, intent(in) :: meqn,mbc,mx,maux
      real(kind=8), intent(inout) ::   q(meqn,1-mbc:mx+mbc)
      real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
      real(kind=8), intent(in) :: xlower, dx, t, dt

      real(kind=8)  q0(meqn,1-mbc:mx+mbc)
      real(kind=8) psi(mx+2)
     
      real(kind=8)   rk_stage(1:mx,4)
      real(kind=8) :: delt
      integer ::  i,k,ii,nstep,rk_order
      real(kind=8) tol,x
      integer(kind=4) INFO
     
      INTEGER            LDB
      INTEGER            IPIV( mx+2 )
      DOUBLE PRECISION   D(mx+2), DL(mx+1), DU(mx+1), DU2(1:mx)             
      
      logical verbose

      verbose=.False.
           
!      ----------------------------------------------------------------
!
!      Boussinesq type-------------------------------------------------
    
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

     rk_order = 2
    
     do ii=1,nstep

          call set_diag(mx,meqn,mbc,dx,q,maux,aux, DL, D, DU)

          call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )
    
          psi = 0.d0
          rk_stage = 0.d0
          q0  = q

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
      
 999  continue
      
  end subroutine src1_bous
  
! =========================================================
      subroutine set_diag(mx,meqn,mbc,dx,q,maux,aux, DL, D, DU)
! =========================================================
      use bous_module
      use geoclaw_module, only: sea_level

      implicit none
     
      integer(kind=4) mx,meqn,mbc,maux,mxy
      integer(kind=4) i,k
      real(kind=8) dx

      real(kind=8)   q(meqn,1-mbc:mx+mbc)
      real(kind=8) aux(maux,1-mbc:mx+mbc)

      DOUBLE PRECISION   D(mx+2), DL(mx+1), DU(mx+1)                
      
      real(kind=8), dimension(1-mbc:mx+mbc) :: h0, hh, hhh
      
      h0 = sea_level - aux(1,:)
     
        do i=1-mbc,mx+mbc
           hh(i)=(max(0.,h0(i)))**2
           hhh(i)=(max(0.,h0(i)))**3
        enddo
      
        D = 1.d0
        DU= 0.d0
        DL= 0.d0

   
        do i=1,mx
        
            if ((minval(q(1,i-2:i+2)) < 0.1d0) &
                .or. (minval(h0(i-2:i+2)) < 0.1d0)) then

            !if ((q(1,i)<0.1d0) .or. (h0(i)<0.1d0) .or. &
            !    (minval(h0(i-2:i+2))<0.1d0) .or. &
            !    (maxval(q(1,i-2:i+2))<0.1d0) .or. &
            !    (minval(h0(i-1:i+1))<=0.1d0)) then

              ! do nothing
        
            else

              D(i+1) = 1.d0 + 2.d0*(B_param+.5d0)*hh(i)/dx**2 &
                    -2.d0/6.d0*hhh(i)/h0(i)/dx**2
     
              DU(i+1)=-(B_param+.5d0)*hh(i)/dx**2 &
                    +1.d0/6.d0*hhh(i)/h0(i+1)/dx**2
     
              DL(i)=-(B_param+.5d0)*hh(i)/dx**2 &
                    +1.d0/6.d0*hhh(i)/h0(i-1)/dx**2
     
            endif
            
        enddo
             
      return
      end subroutine set_diag
   
!======================================================================
      subroutine set_psi(mx,meqn,mbc,dx,q,maux,aux,psi)

      use geoclaw_module, only: g => grav, sea_level
      use bous_module
     
      implicit none

      integer, intent(in) :: mx,meqn,mbc,maux
      real(kind=8), intent(in) :: dx
     
      real(kind=8), intent(in) ::  q(meqn,1-mbc:mx+mbc)
      real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
      real(kind=8), intent(out) :: psi(mx+2)

      real(kind=8)  tol
      real(kind=8)  depth, depth0, depth1
      real(kind=8), dimension(1-mbc:mx+mbc) :: h0, hh, hhh, eta, hu2
      real(kind=8), dimension(1-mbc:mx+mbc) :: hetax,detax,s1,s1_h
      integer :: i,j,k,kk,iL,iR

      ! phase out dispersive term between depth0 and depth1:
      depth0 = 180.d0  ! pure SWE in shallower water than this
      depth1 = 190.d0  ! full Bouss in deeper water than this
    
      h0 = sea_level - aux(1,:)
     
         do i=1-mbc,mx+mbc
           if (q(1,i).gt.1d-4) then
              hu2(i)= q(2,i)**2/q(1,i)
             else
              hu2(i)= 0.d0
             endif
           eta(i)= q(1,i) - h0(i)
           hh(i)=(max(0.d0,h0(i)))**2
           hhh(i)=(max(0.d0,h0(i)))**3
        enddo
     
     hetax = 0.d0
     detax = 0.d0

         do i= 1,mx
         if (minval(h0(i-1:i+1)) > 0.d0) then
            if (i==1) then
                hetax(i)= q(1,i)*(eta(i+1)-eta(i))/dx
                detax(i)= h0(i)*(eta(i+1)-eta(i))/dx
            elseif (i==mx) then
                hetax(i)= q(1,i)*(eta(i)-eta(i-1))/dx
                detax(i)= h0(i)*(eta(i)-eta(i-1))/dx
            else
                hetax(i)= q(1,i)*(eta(i+1)-eta(i-1))/dx/2.
                detax(i)= h0(i)*(eta(i+1)-eta(i-1))/dx/2.
            endif
         endif
        enddo
     
     s1=0.d0
     s1_h=0.d0

     do i=1,mx
          if (i==1) then
             s1(i)= (hu2(i+1)-hu2(i))/dx
          elseif (i==mx) then
             s1(i)= (hu2(i)-hu2(i-1))/dx
          else
             s1(i)= (hu2(i+1)-hu2(i-1))/2.d0/dx
          endif
         
          s1(i)= s1(i)+g*hetax(i)
      
           if (h0(i) > 0.1d0) then 
              s1_h(i)=s1(i)/h0(i)
           else
              s1_h(i)=0.d0
           endif
     enddo
     
      tol = 1d-8


      psi=0.d0

        do i=1,mx

           k = i+1

           if ((h0(i) < depth0) .or. (q(1,i) < depth0)) then
              ! no dispersive term:
              psi(k) = 0.d0

           else
            
              if (i==1) then
                 
                     psi(k)=(B_param+.5d0)*hh(i)* &
                      (s1(i+2)-2.d0*s1(i+1)+s1(i))/dx**2. &
                     -B_param*g*hh(i)* &
                      (detax(i+2)-2.d0*detax(i+1)+detax(i))/dx**2 &
                     -hhh(i)/6.d0*(s1_h(i+2)-2.d0*s1_h(i+1) &
                     +s1_h(i))/dx**2
              
              elseif (i==mx) then

                     psi(k)=(B_param+.5d0)*hh(i)* &
                      (s1(i)-2.d0*s1(i-1)+s1(i-2))/dx**2. &
                     -B_param*g*hh(i)* &
                      (detax(i)-2.d0*detax(i-1)+detax(i-2))/dx**2 &
                     -hhh(i)/6.d0*(s1_h(i)-2.d0*s1_h(i-1) &
                     +s1_h(i-2))/dx**2
                          
              else
               
                  psi(k)=(B_param+.5d0)*hh(i)* &
                   (s1(i+1)-2.d0*s1(i)+s1(i-1))/dx**2. &
                  -B_param*g*hh(i)* &
                   (detax(i+1)-2.d0*detax(i)+detax(i-1))/dx**2 &
                  -hhh(i)/6.d0*(s1_h(i+1)-2.d0*s1_h(i) &
                  +s1_h(i-1))/dx**2
               
              endif
           endif


        if (h0(i) .lt. depth1) then
            ! reduce psi linearly between depth0 and depth1:
            psi(k) = psi(k) * (h0(i)-depth0)/(depth1-depth0)
            endif
                
        enddo

 994  continue

      return
      end subroutine set_psi

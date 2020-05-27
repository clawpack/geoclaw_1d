module bouss_module
    
    implicit none
    logical      :: bouss  ! Turn on the dispersive terms
    real(kind=8) :: B_param 
    real(kind=8) :: sw_depth0, sw_depth1
    integer :: bc_xlo, bc_xhi
    
    ! for tridiagonal solver:
    integer, allocatable, dimension(:) :: IPIV
    real(kind=8), allocatable, dimension(:) :: D, DL, DU, DU2    
    
    save


contains

!======================================================================

    subroutine set_bouss(mx,mbc,mthbc)
    
    ! Set things up for Boussinesq solver, in particular
    ! create and factor tridiagonal matrix for implicit solves.
    
    ! Boussinesq equation parameter B_param is now set in setprob    


    use geoclaw_module, only: sea_level
    use grid_module, only: xcell, zgrid, dxm, dxc, cm, cp, c0
    
    implicit none
    integer, intent(in) :: mx, mbc, mthbc(2)
    real(kind=8) :: DLi, zcell
    integer :: i
    real(kind=8), dimension(1-mbc:mx+mbc) :: h0, h02, h03
    integer(kind=4) :: INFO


    ! Boundary conditions to impose in computing Boussinesq update:
    
    if ((mthbc(1)==2) .or. (mthbc(2)==2)) then
        write(6,*) '*** Periodic BCs not supported in bouss_module'
        stop
    endif

    ! Use Dirichlet BC with value 0 for correction in ghost cells by default:
    bc_xlo = 0
    bc_xhi = 0
    
    ! Check if wall BCs specified (also be used for radial symmetry at r=0):
    
    if (mthbc(1)==3) then
        ! For wall boundary conditions at left boundary:
        bc_xlo = 3
    endif

    if (mthbc(2)==3) then
        ! For wall boundary conditions at right boundary:
        bc_xhi = 3
    endif

    ! To try out Neumann BCs:
    ! This doesn't seem to work well, so not a general option
    !bc_xlo = 1
    !bc_xhi = 1
    
    
    !------------------------------------
    ! Form tridiagonal matrix and factor
    ! Done once at start and then used in each time step
    
    !h0 = sea_level - aux(1,:)  # resting depth

    do i=1-mbc,mx+mbc
        zcell = 0.5d0*(zgrid(i) + zgrid(i+1))
        h0(i) = sea_level - zcell
        h02(i)=(max(0.,h0(i)))**2
        h03(i)=(max(0.,h0(i)))**3
    enddo

    ! allocate tridiagonal and initialize to identity matrix:
    allocate(D(mx+2), DL(mx+1), DU(mx+2), DU2(mx+2), IPIV(mx+2))
    D = 1.d0
    DU= 0.d0
    DL= 0.d0
    
    ! First and last rows (rows 1 and mx+2) of matrix corresponds to BCs
    ! Row i+1 of matrix corresponds to equation for cell i

    do i=1,mx
    
        ! Modify row i+1 (equation for i'th grid cell)
        ! unless this is a cell where SWE are used, then leave as row of I
        ! and set RHS to 0, so no modification to SWE result.
        
        if ((h0(i) > sw_depth0)) then
            
          ! Replace this row of identity matrix with dispersion terms
          ! Note that cm(i)*w(i-1) + c0(i)*w(i) + cp(i)*w(i+1) gives
          ! approximation to d^2 w = w_{xx} in plane wave case, or
          !                  d^2 w = (1/r * (r*w)_r)_r  in radial case
          ! Also h02 is h0**2 and h03 is h0**3:

          D(i+1) = 1.d0 - c0(i)*((B_param+.5d0)*h02(i) &
                - 1.d0/6.d0*h03(i)/h0(i))

          DU(i+1)= -cp(i)*((B_param+.5d0)*h02(i) &
                - 1.d0/6.d0*h03(i)/h0(i+1))

          DL(i)= -cm(i)*((B_param+.5d0)*h02(i) &
                - 1.d0/6.d0*h03(i)/h0(i-1))

        endif
        
    enddo

    ! left boundary
    ! No change ==> Dirichlet value 0 in ghost cell, Q_0 = 0
    ! For other BCs, change superdiagonal in row 1:

    if (bc_xlo==1) then
        ! for Neumann BC at left: impose Q_0 = Q_1
        DU(1) = -1.d0     
    endif
    if (bc_xlo==3) then
        ! for wall-reflecting BC at left: impose Q_0 = -Q_1
        DU(1) = 1.d0     
    endif


    ! right boundary
    ! No change ==> Dirichlet value 0 in ghost cell, Q_{mx+1} = 0
    ! For other BCs, change subdiagonal in row mx+2:

    if (bc_xhi==1) then
        ! for Neumann at right: impose Q_mx = Q_{mx+1}
        DL(mx+1) = -1.d0
    endif
    if (bc_xhi==3) then
        ! for wall at right: impose Q_mx = - Q_{mx+1}
        DL(mx+1) = 1.d0
    endif

    ! factor the tridiagonal matrix:    
    call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )    
    
    if (.true.) then
        ! DEBUG:
        write(66,*) 'D matrix:'
        do i=1,mx+2
            if (i>1) then
                DLi = DL(i-1)
            else
                DLi = -9999.  ! since no subdiag in first row
            endif
            write(66,661) DLi, D(i), DU(i)
  661       format(3e16.6)
        enddo
    endif


    end subroutine set_bouss
    
!======================================================================

    subroutine solve_tridiag(mx,meqn,mbc,dx,q,maux,aux,psi)

    ! Set up right hand side for tridagonal system and then solve it
    ! psi first holds RHS and then the solution to be passed back.
    
    ! Note that matrix has already been factored in set_bouss, so the
    ! arrays D, DL, DU, DU2, IPIV are already set up for the solver.
    
    use geoclaw_module, only: g => grav, sea_level, dry_tolerance
    use grid_module, only: xcell, dxm, dxc, cm, cp, c0, radial

    implicit none

    integer, intent(in) :: mx,meqn,mbc,maux
    real(kind=8), intent(in) :: dx

    real(kind=8), intent(in) ::  q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(out) :: psi(mx+2)

    real(kind=8)  depth
    real(kind=8), dimension(1-mbc:mx+mbc) :: h0, h02, h03, eta, hu2
    real(kind=8), dimension(1-mbc:mx+mbc) :: h_eta_x,h0_eta_x,s1,s1_h0
    real(kind=8)  h0_eta_xxx,s1_xx,s1_h0_xx
    integer :: i,j,k,kk,iL,iR
    integer(kind=4) :: INFO
    integer :: LDB

    h0 = sea_level - aux(1,:)

    do i=1-mbc,mx+mbc
        ! compute hu2 = h*u**2
        if (q(1,i).gt.dry_tolerance) then
            hu2(i)= q(2,i)**2/q(1,i)
        else
            hu2(i)= 0.d0
        endif
        
        !hu2(i) = 0.d0  ! Debug: turning off nonlinear terms
        
        ! h0 is resting depth, h02 = h0**2 and h03 = h0**3:
        eta(i)= q(1,i) - h0(i)
        h02(i)=(max(0.d0,h0(i)))**2
        h03(i)=(max(0.d0,h0(i)))**3
        
    enddo

    ! compute these quantities, extending to one set of ghost cells:
    ! h_eta_x  = h * eta_x
    ! h0_eta_x = h0 * eta_x
    
    ! defaults to 0 if h0==0 at one of the stencil points:
    h_eta_x = 0.d0
    h0_eta_x = 0.d0
        
    do i= 0,mx+1
        if (minval(h0(i-1:i+1)) > 0.d0) then
            ! if h0 > 0 at all stencil points.  Note dxc(i) = 2*dx if uniform:
            h_eta_x(i)= q(1,i)*(eta(i+1)-eta(i-1))/dxc(i)
            h0_eta_x(i)= h0(i)*(eta(i+1)-eta(i-1))/dxc(i)
        endif
    enddo
     

    
    do i=0,mx+1
    
        ! compute approximation to
        !   s1 = (h*u**2)_x + g*h*eta_x   (plane wave) or
        !   s1 = (1/r)*(r*h*u**2)_r + r*g*h*eta_r   (radial)
        !      = (h*u**2)_r + (1/r)*h*u**2 + r*g*h*eta_r
    
        s1(i) =  (hu2(i+1) - hu2(i-1)) / dxc(i) + g*h_eta_x(i)
                 
        if (radial) then
            ! add term for (1/r)*h*u**2:
            s1(i) =  s1(i) + hu2(i) / xcell(i) 
        endif
        
        ! compute s1_h0 = s1 / h0:
        
        if (h0(i) > 0.1d0) then 
            s1_h0(i)=s1(i)/h0(i)
        else
            s1_h0(i)=0.d0
        endif
    enddo
     

    ! Compute right hand side psi:
    
    psi = 0.d0

    do i=1,mx

       k = i+1  ! i'th equation corresponds to row k=i+1 of RHS

       if (h0(i) <= sw_depth0) then
            ! no dispersive term:
            psi(k) = 0.d0

        else
            ! Note that cm(i)*w(i-1) + c0(i)*w(i) + cp(i)*w(i+1) gives
            ! approximation to d^2 w = w_{xx} 
            ! or in radial case to d^2 w = (1/r * (rw)_r)_r

            ! compute s1_xx, approximates d^2 s1:

            s1_xx = cp(i)*s1(i+1) + c0(i)*s1(i) + cm(i)*s1(i-1)

            ! compute s1_h0_xx, approximates d^2 s1_h0 = d^2(s1 / h0):
            
            s1_h0_xx = cp(i)*s1_h0(i+1) + c0(i)*s1_h0(i) + cm(i)*s1_h0(i-1)
        
            ! compute h0_eta_xxx, approximates d^2 h0_eta_x:
            
            h0_eta_xxx = cp(i)*h0_eta_x(i+1) + c0(i)*h0_eta_x(i) + &
                    cm(i)*h0_eta_x(i-1)
                    

            ! Right-hand side:
            
            psi(k) = -(B_param+.5d0) * h02(i) * s1_xx &
                        + h03(i)/6.d0 * s1_h0_xx &
                        + B_param * g * h02(i) * h0_eta_xxx
              
                               
        endif


        if ((h0(i) > sw_depth0) .and. (h0(i) < sw_depth1)) then
            ! reduce psi linearly between sw_depth0 and sw_depth1:
            psi(k) = psi(k) * (h0(i)-sw_depth0)/(sw_depth1-sw_depth0)
        endif
            
    enddo
    
    if (.true.) then
        write(66,*) 'RHS:'
        do i=1,mx+2
            write(66,*) psi(i)
        enddo
    endif
    
    ! -------------------------
    ! solve tridiagonal system.  
    ! psi is RHS on entry, solution on exit:
    
    LDB = mx+2  ! leading dimension of right-hand side matrix psi
    
    call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2, &
                 IPIV, psi, LDB,INFO )
                                
    if (.true.) then
        write(66,*) 'Solution psi: '
        do i=1,mx+2
            write(66,*) psi(i)
        enddo
    endif

    return
end subroutine solve_tridiag


end module bouss_module

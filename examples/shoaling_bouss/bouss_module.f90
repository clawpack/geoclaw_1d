module bouss_module
    
    implicit none
    logical :: bouss  ! Turn on the dispersive terms
    logical :: update_bouss_matrix ! reform and factor each time step?
    real(kind=8) :: B_param 
    real(kind=8) :: sw_depth0, sw_depth1, bouss_wave_factor
    integer :: bc_xlo, bc_xhi


    ! for higher-order derivatives:
    real(kind=8), allocatable, dimension(:) :: dxm, dxc, cm, cp, c0

    logical, allocatable, dimension(:) :: id_row, bouss_cell
    
    ! for tridiagonal solver:
    integer, allocatable, dimension(:) :: IPIV
    real(kind=8), allocatable, dimension(:) :: D, DL, DU, DU2    

    save


contains

!======================================================================

    subroutine set_bouss(mx,mbc,mthbc,meqn)
    
    ! Set things up for Boussinesq solver, in particular
    ! create and factor tridiagonal matrix for implicit solves.

    use geoclaw_module, only: sea_level, coordinate_system
    use grid_module, only: xcell, zcell, grid_type
    
    implicit none
    integer, intent(in) :: mx, mbc, mthbc(2), meqn
    !real(kind=8), intent(in), optional :: q(meqn, 1-mbc,mx+mbc)
    integer :: i, iunit
    character*25 fname
    real(kind=8) :: r

    iunit = 7
    fname = 'bouss.data'
 !  # open the unit with new routine from Clawpack 4.4 to skip over
 !  # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(7,*) bouss
    read(7,*) B_param
    read(7,*) sw_depth0
    read(7,*) sw_depth1
    
    bouss_wave_factor = 0.8d0
    update_bouss_matrix = .false.
    
    allocate(dxm(0:mx+1),dxc(0:mx+1),cm(0:mx+1),cp(0:mx+1),c0(0:mx+1))
    allocate(id_row(0:mx+1), bouss_cell(1:mx))

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

    if (mthbc(1)==0) then
        ! For wavemaker BC at left
        bc_xlo = 3
    endif

    ! To try out Neumann BCs:
    ! This doesn't seem to work well, so not a general option
    !bc_xlo = 1
    !bc_xhi = 1
    
    
    !------------------------------------
    ! coefficients needed for second-order derivatives:

    do i=0,mx+2
        dxm(i) = xcell(i) - xcell(i-1)
    enddo

    do i=0,mx+1
        if (grid_type == 0) then
            ! uniform grid
            dxc(i) = 2*dxm(i)
            cm(i) = 1.d0/dxm(i)**2
            c0(i) = -2.d0/dxm(i)**2
            cp(i) = 1.d0/dxm(i)**2
        else
            dxc(i) = (xcell(i+1) - xcell(i-1))  ! = 2*dx for uniform grid!
            cm(i) = 2.d0 / (dxm(i)*dxc(i))
            cp(i) = 2.d0 / (dxm(i+1)*dxc(i))
            c0(i) = -(cm(i) + cp(i))
        endif

        if (coordinate_system == -1) then
            ! x = radial coordinate r >= 0
            ! include factors for ((1/r)*(r*q)_r)_r 
            ! using form q_{rr} + (1/r)q_r - (1/r**2)q
            r = xcell(i)
            cm(i) = cm(i) - 1.d0/(r*dxc(i))
            cp(i) = cp(i) + 1.d0/(r*dxc(i))
            c0(i) = c0(i) - 1.d0/(r**2)
        else if (coordinate_system == 2) then
            write(6,*) '*** latitude coordinates not yet implemented in Bouss'
            stop
        endif
    enddo
    
    ! initial setup and factorization of matrix:
    call factor_bouss_matrix(mx,mbc,meqn)
    
    end subroutine set_bouss
    

    subroutine factor_bouss_matrix(mx,mbc,meqn,q)
    
    ! Form tridiagonal matrix and factor.
    ! If optional argument q is present, set up matrix and factor based
    ! on current depth q(1,:), otherwise based only on initial h0(:).
    
    use geoclaw_module, only: sea_level, coordinate_system
    use grid_module, only: xcell, zcell, grid_type
    
    implicit none
    integer, intent(in) :: mx, mbc, meqn
    !integer, intent(in), optional :: meqn
    real(kind=8), intent(in), optional :: q(meqn, 1-mbc:mx+mbc)
    real(kind=8) :: DLi, r, eta
    integer :: i
    real(kind=8), dimension(1-mbc:mx+mbc) :: h0, h02, h03
    integer(kind=4) :: INFO
    
    !------------------------------------
    
    do i=1-mbc,mx+mbc
        h0(i) = sea_level - zcell(i)
  666   format(i3,3e16.6)
        h02(i)=(max(0.,h0(i)))**2
        h03(i)=(max(0.,h0(i)))**3
    enddo

    if (.not. present(q)) then
        ! allocate tridiagonal
        allocate(D(mx+2), DL(mx+1), DU(mx+2), DU2(mx+2), IPIV(mx+2))
    endif
    
    ! initialize to identity matrix:
    D = 1.d0
    DU= 0.d0
    DL= 0.d0
    
    ! First and last rows (rows 1 and mx+2) of matrix corresponds to BCs
    ! Row i+1 of matrix corresponds to equation for cell i

    id_row(1) = .false.
    id_row(mx+2) = .false.
    
    do i=1,mx
    
        ! Modify row i+1 (equation for i'th grid cell) to form (I-d^2) operator,
        ! unless this is a cell where SWE are used, then leave as row of I
        ! and set RHS to 0, so no modification to SWE result.
        
        if (present(q)) then
            ! if updating matrix each step, use current eta to determine:
            eta = q(1,i) - h0(i)    
            id_row(i+1) = (h0(i) <= sw_depth0) .or. &
                        (abs(eta) >= bouss_wave_factor*h0(i))
            
        else
            ! only have initial depth to work with:
            id_row(i+1) = (h0(i) <= sw_depth0)
        endif
        
        if (.not. id_row(i+1)) then
            
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

    
    if (.false.) then
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
        write(66,*) '===== end of D'
    endif


    ! factor the tridiagonal matrix:    
    call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )

    end subroutine factor_bouss_matrix
    
    
!======================================================================

    subroutine solve_tridiag(mx,meqn,mbc,dx,q,maux,aux,psi)

    ! Set up right hand side for tridagonal system and then solve it
    ! psi first holds RHS and then the solution to be passed back.
    
    ! Note that matrix has already been factored in set_bouss, so the
    ! arrays D, DL, DU, DU2, IPIV are already allocated, and are
    ! set up for the solver in the case update_bouss_matrix == .false.
    
    ! If update_bouss_matrix is .false., the initial matrix based on h0
    ! and its factorization computed in set_bouss is used.
    ! If update_bouss_matrix is .true., calculate and factor a new matrix
    ! each time step based on the current depth q(1,:).
    
    use geoclaw_module, only: g => grav, sea_level, dry_tolerance
    use geoclaw_module, only: coordinate_system
    use grid_module, only: xcell

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
    integer :: sw_row_count, reset_count
    
    if (update_bouss_matrix) then
        call factor_bouss_matrix(mx,mbc,meqn,q)
    endif

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
        !s1(i) =  g*h_eta_x(i)  ! debug
        !s1(i) =  (hu2(i+1) - hu2(i-1)) / dxc(i) 
                 
        if (coordinate_system == -1) then
            ! radial: add term for (1/r)*h*u**2:
            s1(i) =  s1(i) + hu2(i) / xcell(i) 
        else if (coordinate_system == 2) then
            ! latitude
            write(6,*) '*** latitude not yet implemented in bouss_module'
            stop
        endif
        
        ! compute s1_h0 = s1 / h0:
        
        if (h0(i) > dry_tolerance) then
            s1_h0(i)=s1(i)/h0(i)
        else
            s1_h0(i)=0.d0
        endif
    enddo
     

    ! Compute right hand side psi:
    
    psi = 0.d0

    do i=1,mx

       k = i+1  ! i'th equation corresponds to row k=i+1 of RHS
       
       if (id_row(k)) then
            psi(k) = 0.d0
      
        else
            ! row of matrix is not identity, but possibly flag this cell
            ! to set solution psi(i+1) = 0. after solving system:
            bouss_cell(i) = ((abs(eta(i)) <= bouss_wave_factor*h0(i)))
                        
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
            if (.false.) then
                write(66,*) '+++i,rhs(i+1): ',i, psi(k)
                write(66,*) (B_param+.5d0) * h02(i), s1_xx
                write(66,*) h03(i)/6.d0,  s1_h0_xx 
                write(66,*) B_param * g * h02(i) * h0_eta_xxx
                write(66,*) B_param,  g,  h02(i),  h0_eta_xxx
                write(66,*) s1(i)
                write(66,*) h0_eta_x(i)
                write(66,*) eta(i)
            endif
            
            if ((h0(i) > sw_depth0) .and. (h0(i) < sw_depth1)) then
                ! reduce psi linearly between sw_depth0 and sw_depth1:
                psi(k) = psi(k) * (h0(i)-sw_depth0)/(sw_depth1-sw_depth0)
            endif          
                               
        endif  ! if not id_row
    
    enddo
    
    if (.false.) then
        write(66,*) 'RHS:'
        do i=1,mx+2
            write(66,*) '+++i,rhs: ',i, psi(i)
        enddo
    endif
    
    ! -------------------------
    ! solve tridiagonal system.  
    ! psi is RHS on entry, solution on exit:
    
    LDB = mx+2  ! leading dimension of right-hand side matrix psi
    
    call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2, &
                 IPIV, psi, LDB,INFO )
    
    do i=1,mx
       if (.not. bouss_cell(i)) then
           ! do not apply bouss correction in this cell
           ! i'th equation corresponds to row k=i+1 of solution
           psi(i+1) = 0.d0
       endif
    enddo
                             
    if (.false.) then 
        ! version with debug output
        reset_count = 0  
        sw_row_count = 0                  
        do i=1,mx
            k = i+1  ! i'th equation corresponds to row k=i+1 of solution
            if (.not. bouss_cell(i)) then
                sw_row_count = sw_row_count + 1
                if (abs(psi(k))>0.1d0) then
                    !write(66,*) '+++ sw_row, k=',k,'  psi(k)=',psi(k)
                    reset_count = reset_count + 1
                endif
                psi(k) = 0.d0
            endif
        enddo
        
        if (reset_count > 0) then
            write(66,*) '+++ sw_row_count = ',sw_row_count, ' reset = ',reset_count
        endif
    endif
            
    if (.false.) then
        write(66,*) 'Solution psi: '
        do i=1,mx+2
            write(66,*) '+++ i,psi: ',i,psi(i)
        enddo
    endif

    return
end subroutine solve_tridiag


end module bouss_module


module bous_module

    ! Module parameters:

    implicit none
    logical      :: bouss  ! Turn on the dispersive terms
    real(kind=8) :: B_param 
    real(kind=8) :: sw_depth0, sw_depth1
    integer :: bc_xlo, bc_xhi
    
    save


contains

    subroutine set_bous(mthbc)


    ! Set Bparam and bc choices for Boussinesq implicit solver

    ! Eventually allow more things to be specified in setrun.py
    ! and read them in here.

    implicit none
    integer, intent(in) :: mthbc(2)

    ! Boussinesq equation parameter B, now set in setprob
    !Bparam = 1.d0/15.d0
    

    ! Boundary conditions to impose in computing Boussinesq update:
    if ((mthbc(1)==2) .or. (mthbc(2)==2)) then
        write(6,*) '*** Periodic BCs not supported in bouss_module'
        stop
    endif


    ! Use Dirichlet BC with value 0 for correction in ghost cells by default:
    bc_xlo = 0
    bc_xhi = 0
    
    if (mthbc(1)==3) then
        ! For wall boundary conditions at left boundary:
        bc_xlo = 3
    endif

    if (mthbc(2)==3) then
        ! For wall boundary conditions at right boundary:
        bc_xhi = 3
    endif

    ! To try out Neumann BCs:
    !bc_xlo = 1
    !bc_xhi = 1

    end subroutine set_bous
    

end module bous_module


subroutine setprob

   !set gauges and other parameters for 1D geoclaw
   ! for other problem set-up features, copy to your directory and modify

   use gauges_module, only: set_gauges
   use geoclaw_module, only: set_geo
   use bous_module, only: bouss, B_param, sw_depth0, sw_depth1

   use grid_module, only: read_grid, radial

   implicit none
   integer :: ndim, iunit
   character*25 fname
   common /comsrc/ ndim

   integer :: i

   call set_gauges(.false., 2) ! no restart, nvar=2
   call set_geo()


   iunit = 7
   fname = 'setprob.data'
!  # open the unit with new routine from Clawpack 4.4 to skip over
!  # comment lines starting with #:
   call opendatafile(iunit, fname)

   read(7,*) bouss
   read(7,*) B_param
   read(7,*) sw_depth0
   read(7,*) sw_depth1

   read(7,*) radial

   call read_grid('grid.data')

end subroutine setprob

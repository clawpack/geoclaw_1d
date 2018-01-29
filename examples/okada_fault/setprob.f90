
subroutine setprob

   !set gauges and other parameters for 1D geoclaw
   ! for other problem set-up features, copy to your directory and modify

   use setprob_module, only: imax0,hmax,mx_grid_max
   use gauges_module, only: set_gauges
   use geoclaw_module, only: set_geo
   use bous_module, only: bouss, B_param, sw_depth0, sw_depth1

   use grid_module, only: read_grid, radial

   implicit none
   integer :: ndim, iunit
   character*25 fname
   common /comsrc/ ndim

   integer :: i

   call set_gauges()
   call set_geo()

   imax0 = 1   ! where to start monitoring hmax
   write(6,*) 'monitoring hmax for i > imax0 = ',imax0

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

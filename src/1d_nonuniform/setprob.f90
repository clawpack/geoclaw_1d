
subroutine setprob

   !set gauges and other parameters for 1D geoclaw
   ! for other problem set-up features, copy to your directory and modify

   use gauges_module, only: set_gauges
   use geoclaw_module, only: set_geo

   use grid_module, only: read_grid, radial

   implicit none
   integer :: ndim, iunit
   common /comsrc/ ndim

   integer :: i

   call set_gauges(.false., 2) ! no restart, nvar=2
   call set_geo()

   call read_grid('grid.data')

end subroutine setprob

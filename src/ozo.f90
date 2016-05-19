program zo
  use mod_wrf_file
  use mod_time_step_loop
  implicit none

  character ( * ), parameter :: wrf_name = '../../Documents/wrf/data/wrf_omega_100km.nc'
  integer, parameter :: time_1 = 118, time_n = 120
  type ( wrf_file ) :: wrfin_file

  wrfin_file   = open_wrf_file   ( wrf_name )
  call time_step_loop ( wrfin_file, time_1, time_n )
  call close_wrf_file ( wrfin_file )


end program zo


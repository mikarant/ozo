program zo
  use mod_wrf_file
  use mod_time_step_loop
  implicit none

  character*140 :: infile, outfile
  character :: mode
  real :: alfa, toler
  integer :: time_1, time_n
  logical :: calc_htends,calc_omegas
  type ( wrf_file ) :: wrfin_file, omegafile

  namelist/PARAM/infile,outfile,alfa,toler,time_1,time_n,mode,calc_htends,&
       calc_omegas
  read(*,nml=PARAM)

  wrfin_file = open_wrf_file ( infile )
  if (calc_omegas) then
     omegafile = create_out_file ( outfile, wrfin_file, mode )
  else
     omegafile = open_out_file ( outfile )
  end if

  call time_step_loop ( wrfin_file, omegafile, time_1, time_n, alfa, toler, &
                        mode, calc_htends, calc_omegas)
  call close_wrf_file ( wrfin_file )
  call close_wrf_file ( omegafile )


end program zo


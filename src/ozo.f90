program zo
  use mod_wrf_file
  use mod_time_step_loop
  implicit none

  character*140 :: infile, outfile
  character :: mode
  real :: alfa, toler
  integer :: time_1, time_n, ny1, ny2
  logical :: calc_omegas, calc_b, debug
  type ( wrf_file ) :: wrfin_file, out_file

  namelist/PARAM/infile,outfile,alfa,toler,ny1,ny2,time_1,time_n,&
       mode,calc_omegas,debug
  read(*,nml=PARAM)
  calc_b=.false.

  if(mode.eq.'G')write(*,*)'Generalized omega equation'   
  if(mode.eq.'Q')write(*,*)'Quasi-geostrophic omega equation'   
  if(mode.eq.'T')write(*,*)'Generalized test version'   
  if(mode.eq.'t')write(*,*)'Quasigeostrophic test version'   
  if(mode.ne.'G'.and.mode.ne.'Q'.and.mode.ne.'T'.and.mode.ne.'t')then
     write(*,*)'Unknown mode of operation. Aborting'
     stop
  endif

  wrfin_file = open_wrf_file ( infile )
  if (calc_omegas) then
     out_file = create_out_file ( outfile, wrfin_file, mode, calc_b )
  else
     out_file = open_out_file ( outfile )
  end if

  call time_step_loop ( wrfin_file, out_file, time_1, time_n, alfa, toler, &
                        ny1, ny2, mode, calc_omegas, calc_b, debug)
  call close_wrf_file ( wrfin_file )
  call close_wrf_file ( out_file )


end program zo

program zo
  use mod_const
  use mod_wrf_file
  use mod_time_step_loop
  implicit none

  logical :: debug,calc_b
  type (wrf_file) :: wrfin_file, out_file
  type (parameters) :: param

  call read_parameters(param)
  calc_b=.false.

  if(param % mode.eq.'G')write(*,*)'Generalized omega equation'
  if(param % mode.eq.'Q')write(*,*)'Quasi-geostrophic omega equation'
  if(param % mode.eq.'T')write(*,*)'Generalized test version'
  if(param % mode.eq.'t')write(*,*)'Quasigeostrophic test version'
  if(param % mode.ne.'G'.and.param % mode.ne.'Q'.and.param % mode.ne. &
       'T'.and.param % mode.ne.'t')then
     write(*,*)'Unknown mode of operation. Aborting'
     stop
  endif

  write(*,*) 'Namelist parameters:'
  write(*,*) '--------------------'
  write(*,*) param
  write(*,*) '--------------------'


  wrfin_file = open_wrf_file ( param % infile )
  if (param % calc_omegas) then
     out_file = create_out_file ( param % outfile, wrfin_file, param % mode, &
          param % calc_b, param % forc )
  else
     out_file = open_out_file ( param % outfile )
  end if

  call time_step_loop ( wrfin_file, out_file, param, calc_b, debug)
  call close_wrf_file ( wrfin_file )
  call close_wrf_file ( out_file )


end program zo

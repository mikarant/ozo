module mod_const
  implicit none

  type parameters
     character*140 :: infile, outfile
     character :: mode
     real :: alfa, toler
     integer :: time_1, time_n, ny1, ny2
     logical :: calc_div, calc_b, calc_omegas, forc, ellipticity_correction
  end type parameters

  !   Threshold values to keep the generalized omega equation elliptic.
  real,parameter :: sigmamin=2e-7,etamin=2e-6

  real, parameter :: &
       r = 287.058, &       ! Gas constant for dry air
       cp = 1004., &        ! Specific heat of dry air
       g=9.80665            ! Gravitational acceleration

end module mod_const

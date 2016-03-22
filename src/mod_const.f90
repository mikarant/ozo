module mod_const
  implicit none

! This module contains only constants used in programs.
  
  real, parameter :: a = 6371e3,&        ! Earth radius
                     pi = 3.1415926536,& ! Pi
                     r = 287.,&          ! Gas constant for dry air
                     cp = 1004.,&        ! Specific heat of dry air
                     !omg = 7.292e-5,&    ! Angular speed of rotation of Earth
                     g=9.80665           ! Gravitational acceleration

end module mod_const

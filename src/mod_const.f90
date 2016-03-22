module mod_const
  implicit none

! This module contains only constants used in programs.
  
  real, parameter :: r = 287.,&          ! Gas constant for dry air
                     cp = 1004.,&        ! Specific heat of dry air
                     g=9.80665           ! Gravitational acceleration

end module mod_const

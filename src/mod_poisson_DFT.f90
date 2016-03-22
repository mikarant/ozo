module mod_poisson_DFT
  implicit none
  
contains

  subroutine poisson_solver_2D(f, dx, dy, phi, bd_ay, bd_by)
    use mod_poisson_green
    implicit none

    real, dimension ( :, : ), intent ( in ) :: f
    real, dimension ( :, : ), intent ( out ) :: phi
    double precision, dimension ( : ), intent ( in ), optional :: &
         bd_ay, bd_by
    real, intent ( in ) :: dx, dy
    character ( 3 ), parameter :: solver="DFT"
        
    select case ( solver )
       case ( "DFT" )
          call poisson_solver_DFT( f, dx, dy, bd_ay, bd_by, phi )
       case ( "Gre" )
          call poisson_solver_Green ( f, dx, dy, phi )
       end select
       
  end subroutine poisson_solver_2D
 
  subroutine poisson_solver_DFT ( rho, dx, dy, bd_ay, bd_by, phi )
    use mkl_poisson
    use mod_poisson_interp
    real, intent ( in ) :: dx, dy
    real, dimension ( :, : ), intent ( in )  :: rho
    real, dimension ( :, : ), intent ( out ) :: phi
    integer :: ipar(128), stat, nx, ny
    double precision :: q, ax, bx, ay, by
    double precision, dimension ( : ), intent ( in ) :: bd_ay, bd_by 
    double precision, dimension ( : ), allocatable :: &
         dpar, bd_ax, bd_bx
    character ( 4 ), parameter :: BCTYPE="PPDD"
    double precision, dimension ( :, : ), allocatable :: f_vc
    type(DFTI_DESCRIPTOR), pointer :: xhandle
    integer :: nlon, nlat
    double precision :: average
    integer, parameter :: bc = 9
    logical, parameter :: shift=.false.,interp=.false.

    nlon = size ( rho, 1)
    nlat = size ( rho, 2)
    nx = nlon
    ny = nlat
    allocate ( dpar ( 13 * nx / 2 + 7 ) )
    allocate ( bd_ax ( ny + 1 ), bd_bx ( ny + 1 ) )
    allocate ( f_vc ( nx + 1, ny + 1 ) )
    ax = 0.0e0
    bx = nx * dx
    ay = 0.0e0
    by = ny * dy
    Q = 0.0e0
    bd_ax = 0.0e0
    bd_bx = 0.0e0
    IPAR = 0
    f_vc = 0.0e0

    if ( interp ) then
       call interpolate_mc_to_vc( rho, f_vc, BCTYPE )
    else
       call shift_mc_to_vc ( rho, f_vc )
    end if

    f_vc = -1.0e0 * f_vc

    if ( shift ) then
       average = sum ( f_vc ( 1:nlon + 1, 1:nlat + 1 ) ) / &
            (nlon + 1) / (nlat + 1)
       f_vc = f_vc - average
    end if

    call d_init_Helmholtz_2D( AX, BX, AY, BY, NX, NY, BCTYPE, Q, &
         IPAR, dPAR, STAT )
    call d_commit_Helmholtz_2D( f_vc, bd_ax, bd_bx, bd_ay, bd_by, &
         xhandle, ipar, dpar, stat)
    call d_Helmholtz_2D(f_vc, bd_ax, bd_bx, bd_ay, bd_by, &
         xhandle, ipar, dpar, stat)
    call free_Helmholtz_2D(xhandle, ipar, stat)

    if ( interp ) then
       call interpolate_vc_to_mc( f_vc, phi )
    else
       call shift_vc_to_mc ( f_vc, phi )
    end if

  end subroutine poisson_solver_DFT

  subroutine shift_mc_to_vc ( f_mc, f_vc )
    implicit none

    real, dimension( :, : ), intent(in) :: f_mc
    double precision, dimension( :, : ), intent(out) :: f_vc
    integer :: nlon, nlat

    nlon = size( f_mc, 1)
    nlat = size( f_mc, 2)

    f_vc ( 1 : nlon, 1 : nlat ) = f_mc ( :, : )
    f_vc ( nlon + 1, : ) = f_vc ( 1, : )
    f_vc ( :, nlat + 1 ) = 0.0e0

  end subroutine shift_mc_to_vc

  subroutine shift_vc_to_mc ( f_vc, f_mc )
    implicit none

    double precision, dimension( :, : ), intent(in) :: f_vc
    real, dimension( :, : ), intent(out) :: f_mc

    integer :: nlon, nlat

    nlon = size( f_mc, 1)
    nlat = size( f_mc, 2)

    f_mc ( :, : ) = f_vc ( 1 : nlon, 1 : nlat )

  end subroutine shift_vc_to_mc


end module mod_poisson_DFT

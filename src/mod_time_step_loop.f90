module mod_time_step_loop
  use mod_wrf_file
  use mod_subrs
  implicit none

contains

  subroutine time_step_loop ( wrfin_file, time_1, time_n )
    type ( wrf_file ) :: wrfin_file
    integer :: time_1, time_n
    real, dimension ( :, :, : ), pointer :: T, u, v, z
    real, dimension ( :, :, : ), allocatable :: &
         dT_dt, du_dt, dv_dt, dz_dt, fx, fy, q, w
    real, dimension ( :, : ), allocatable :: mu_inv, p_sfc
    real, dimension ( :, :, :, : ), allocatable :: omegas, hTends

    integer :: time

    associate ( &
         nlon   => wrfin_file % dims ( 1 ), &
         nlat   => wrfin_file % dims ( 2 ), &
         nlev   => wrfin_file % dims ( 3 ), &
         dx     => wrfin_file % dx, &
         dy     => wrfin_file % dy, &
         p_levs => wrfin_file % pressure_levels, &
         corpar => wrfin_file % corpar )

      allocate ( hTends ( nlon, nlat, nlev, n_terms ) )

      call read_T_u_v_z ( wrfin_file, time_1 - 2 )
      call read_T_u_v_z ( wrfin_file, time_1 - 1 )
      do time = time_1, time_n
         print *, 'Time step: ', time
         call read_T_u_v_z ( wrfin_file, time )
         mu_inv = 1.0 / real2d ( wrfin_file, time, [ 'MU ', 'MUB' ]  )
         q      = diabatic_heating ( wrfin_file, time, mu_inv )
         fx     = friction ( wrfin_file, time, 'U', mu_inv )
         fy     = friction ( wrfin_file, time, 'V', mu_inv )
         p_sfc  = real2d ( wrfin_file, time, [ 'PSFC' ]  )
         w      = real3d ( wrfin_file, time, [ 'WW' ]  )
         omegas = read_omegas ( wrfin_file, time )
         call calculate_tendencies ( omegas, T, u, v, w, z, p_levs, &
              dx, dy, corpar, q, &
              fx,fy, dz_dt,  du_dt, &
              dv_dt, dT_dt, p_sfc, hTends )

         call write_tendencies ( wrfin_file, time, hTends )
         call write3d ( wrfin_file, time, ztend_name, dz_dt )

      end do

    end associate

  contains

    subroutine read_T_u_v_z ( file, time )
      type ( wrf_file ) :: file
      integer, intent ( in ) :: time
      real, dimension ( :, :, :, :, : ), allocatable, target, save :: dtvars
      character ( 3 ), dimension ( 4 ), parameter :: &
           dtvar_names = [ 'TT ', 'UU ', 'VV ', 'GHT' ]
      integer, save :: prev = 1, curr = 2, next = 3
      real :: dt2inv
      integer :: i
      if ( .not. allocated ( dtvars ) ) then
         associate ( &
              nlon => file % dims ( 1 ), &
              nlat => file % dims ( 2 ), &
              nlev => file % dims ( 3 ) )
           allocate ( dtvars ( nlon, nlat, nlev, 4, 3 ) )
         end associate
      end if
      prev = curr
      curr = next
      next = mod ( next, 3 ) + 1
      T => dtvars ( :, :, :, 1, curr )
      u => dtvars ( :, :, :, 2, curr )
      v => dtvars ( :, :, :, 3, curr )
      z => dtvars ( :, :, :, 4, curr )
      do i = 1, 4
         dtvars ( :, :, :, i, next ) = real3d ( file, time + 1, &
              [ trim ( dtvar_names ( i ) ) ] )
      end do
      if ( time .gt. lbound ( file % times, 1 ) ) then
         dt2inv = 1.0 / ( file % times ( time + 1 ) - file % times ( time - 1 ) )
         dT_dt = ( dtvars ( :, :, :, 1, next ) - dtvars ( :, :, :, 1,prev ) ) * dt2inv
         du_dt = ( dtvars ( :, :, :, 2, next ) - dtvars ( :, :, :, 2,prev ) ) * dt2inv
         dv_dt = ( dtvars ( :, :, :, 3, next ) - dtvars ( :, :, :, 3,prev ) ) * dt2inv
         dz_dt = ( dtvars ( :, :, :, 4, next ) - dtvars ( :, :, :, 4,prev ) ) * dt2inv
      end if
    end subroutine read_T_u_v_z

    function read_omegas ( file, time )
      type ( wrf_file ) :: file
      integer, intent ( in ) :: time
      real, dimension ( :, :, :, : ), allocatable :: read_omegas
      integer :: i
      associate ( nlon => file % dims ( 1 ), &
           nlat => file % dims ( 2 ), &
           nlev => file % dims ( 3 ), &
           nterms => size ( omega_term_names ) )
        allocate ( read_omegas ( nlon, nlat, nlev, nterms ) )
        do i = 1, nterms
           read_omegas ( :, :, :, i ) &
                = real3d ( file, time, [ trim ( omega_term_names ( i ) ) ] )
        end do
      end associate
    end function read_omegas

    function diabatic_heating ( file, time, mu_inv ) result ( q )
      type ( wrf_file ), intent ( in ) :: file
      integer, intent ( in ) :: time
      real, dimension ( :, : ), intent ( in ) :: mu_inv
      real, dimension ( :, :, : ), allocatable :: q
      integer :: lev
      real :: a_lev
      allocate ( q ( file % dims ( 1 ), file % dims ( 2 ), &
           file % dims ( 3 ) ) )
      associate ( plev => file % pressure_levels )
        select case ( file % wrf_cu_phys )
        case ( 1 )
           q = real3d ( file, time, &
                [ 'RTHSHTEN', 'RTHCUTEN', 'RTHRATEN', 'RTHBLTEN' ] )
        end select
        do lev = 1, size ( q, 3 )
           q ( :, :, lev ) = q ( :, :, lev ) * mu_inv
        end do
        q = q + real3d ( file, time, [ 'H_DIABATIC' ] )
        do lev = 1, size ( plev )
           a_lev = ( plev ( lev ) / 100000.0 ) ** ( 287.0 / 1004.0 )
           q ( :, :, lev ) = q ( :, :, lev ) * a_lev
        end do
      end associate
    end function diabatic_heating

    function friction ( file, time, direction, mu_inv )
      type ( wrf_file ), intent ( in ) :: file
      integer, intent ( in ) :: time
      character ( 1 ), intent ( in ) :: direction
      real, dimension ( :, : ), intent ( in ) :: mu_inv
      real, dimension ( :, :, : ), allocatable :: friction
      integer :: lev
      allocate ( friction ( file % dims ( 1 ), file % dims ( 2 ), &
           file % dims ( 3 ) ) )
      select case ( file % wrf_cu_phys )
      case ( 1 )
         friction = real3d ( file, time, &
              [ 'R'//direction//'SHTEN', 'R'//direction//'CUTEN', &
              'R'//direction//'BLTEN' ] )
      end select
      do lev = 1, size ( friction, 3 )
         friction ( :, :, lev ) = friction ( :, :, lev ) * mu_inv
      end do
    end function friction

    subroutine write_tendencies ( file, time, hTends )
      type ( wrf_file ), intent ( in ) :: file
      integer, intent ( in ) :: time
      real, dimension ( :, :, :, : ) :: hTends
      integer :: t, varid
      associate ( &
           ncid => file % ncid, &
           nlon => file % dims ( 1 ), &
           nlat => file % dims ( 2 ), &
           nlev => file % dims ( 3 ) )
        do t = 1, size ( htend_term_names )
           call check ( nf90_inq_varid ( ncid, &
                trim ( htend_term_names ( t ) ), varid ) )
           call check ( nf90_put_var ( ncid, varid, &
                hTends ( :, :, :, t ), &
                [ 1, 1, 1, time ], &
                [ nlon, nlat, nlev, 1 ] ) )
        end do
      end associate
    end subroutine write_tendencies

    subroutine write3d ( file, time, name, data )
      type ( wrf_file ), intent ( in ) :: file
      integer, intent ( in ) :: time
      character(*) :: name
      real, dimension ( :, :, : ) :: data
      integer :: varid

      associate ( &
           ncid => file % ncid, &
           nlon => file % dims ( 1 ), &
           nlat => file % dims ( 2 ), &
           nlev => file % dims ( 3 ) )
        
        call check ( nf90_inq_varid ( ncid, &
             trim ( name ), varid ) )
        call check ( nf90_put_var ( ncid, varid, &
             data ( :, :, : ), &
             [ 1, 1, 1, time ], &
             [ nlon, nlat, nlev, 1 ] ) )

      end associate

    end subroutine write3d
  end subroutine time_step_loop

end module mod_time_step_loop

module mod_poisson_interp
  implicit none 

contains
  
  subroutine interpolate_mc_to_vc( f_mc, f_vc, BCTYPE )
    implicit none

    real, dimension( :, : ), intent(in) :: f_mc
    double precision, dimension( :, : ), intent(inout) :: f_vc
    character ( 4 ) ,intent ( in ) :: BCTYPE
    integer :: nlon, nlat, i, j

    nlon = size( f_mc, 1)
    nlat = size( f_mc, 2)

    do i = 2, nlon
       do j = 2, nlat
          f_vc ( i, j ) = ( f_mc ( i-1, j-1 ) + &
               f_mc ( i-1, j ) + &
               f_mc ( i, j-1 ) + &
               f_mc ( i, j ) ) / 4
       enddo
    enddo

    do j = 2, nlat
       f_vc ( 1, j ) = ( f_mc ( 1, j-1 ) + &
            f_mc ( 1, j ) + &
            f_mc ( nlon, j-1 ) + &
            f_mc ( nlon, j ) ) / 4
       f_vc ( nlon + 1, j ) = f_vc ( 1, j )
    enddo

    select case ( BCTYPE )

    case ( "PPPP" ) ! Periodic y-dimension
       
       do i = 2, nlon
          f_vc ( i, 1 ) = ( f_mc ( i, 1) + &
               f_mc ( i-1, 1 ) + &
               f_mc ( i, nlat ) + &
               f_mc ( i-1, nlat ) ) / 4
          f_vc ( i, nlat + 1 ) = f_vc ( i, 1)
       enddo

       f_vc ( 1, 1 ) = ( f_mc ( 1, 1 ) + &
            f_mc ( nlon, 1 ) + &
            f_mc ( 1, nlat ) + &
            f_mc ( nlon, nlat ) ) / 4
       f_vc ( nlon +1 , 1 ) = f_vc ( 1, 1 )
       f_vc ( 1, nlat + 1 ) = f_vc ( 1, 1 )
       f_vc ( nlon + 1, nlat + 1 ) = f_vc ( 1, 1 )

    case ( "PPDD", "PPNN" ) ! Diriclet y-dimension
       
       do i = 2, nlon
          f_vc ( i, 1 ) = ( f_mc ( i, 1 ) + &
               f_mc ( i-1, 1 ) ) / 2
          f_vc ( i, nlat + 1 ) = ( f_mc ( i, nlat ) + &
               f_mc ( i-1, nlat ) ) / 2
       enddo
       
       f_vc ( 1, 1 ) =  ( f_mc ( 1, 1 ) + &
            f_mc ( nlon, 1 ) ) / 2
       f_vc ( nlon + 1, 1 ) = f_vc ( 1, 1 )
       f_vc ( 1, nlat + 1 ) =  ( f_mc ( 1, nlat ) + &
            f_mc ( nlon, nlat ) ) / 2
       f_vc ( nlon + 1, nlat + 1 ) = f_vc ( 1, nlat + 1 )

    end select

  end subroutine interpolate_mc_to_vc

  subroutine interpolate_vc_to_mc( f_vc, f_mc )
    implicit none

    double precision, dimension ( :, : ), intent ( in ) :: f_vc
    real, dimension ( :, : ), intent ( out ) :: f_mc

    integer :: nlon, nlat, i, j

    nlon = size( f_mc, 1)
    nlat = size( f_mc, 2)

    do i = 1, nlon
       do j = 1, nlat
          f_mc ( i, j ) = ( f_vc ( i, j ) + &
               f_vc ( i + 1, j ) + &
               f_vc ( i, j + 1 ) + &
               f_vc ( i + 1, j + 1 ) ) / 4
       enddo
    enddo

  end subroutine interpolate_vc_to_mc

end module mod_poisson_interp

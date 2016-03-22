module mod_poisson
  implicit none
  
contains

  subroutine poisson_solver_2D(f, BCTYPE, shift, solver, dx, dy, phi, &
        bd_ay, bd_by)
    implicit none

    real, dimension ( :, : ), intent ( in ) :: f
    real, dimension ( :, : ), intent ( out ) :: phi
    double precision, dimension ( : ), intent ( in ), optional :: &
         bd_ay, bd_by
    real, intent ( in ) :: dx, dy
    logical, intent ( in ) :: shift
    character ( 4 ), intent ( in ) :: BCTYPE
    character ( 3 ), intent ( in ) :: solver
        
    select case ( solver )
       case ( "DFT" )
          call poisson_solver_DFT( f, dx, dy, bd_ay, bd_by, BCTYPE, shift, phi )
       case ( "Gre" )
          call poisson_solver_Green ( f, dx, dy, phi )
       end select
       
  end subroutine poisson_solver_2D
 
  subroutine poisson_solver_DFT ( rho, dx, dy, bd_ay, bd_by, BCTYPE, &
       shift, phi )
    use mkl_poisson
    real, intent ( in ) :: dx, dy
    real, dimension ( :, : ), intent ( in )  :: rho
    real, dimension ( :, : ), intent ( out ) :: phi
    integer :: ipar(128), stat, nx, ny
    double precision :: q, ax, bx, ay, by
    double precision, dimension ( : ), intent ( in ) :: bd_ay, bd_by 
    double precision, dimension ( : ), allocatable :: &
         dpar, bd_ax, bd_bx
    character ( 4 ) ,intent ( in ) :: BCTYPE
    double precision, dimension ( :, : ), allocatable :: f_vc
    type(DFTI_DESCRIPTOR), pointer :: xhandle
    integer :: nlon, nlat
    double precision :: average
    integer, parameter :: bc = 9
    logical, intent ( in ) :: shift

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

!    call interpolate_mc_to_vc( rho, f_vc, BCTYPE )
    call shift_mc_to_vc ( rho, f_vc )

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

!    call interpolate_vc_to_mc ( f_vc , phi )
    call shift_vc_to_mc ( f_vc, phi )

  end subroutine poisson_solver_DFT

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

  subroutine poisson_solver_Green ( f, dx, dy, hgttend )
!
!   New version of "height_tendency", using precalculated response kernel
!   (JR 10.12.2015)
 
    implicit none
    real, dimension ( :, :, : ), allocatable, save  :: response
    real, dimension ( :, : ), intent ( in ) :: f
    real,                     intent ( in ) :: dx, dy
    real, dimension ( :, : ), intent ( out ) :: hgttend
    integer :: i,j
    integer :: ind, ii, jj, nlon, nlat

    nlon=size(f,1)
    nlat=size(f,2)

    if(.not. allocated (response))then
       allocate(response(nlon/2+1,nlat,(nlat+1)/2))
       call responsekernel(response,nlon,nlat,dx,dy)
    end if
    
    hgttend=0.0e0

    do i=1,nlon
       do j=2,nlat-1
          do ii=1,nlon

             ind = 1 + min ( abs(i-ii), abs(i-ii+nlon), abs(i-ii-nlon) )

             do jj=1,nlat/2
                hgttend(i,j) = hgttend(i,j) + f(ii,jj) &
                     * response(ind,j,jj)
             enddo

             do jj=nlat/2+1,nlat
                hgttend(i,j) = hgttend(i,j) + f(ii,jj) &
                     *response(ind,nlat+1-j,nlat+1-jj)
             enddo

          enddo
       enddo
    enddo

  end subroutine poisson_solver_Green

  subroutine responsekernel(response,nlon,nlat,dx,dy) 

    implicit none
    real,dimension(:,:,:),intent(inout) :: response
    integer,intent(in) :: nlon,nlat
    real,intent(in) :: dx,dy

    real :: forcing(nlon,nlat,1),res1(nlon,nlat,1)
    integer :: i,j,j0

    do j0=1,(nlat+1)/2
       forcing=0.
       forcing(1,j0,1)=1.
       call invert_lapl(forcing,nlon,nlat,1,dx,dy,res1) 
       do i=1,nlon/2+1
          do j=1,nlat
             response(i,j,j0)=res1(i,j,1)
          enddo
       enddo
       write(*,*)'Calculated response, j0=',j0
    enddo

  end subroutine responsekernel

   subroutine invert_lapl(forcing,nlon,nlat,nlev,dx,dy,phi)
!
!   Former "height_tendency". Only name of subroutine changed (JR 10.12.15)
!
!   Calculation of inverse of laplacian with iteration and by using function laplace2_cart.
!   Southern and northern boundaries are left as zero.
!
    implicit none
    real,dimension(:,:,:),intent(in) :: forcing
    real,intent(in) :: dx,dy
    integer,intent(in) :: nlev
    real,dimension(:,:,:),intent(inout) :: phi

    integer :: nlon,nlat,l,i,j,k,itermax
    real :: coeff(nlon,nlat,nlev),lapl2(nlon,nlat,nlev)
     
    phi=0.
    itermax=1000

    do l=1,itermax
       call laplace2_cart(phi,lapl2,coeff,dx,dy)
       do k=1,nlev
          do j=2,nlat-1
             do i=1,nlon
                phi(i,j,k)=(forcing(i,j,k)-lapl2(i,j,k))/coeff(i,j,k)
             enddo
          enddo
       enddo
    enddo
    
  end subroutine invert_lapl

  subroutine laplace2_cart(f,lapl2,coeff,dx,dy)
!
!     As laplace_cart
!       - jätetään kussakin pisteessä paikallisen arvon osuus pois
!       - lasketaan paikallisen arvon kerroin coeff
!
    real,dimension(:,:,:),intent(in) :: f
    real,    intent ( in )  :: dx, dy
    real,dimension(:,:,:), intent ( out ) :: lapl2, coeff
    integer :: nlon, nlat, nlev
 
    nlon=size(f,1)
    nlat=size(f,2)
    nlev=size(f,3)

    ! x-direction
    lapl2 ( 2 : nlon - 1, :, : ) = f( 1: nlon - 2, :, : ) + f ( 3: nlon, :, : )
    lapl2 ( 1, :, : )    = f( nlon,     :, : ) + f ( 2, :, : )
    lapl2 ( nlon, :, : ) = f( nlon - 1, :, : ) + f ( 1, :, : )
    lapl2 = lapl2 / ( dx * dx )

    ! y-directon
    lapl2 ( :, 2 : nlat -1, : ) = lapl2 ( :, 2 : nlat -1, : ) &
         + ( f ( :, 1 : nlat -2, : ) + f ( :, 3 : nlat, :) ) / ( dy * dy )

    coeff ( :, 2 : nlat -1, : ) = -2.0 / (dx * dx ) -2.0 / (dy * dy )
    coeff ( :, 1,    : ) = -2.0 / (dx * dx )
    coeff ( :, nlat, : ) = -2.0 / (dx * dx )
 
  end subroutine laplace2_cart

end module mod_poisson

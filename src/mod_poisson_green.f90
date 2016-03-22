module mod_poisson_green
  implicit none

contains

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


end module mod_poisson_green

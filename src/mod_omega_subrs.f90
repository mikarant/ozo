module mod_omega_subrs
  use mod_const
  use mod_common_subrs
  implicit none

contains

!
!****************** SUBROUTINES **************************************************
!
  subroutine QG_test(omegaan,sigma,feta,dx,dy,dlev,ftest)
!   Forcing for quasigeostrophic test case ('t')
!   In essence: calculating the LHS from the WRF omega (omegaan) 
!   and substituting it to the RHS

    implicit none

    real,dimension(:,:,:,:),intent(in) :: sigma,feta
    real,dimension(:,:,:),intent(in) :: omegaan
    real,intent(in) :: dx,dy,dlev
    real,dimension(:,:,:,:),intent(out) :: ftest 

    real,dimension(:,:,:),allocatable :: df2dp2,lapl
    integer :: nlon,nlat,nlev

    nlon=size(omegaan,1)
    nlat=size(omegaan,2)
    nlev=size(omegaan,3)
    allocate(df2dp2(nlon,nlat,nlev),lapl(nlon,nlat,nlev))

    call p2der(omegaan,dlev,df2dp2)
    call laplace_cart(omegaan,lapl,dx,dy)

    ftest(:,:,:,1)=sigma(:,:,:,1)*lapl+feta(:,:,:,1)*df2dp2 

  end subroutine QG_test

  subroutine gen_test(sigmaraw,omegaan,zetaraw,corpar,dx,dy,dlev,ftest)

!   Forcing for the general test case
!   In essence: calculating the LHS from the WRF omega (omegaan) 
!   and substituting it to the RHS

    implicit none

    real,dimension(:,:,:),intent(in) :: sigmaraw,omegaan,zetaraw
    real,dimension(:),intent(in) :: corpar
    real,intent(in) :: dx,dy,dlev
    real,dimension(:,:,:,:),intent(out) :: ftest 
    real,dimension(:,:,:),allocatable :: dum2,dum1,dum3,dum4,dum5,dum6
    real,dimension(:,:,:),allocatable :: dvdp,dudp

    integer :: nlon,nlat,nlev,j
    nlon=size(omegaan,1)
    nlat=size(omegaan,2)
    nlev=size(omegaan,3)

    allocate(dum2(nlon,nlat,nlev),dum1(nlon,nlat,nlev))
    allocate(dum3(nlon,nlat,nlev),dum4(nlon,nlat,nlev))
    allocate(dum5(nlon,nlat,nlev),dum6(nlon,nlat,nlev))
    allocate(dvdp(nlon,nlat,nlev),dudp(nlon,nlat,nlev))

    dum2=sigmaraw*omegaan
          
    call laplace_cart(dum2,dum1,dx,dy)
    call p2der(omegaan,dlev,dum3)
    call p2der(zetaraw,dlev,dum5)
    
    do j=1,nlat       
       dum2(:,j,:)=(corpar(j)+zetaraw(:,j,:))*corpar(j)*dum3(:,j,:)
    enddo
    
    do j=1,nlat       
       dum3(:,j,:)=-corpar(j)*omegaan(:,j,:)*dum5(:,j,:)
    enddo
    
    call xder_cart(omegaan,dx,dum4) 
    call yder_cart(omegaan,dy,dum5) 
    
    do j=1,nlat
       dum6(:,j,:)=-corpar(j)*(dvdp(:,j,:)*dum4(:,j,:)-dudp(:,j,:)*dum5(:,j,:))
    enddo
    call pder(dum6,dlev,dum4) 
    
    ftest(:,:,:,1)=dum1+dum2+dum3+dum4 

  end subroutine gen_test

  subroutine coarsen3D(f,g,nlon1,nlat1,nlev1,nlon2,nlat2,nlev2)
!
!      Averages the values of field f (grid size nlon1 x nlat1 x nlev1) over larger
!      grid boxes (grid size nlon2 x nlat2 x nlev2), to field g

!      To keep the algorithm
!      simple, only 0/1 weights are used -> works only well if
!      div(nlon1/nlon2)=div(nlat1/nlat2)=div(nlev1/nlev2)
!
    implicit none
    integer,intent(in) :: nlon1,nlat1,nlon2,nlat2,nlev1,nlev2
    real,dimension(nlon1,nlat1,nlev1),intent(in) :: f
    real,dimension(nlon2,nlat2,nlev2),intent(out) :: g
    integer :: i,i2,j,j2,k,k2,imin,imax,jmin,jmax,kmin,kmax
    real :: fsum
    
    do i2=1,nlon2
       imin=nint((i2-1)*real(nlon1)/real(nlon2)+1)
       imax=nint(i2*real(nlon1)/real(nlon2))
       do j2=1,nlat2
          jmin=nint((j2-1)*real(nlat1)/real(nlat2)+1)
          jmax=nint(j2*real(nlat1)/real(nlat2))
          do k2=1,nlev2
             kmin=nint((k2-1)*real(nlev1)/real(nlev2)+1)
             kmax=nint(k2*real(nlev1)/real(nlev2))
             fsum=0.
             do i=imin,imax
                do j=jmin,jmax
                   do k=kmin,kmax
                      fsum=fsum+f(i,j,k)
                   enddo
                enddo
             enddo
             g(i2,j2,k2)=fsum/((imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1))
          enddo
       enddo
    enddo
    
  end subroutine coarsen3D

  subroutine finen3D(f,g,nlon1,nlat1,nlev1,nlon2,nlat2,nlev2)
!
!      Distributes the values of field f (grid size nlon2 x nlat2 x nlev2) to a
!      finer grid g (grid size nlon1 x nlat1 x nlev1), assuming that f is
!      constant in each grid box of the original grid 
!
!      ** PERHAPS THIS SHOULD BE REPLACED BY BILINEAR INTERPOLATION 
!      ** TO AVOID ARTIFICIAL JUMPS

    implicit none
    integer,intent(in) :: nlon1,nlat1,nlev1,nlon2,nlat2,nlev2
    real,intent(in) :: f(nlon2,nlat2,nlev2)
    real,intent(out) :: g(nlon1,nlat1,nlev1)
    integer :: i,i2,j,j2,k,k2,imin,imax,jmin,jmax,kmin,kmax

    do i2=1,nlon2
       imin=nint((i2-1)*real(nlon1)/real(nlon2)+1)
       imax=nint(i2*real(nlon1)/real(nlon2))
       do j2=1,nlat2
          jmin=nint((j2-1)*real(nlat1)/real(nlat2)+1)
          jmax=nint(j2*real(nlat1)/real(nlat2))
          do k2=1,nlev2
             kmin=nint((k2-1)*real(nlev1)/real(nlev2)+1)
             kmax=nint(k2*real(nlev1)/real(nlev2))
             do i=imin,imax
                do j=jmin,jmax
                   do k=kmin,kmax
                      g(i,j,k)=f(i2,j2,k2)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    
  end subroutine finen3D

  subroutine gwinds(z,dx,dy,corpar,u,v)
!
!      Calculation of geostrophic winds (u,v) from z. At the equator, mean of
!      the two neighbouring latitudes is used (should not be a relevant case).
!
    implicit none

    real,dimension(:,:,:),intent(in) :: z
    real,dimension(:,:,:),intent(out) :: u,v
    real,dimension(:),intent(in) :: corpar
    real,intent(in) :: dx,dy
    integer :: nlon,nlat,nlev,i,j,k
    real,dimension(:,:,:),allocatable :: dzdx,dzdy

    nlon=size(u,1)
    nlat=size(u,2)
    nlev=size(u,3)
    allocate(dzdx(nlon,nlat,nlev),dzdy(nlon,nlat,nlev))

    call xder_cart(z,dx,dzdx) 
    call yder_cart(z,dy,dzdy) 

    do k=1,nlev
       do j=1,nlat
          if(abs(corpar(j)).gt.1e-7)then
             do i=1,nlon
                u(i,j,k)=-g*dzdy(i,j,k)/corpar(j) 
                v(i,j,k)=g*dzdx(i,j,k)/corpar(j) 
             enddo
          endif
       enddo
       do j=1,nlat
          if(abs(corpar(j)).lt.1e-7)then
             do i=1,nlon
                u(i,j,k)=(u(i,j+1,k)+u(i,j-1,k))/2.
                v(i,j,k)=(v(i,j+1,k)+v(i,j-1,k))/2.
             enddo
          endif
       enddo
    enddo

  end subroutine gwinds

  subroutine modify(sigmaraw,sigmamin,etamin,zetaraw,&
                        corpar,dudp,dvdp,sigma,feta,zeta)      
!
!      Modifying stability and vorticity to keep the LHS of the genearlized
!      omega equation elliptic
!
    implicit none

    real,dimension(:,:,:),intent(in) :: sigmaraw,zetaraw,dudp,dvdp
    real,dimension(:),intent(in) :: corpar
    real,dimension(:,:,:),intent(inout) :: zeta,feta,sigma
    real,intent(in) :: sigmamin,etamin
    integer :: nlon,nlat,nlev,i,j,k
    
    nlon=size(sigmaraw,1)
    nlat=size(sigmaraw,2)
    nlev=size(sigmaraw,3)
     
    do k=1,nlev
       do j=1,nlat
          do i=1,nlon
             sigma(i,j,k)=max(sigmaraw(i,j,k),sigmamin)
             zeta(i,j,k)=0.
             ! Northern Hemisphere
             if(corpar(j).gt.1e-7)then
                zeta(i,j,k)=max(zetaraw(i,j,k),etamin+corpar(j)/ &
                            (4*sigma(i,j,k))*(dudp(i,j,k)**2.+dvdp(i,j,k)**2.) &
                            -corpar(j))
             endif
             ! Southern Hemisphere
             if(corpar(j).lt.-1e-7)then
                zeta(i,j,k)=min(zetaraw(i,j,k),-etamin+corpar(j)/ &
                            (4*sigma(i,j,k))*(dudp(i,j,k)**2.+dvdp(i,j,k)**2.) &
                            -corpar(j))
             endif
             feta(i,j,k)=(zeta(i,j,k)+corpar(j))*corpar(j)
          enddo
       enddo
    enddo

  end subroutine modify

  subroutine aave(f,res)                      
!
!     Calculation of area mean (res) of field f in cartesian coordinates.
!     Simplest possible way.
!
    implicit none

    real,dimension(:,:),intent(in) :: f
    real :: res,sum,wsum
    integer :: i,j,nlon,nlat
    nlon=size(f,1)
    nlat=size(f,2)

    sum=0
    wsum=0
    do j=1,nlat
       do i=1,nlon
          sum=sum+f(i,j)
          wsum=wsum+1.
       enddo
    enddo
    res=sum/wsum

  end subroutine aave

  subroutine fvort(u,v,zeta,corpar,dx,dy,dp,mulfact,fv)
!
!     Calculation of vorticity advection forcing
!     Input: u,v,zeta
!     Output: stored in "fv"
!
    implicit none
    real,dimension(:,:,:),intent(in) :: u,v,zeta,mulfact
    real,dimension(:,:,:),intent(out) :: fv
    real,dimension(:),intent(in) :: corpar
    real,dimension(:,:,:),allocatable :: eta,adv,dadvdp
    real,intent(in) :: dx,dy,dp
    integer :: i,j,k,nlon,nlat,nlev

    nlon=size(u,1)
    nlat=size(u,2)
    nlev=size(u,3)
    allocate(eta(nlon,nlat,nlev),adv(nlon,nlat,nlev))
    allocate(dadvdp(nlon,nlat,nlev))

    do k=1,nlev
       do j=1,nlat
          do i=1,nlon
             eta(i,j,k)=zeta(i,j,k)+corpar(j)
          enddo
       enddo
    enddo

    call advect_cart(u,v,eta,dx,dy,adv)
    adv=adv*mulfact
    call pder(adv,dp,dadvdp) 
!    print*,"fvort"
!    write(*,"(20E11.3E2)")sum(dadvdp)
    do k=1,nlev
       do j=1,nlat      
          do i=1,nlon
             fv(i,j,k)=corpar(j)*dadvdp(i,j,k)
          enddo
       enddo
    enddo

  end subroutine fvort

  subroutine ftemp(u,v,t,lev,dx,dy,mulfact,ft)
!
!     Calculation of temperature advection forcing
!     Input: u,v,t
!     Output: stored in "adv" (bad style ...)
!
    implicit none
    real,dimension(:,:,:),intent(in) :: u,v,t,mulfact
    real,dimension(:,:,:),intent(out) :: ft
    real,dimension(:),intent(in) :: lev
    real,intent(in) :: dx,dy
    real,dimension(:,:,:),allocatable :: adv,lapladv

    integer :: i,j,k,nlon,nlat,nlev

    nlon=size(u,1)
    nlat=size(u,2)
    nlev=size(u,3)
    allocate(adv(nlon,nlat,nlev),lapladv(nlon,nlat,nlev))

    call advect_cart(u,v,t,dx,dy,adv)
    adv=adv*mulfact
    call laplace_cart(adv,lapladv,dx,dy)         

    do k=1,nlev
       do j=1,nlat      
          do i=1,nlon
             ft(i,j,k)=lapladv(i,j,k)*r/lev(k)
          enddo
       enddo
    enddo

  end subroutine ftemp

  subroutine ffrict(fx,fy,corpar,dx,dy,dp,mulfact,ff)
!
!     Calculation of friction forcing
!     Input: fx,fy = x and y components of "friction force"
!     Output: res 
!
    implicit none
    real,dimension(:,:,:),intent(in) :: fx,fy,mulfact
    real,dimension(:),intent(in) :: corpar
    real,dimension(:,:,:),intent(out) :: ff
    real,dimension(:,:,:),allocatable :: fcurl,dcurldp
    real,intent(in) :: dx,dy,dp
    integer :: i,j,k,nlon,nlat,nlev

    nlon=size(fx,1)
    nlat=size(fx,2)
    nlev=size(fx,3)
    allocate(fcurl(nlon,nlat,nlev),dcurldp(nlon,nlat,nlev))

    call curl_cart(fx,fy,dx,dy,fcurl)
    fcurl=fcurl*mulfact
    call pder(fcurl,dp,dcurldp) 

    do k=1,nlev
       do j=1,nlat      
          do i=1,nlon
             ff(i,j,k)=-corpar(j)*dcurldp(i,j,k)
          enddo
       enddo
    enddo

  end subroutine ffrict

  subroutine fdiab(q,lev,dx,dy,mulfact,fq)
!
!     Calculation of diabatic heaging forcing
!     Input: q = diabatic temperature tendency (already normalized by cp)
!     Output: stored in "fq"
!
    implicit none
    real,dimension(:,:,:),intent(inout) :: q,mulfact
    real,dimension(:,:,:),intent(out) :: fq
    real,dimension(:),intent(in) :: lev
    real,intent(in) :: dx,dy      
    integer :: i,j,k,nlon,nlat,nlev

    nlon=size(q,1)
    nlat=size(q,2)
    nlev=size(q,3)
    
    q=q*mulfact
    call laplace_cart(q,fq,dx,dy)         

    do k=1,nlev
       do j=1,nlat      
          do i=1,nlon
             fq(i,j,k)=-r*fq(i,j,k)/lev(k)
          enddo
       enddo
    enddo
    
  end subroutine fdiab

  subroutine fimbal(dzetadt,dtdt,corpar,lev,dx,dy,dp,mulfact,fa)
!
!     Calculation of the FA ("imbalance") forcing term
!     Input: dzetadt, dtdt = vorticity & temperature tendencies
!     Output: fa
!
    implicit none
    real,dimension(:,:,:),intent(inout) :: dzetadt,dtdt
    real,dimension(:,:,:),intent(in) ::  mulfact
    real,dimension(:,:,:),intent(out) :: fa
    real,dimension(:),intent(in) :: lev,corpar
    real,intent(in) :: dx,dy,dp
    real,dimension(:,:,:),allocatable :: ddpdzetadt,lapldtdt
    integer i,j,k,nlon,nlat,nlev

    nlon=size(dtdt,1)
    nlat=size(dtdt,2)
    nlev=size(dtdt,3)
    allocate(ddpdzetadt(nlon,nlat,nlev),lapldtdt(nlon,nlat,nlev))

    dzetadt=dzetadt*mulfact
    dtdt=dtdt*mulfact

    call pder(dzetadt,dp,ddpdzetadt) 
    call laplace_cart(dtdt,lapldtdt,dx,dy)         

    do k=1,nlev
       do j=1,nlat      
          do i=1,nlon
             ddpdzetadt(i,j,k)=corpar(j)*ddpdzetadt(i,j,k)
             fa(i,j,k)=ddpdzetadt(i,j,k)+lapldtdt(i,j,k)*r/lev(k)
          enddo
       enddo
    enddo
    
  end subroutine fimbal

  subroutine callsolveQG(rhs,boundaries,omega,nlonx,nlatx,nlevx,&
                         dx,dy,dlev,sigma0,feta,nres,alfa,toler)
!
!   Calling solveQG + writing out omega. Multigrid algorithm.
!
    implicit none
    integer,intent(in) :: nres
    integer,dimension(:),intent(in) :: nlonx,nlatx,nlevx
    real,dimension(:),intent(in) :: dx,dy,dlev
    real,dimension(:,:,:,:),intent(inout) :: rhs,omega
    real,dimension(:,:,:,:),intent(in) :: boundaries,feta
    real,dimension(:,:),intent(in) :: sigma0
    real,intent(in) :: alfa,toler
    real,dimension(:,:,:,:),allocatable :: omegaold
    real,dimension(:,:,:),allocatable :: omega1,dum1,resid
    real :: maxdiff,aomega
    integer :: iter,i,j,k,ires

    integer,parameter :: itermax=1000
    integer,parameter :: ny1=2,ny2=2 ! number of iterations at each grid 
                                     ! resolution when proceeding to coarser 
                                     ! (ny1) and when returning to finer (ny2)
    logical,parameter :: lzeromean=.true. ! Area means of omega are set to zero

    allocate(omegaold(nlonx(1),nlatx(1),nlevx(1),nres))
    allocate(omega1(nlonx(1),nlatx(1),nlevx(1)))
    allocate(dum1(nlonx(1),nlatx(1),nlevx(1)))
    allocate(resid(nlonx(1),nlatx(1),nlevx(1)))
    omega=0.
    omegaold=boundaries
!
!   The whole multigrid cycle is written explicitly here. Better as a separate subroutine?  
!
!----------------------------------------------------------------------------------------
    do iter=1,itermax  ! Each iteration = one (fine->coarse->fine) multigrid cycle
!
!   Loop from finer to coarser resolutions
!
       do ires=1,nres
!         write(*,*)'fine-coarse:iter,ires',iter,ires
          call solveQG(rhs(:,:,:,ires),boundaries(:,:,:,ires),& 
               omega(:,:,:,ires),omegaold(:,:,:,ires),nlonx(ires),nlatx(ires),nlevx(ires),&
               dx(ires),dy(ires),dlev(ires),sigma0(:,ires),feta(:,:,:,ires),&
               ny1,alfa,.true.,resid)
          if(ires.eq.1)omega1(:,:,:)=omega(:,:,:,1)
          if(ires.lt.nres)then
             call coarsen3d(resid,rhs(:,:,:,ires+1),nlonx(ires),nlatx(ires),nlevx(ires),&
                  nlonx(ires+1),nlatx(ires+1),nlevx(ires+1))          
          endif
       enddo
!
!      Loop from coarser to finer resolutions
!
       do ires=nres-1,1,-1
!        write(*,*)'coarse-fine:iter,ires',iter,ires
          call finen3D(omega(:,:,:,ires+1),dum1,nlonx(ires),nlatx(ires),nlevx(ires),& 
               nlonx(ires+1),nlatx(ires+1),nlevx(ires+1))          
!        Without the underrelaxation (coefficient alfa), the solution diverges
          omegaold(:,:,:,ires)=omega(:,:,:,ires)+alfa*dum1(:,:,:)
!         if(ires.eq.1)then
!         write(*,*)'omega',ires,omega(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2,ires)
!         write(*,*)'dum1',ires,dum1(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2)
!         write(*,*)'omegaold',ires,omegaold(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2,ires)
!         endif

          call solveQG(rhs(:,:,:,ires),boundaries(:,:,:,ires),& 
               omega(:,:,:,ires),omegaold(:,:,:,ires),nlonx(ires),nlatx(ires),nlevx(ires),&
               dx(ires),dy(ires),dlev(ires),sigma0(:,ires),feta(:,:,:,ires),&
               ny2,alfa,.false.,resid)
       enddo

       maxdiff=0.
       do k=1,nlevx(1)  
          do j=1,nlatx(1)  
             do i=1,nlonx(1)  
                maxdiff=max(maxdiff,abs(omega(i,j,k,1)-omega1(i,j,k)))
             enddo
          enddo
       enddo
       print*,iter,maxdiff
       if(maxdiff.lt.toler.or.iter.eq.itermax)then
!          write(*,*)'iter,maxdiff',iter,maxdiff
          goto 10
       endif
       
       omegaold=omega
          
    enddo ! iter=1,itermax
10  continue         
!----------------------------------------------------------------------------------------
!
!       Subtract the area mean of omega
!
    if(lzeromean)then
       do k=1,nlevx(1) 
          call aave(omega(:,:,k,1),aomega)
          do j=1,nlatx(1)
             do i=1,nlonx(1)
                omega(i,j,k,1)=omega(i,j,k,1)-aomega
             enddo
          enddo
       enddo
    endif
    
  end subroutine callsolveQG

  subroutine solveQG(rhs,boundaries,omega,omegaold,nlon,nlat,nlev,&
       dx,dy,dlev,sigma0,feta,niter,alfa,lres,resid)
!
!   Solving the QG omega equation using 'niter' iterations.
!
    implicit none

    integer,intent(in) :: nlon,nlat,nlev,niter
    real,dimension(nlon,nlat,nlev),intent(in) :: rhs,boundaries,feta
    real,dimension(nlon,nlat,nlev),intent(inout) :: omegaold,omega,resid
    logical,intent(in) :: lres
    real,intent(in) :: sigma0(nlev),dx,dy,dlev,alfa
    integer :: i,j,k

    do j=1,nlat       
       do i=1,nlon         
          omegaold(i,j,1)=boundaries(i,j,1)
          omegaold(i,j,nlev)=boundaries(i,j,nlev)
       enddo
    enddo
    do k=2,nlev-1
       do i=1,nlon
          omegaold(i,1,k)=boundaries(i,1,k)
          omegaold(i,nlat,k)=boundaries(i,nlat,k)
       enddo
    enddo
    
    omega=omegaold

    do i=1,niter
       call updateQG(omegaold,omega,sigma0,feta,rhs,dx,dy,dlev,alfa)
    enddo
    
    if(lres)then
       call residQG(rhs,omega,sigma0,feta,dx,dy,dlev,resid)
    endif
    
  end subroutine solveQG

  subroutine updateQG(omegaold,omega,sigma,etasq,rhs,dx,dy,dlev,alfa)
!
!   New estimate for the local value of omega, using omega in the 
!   surrounding points and the right-hand-side forcing (rhs)
!
!   QG version: for 'sigma' and 'etasq', constant values from the QG theory
!   are used.
! 
    implicit none

    real,dimension(:,:,:),intent(in) :: rhs,etasq
    real,dimension(:,:,:),intent(inout) :: omegaold,omega
    real,dimension(:),    intent(in) :: sigma
    real,                 intent(in) :: dx,dy,dlev,alfa
    
    real :: maxdiff
    integer :: i,j,k,nlon,nlat,nlev
    real,dimension(:,:,:),allocatable :: lapl2,coeff1,coeff2,coeff,domedp2

    nlon=size(rhs,1)
    nlat=size(rhs,2)
    nlev=size(rhs,3)
    allocate(lapl2(nlon,nlat,nlev),coeff1(nlon,nlat,nlev))
    allocate(coeff2(nlon,nlat,nlev),coeff(nlon,nlat,nlev))
    allocate(domedp2(nlon,nlat,nlev))
    
!   Top and bottom levels: omega directly from the boundary conditions,
!   does not need to be solved.
!
    call laplace2_cart(omegaold,dx,dy,lapl2,coeff1)
    call p2der2(omegaold,dlev,domedp2,coeff2) 

!    write(*,*)'Calculate the coefficients'
!    write(*,*)nlon,nlat,nlev
!    write(*,*)'coeff(nlon,nlat,nlev-1)',coeff(nlon,nlat,nlev-1)
!    write(*,*)'coeff1(nlon,nlat,nlev-1)',coeff1(nlon,nlat,nlev-1)
!    write(*,*)'coeff2(nlon,nlat,nlev-1)',coeff2(nlon,nlat,nlev-1)
!    write(*,*)'sigma(nlon,nlat,nlev-1)',sigma(nlon,nlat,nlev-1)
!    write(*,*)'domedp2(nlon,nlat,nlev-1)',domedp2(nlon,nlat,nlev-1)
!    write(*,*)'lapl2(nlon,nlat,nlev-1)',lapl2(nlon,nlat,nlev-1)
!    write(*,*)'etasq(nlon,nlat,nlev-1)',etasq(nlon,nlat,nlev-1)

    do k=2,nlev-1
       do j=2,nlat-1
          do i=1,nlon
             coeff(i,j,k)=sigma(k)*coeff1(i,j,k)+etasq(i,j,k)*coeff2(i,j,k)
             omega(i,j,k)=(rhs(i,j,k)-sigma(k)*lapl2(i,j,k)-etasq(i,j,k)*domedp2(i,j,k)) &
                  /coeff(i,j,k) 
          enddo
       enddo
    enddo
    
!   write(*,*)'Updating omega'
    maxdiff=0.
    do k=2,nlev-1
       do j=2,nlat-1
          do i=1,nlon
             maxdiff=max(maxdiff,abs(omega(i,j,k)-omegaold(i,j,k)))
             omegaold(i,j,k)=alfa*omega(i,j,k)+(1-alfa)*omegaold(i,j,k)
          enddo
       enddo
    enddo
    
  end subroutine updateQG


  subroutine residQG(rhs,omega,sigma,etasq,dx,dy,dlev,resid)
!
!   Calculating the residual RHS - LQG(omega)
!      
!    Variables:
!
!    omega = approximation for omega
!    sigma = local values of sigma (*after modifying for ellipticity*)
!    feta = f*eta (*after modifying for ellipticity*)
!    f = coriolis parameter
!    d2zetadp = second pressure derivative of relative vorticity 
!    dudp,dvdp = pressure derivatives of wind components
!    rhs = right-hand-side forcing
!
    implicit none

    real,dimension(:,:,:),intent(in) :: rhs,omega,etasq
    real,dimension(:,:,:),intent(out) :: resid
    real,dimension(:),    intent(in) :: sigma
    real,                 intent(in) :: dx,dy,dlev
    real,dimension(:,:,:),allocatable :: laplome,domedp2
    integer :: i,j,k,nlon,nlat,nlev

    nlon=size(rhs,1)
    nlat=size(rhs,2)
    nlev=size(rhs,3)

    allocate(laplome(nlon,nlat,nlev),domedp2(nlon,nlat,nlev))

    call laplace_cart(omega,laplome,dx,dy)         
    call p2der(omega,dlev,domedp2) 
 
    do k=1,nlev
       do j=1,nlat
          do i=1,nlon
             resid(i,j,k)=rhs(i,j,k)-(sigma(k)*laplome(i,j,k)+etasq(i,j,k)*domedp2(i,j,k))
             !if(i.eq.nlon/2.and.j.eq.nlat/2.and.k.eq.nlev/2)write(*,*)'rhs,resid',rhs(i,j,k),resid(i,j,k)
          enddo
       enddo
    enddo
    
  end subroutine residQG

  subroutine callsolvegen(rhs,boundaries,omega,nlon,nlat,nlev,&
       dx,dy,dlev,sigma0,sigma,feta,corpar,d2zetadp,dudp,dvdp,&
       nres,alfa,toler)
!
!      Calling solvegen + writing out omega. Multigrid algorithm
!            
    implicit none

    real,dimension(:,:,:,:),intent(inout) :: rhs,omega
    real,dimension(:,:,:,:),intent(in) :: boundaries,sigma,feta,d2zetadp
    real,dimension(:,:,:,:),intent(in) :: dudp,dvdp
    real,dimension(:,:),    intent(in) :: sigma0,corpar
    real,dimension(:),      intent(in) :: dx,dy,dlev
    real,                   intent(in) :: alfa,toler
    integer,dimension(:),   intent(in) :: nlon,nlat,nlev
    integer,                intent(in) :: nres

    real,dimension(:,:,:,:),allocatable :: omegaold
    real,dimension(:,:,:),allocatable :: dum1,resid,omega1

    real :: maxdiff,aomega
    integer :: ires,iter,i,j,k

    integer,parameter :: itermax=1000
    integer,parameter :: ny1=2,ny2=2 ! number of iterations at each grid 
                                     ! resolution when proceeding to coarser 
                                     ! (ny1) and when returning to finer (ny2)
    logical,parameter :: lzeromean=.true. ! Area means of omega are set to zero                    
 
    allocate(dum1(nlon(1),nlat(1),nlev(1)))
    allocate(resid(nlon(1),nlat(1),nlev(1)))
    allocate(omega1(nlon(1),nlat(1),nlev(1)))

    omega=0.
    omegaold=boundaries
!
!      This far: the whole multigrid cycle is written explicitly here  
!
!------------------------------------------------------------------------------------------
!
    do iter=1,itermax  ! Each iteration = one (fine->coarse->fine) multigrid cycle
!
!      Loop from finer to coarser resolutions
!
       do ires=1,nres
!         write(*,*)'fine-coarse:iter,ires',iter,ires
          call solvegen(rhs(:,:,:,ires),boundaries(:,:,:,ires),&
               omega(:,:,:,ires),omegaold(1,1,1,ires),nlon(ires),&
               nlat(ires),nlev(ires),dx(ires),dy(ires),dlev(ires),&
               sigma0(:,ires),sigma(:,:,:,ires),feta(:,:,:,ires),&
               corpar(:,ires),d2zetadp(:,:,:,ires),dudp(:,:,:,ires),&
               dvdp(:,:,:,ires),ny1,alfa,.true.,resid)
!             write(*,*)'ires,omega',ires,omega(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2,1)
!             write(*,*)'ires,resid',ires,resid(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2)
          if(ires.eq.1)omega1(:,:,:)=omega(:,:,:,1)
          if(ires.lt.nres)then
            call coarsen3d(resid,rhs(:,:,:,ires+1),nlon(ires),nlat(ires),&
                 nlev(ires),nlon(ires+1),nlat(ires+1),nlev(ires+1))          
!             write(*,*)'ires,rhs',ires,rhs(nlon(ires+1)/2,nlat(ires+1)/2,nlev(ires+1)/2,ires)
         endif  
       enddo         
!
!      Loop from coarser to finer resolutions
!
       do ires=nres-1,1,-1
!         write(*,*)'coarse-fine:iter,ires',iter,ires
          call finen3D(omega(:,:,:,ires+1),dum1,nlon(ires),nlat(ires),&
               nlev(ires),nlon(ires+1),nlat(ires+1),nlev(ires+1))          
!      Without the underrelaxation (coeffient alfa) the soultion diverges
          omegaold(:,:,:,ires)=omega(:,:,:,ires)+alfa*dum1(:,:,:)
!         if(ires.eq.1)then
!         write(*,*)'omega',ires,omega(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2,ires)
!         write(*,*)'dum1',ires,dum1(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2)
!         write(*,*)'omegaold',ires,omegaold(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2,ires)
!         endif

          call solvegen(rhs(:,:,:,ires),boundaries(:,:,:,ires),& 
               omega(:,:,:,ires),omegaold(:,:,:,ires),nlon(ires),nlat(ires),&
               nlev(ires),dx(ires),dy(ires),dlev(ires),sigma0(:,ires),&
               sigma(:,:,:,ires),feta(:,:,:,ires),corpar(:,ires),&
               d2zetadp(:,:,:,ires),dudp(:,:,:,ires),dvdp(:,:,:,ires),ny2,alfa,&
               .false.,resid)
       enddo

       maxdiff=0.
       do k=1,nlev(1)  
          do j=1,nlat(1)  
             do i=1,nlon(1)  
                maxdiff=max(maxdiff,abs(omega(i,j,k,1)-omega1(i,j,k)))
             enddo
          enddo
       enddo
!       write(*,*)iter,maxdiff
       if(maxdiff.lt.toler.or.iter.eq.itermax)then
!          write(*,*)'iter,maxdiff',iter,maxdiff
          goto 10
       endif

       omegaold=omega
          
    enddo ! iter=1,itermax
10  continue
!----------------------------------------------------------------------------------

   !       Subtracting area mean of omega
    if(lzeromean)then
       do k=1,nlev(1) 
          call aave(omega(:,:,k,1),aomega)
          do j=1,nlat(1)
             do i=1,nlon(1)
                omega(i,j,k,1)=omega(i,j,k,1)-aomega
             enddo
          enddo
       enddo
    endif
    
  end subroutine callsolvegen

  subroutine solvegen(rhs,boundaries,omega,omegaold,nlon,nlat,nlev,&
       dx,dy,dlev,sigma0,sigma,feta,corpar,d2zetadp,dudp,dvdp,&
       niter,alfa,lres,resid)
!
!      Solving omega iteratively using the generalized LHS operator.
!      'niter' iterations with relaxation coefficient alfa
!               
!      Input:
!
!      rhs = right-hand-side forcing
!      boundaries = boundary conditions 
!      omegaold,omega = old and new omega
!      sigma0 = area means of sigma at each pressure level
!      sigma = local values of sigma (*after modifying for ellipticity*)
!      feta = f*eta (*after modifying for ellipticity*)
!      f = coriolis parameter
!      d2zetadp = second pressure derivative of relative vorticity 
!      dudp,dvdp = pressure derivatives of wind components
!      rhs = right-hand-side forcing
!
!      output:
!
!      omega 
!      resid (if (lres))
!
    implicit none

    integer,intent(in) :: nlon,nlat,nlev,niter
    real,dimension(nlon,nlat,nlev),intent(inout) :: omega,omegaold
    real,dimension(nlon,nlat,nlev),intent(in) :: rhs,feta,sigma,boundaries
    real,dimension(nlon,nlat,nlev),intent(in) :: d2zetadp,dudp,dvdp
    real,dimension(nlon,nlat,nlev),intent(out) :: resid
    real,dimension(nlev),intent(in) :: sigma0
    real,dimension(nlat),intent(in) :: corpar
    real,                intent(in) :: dx,dy,dlev,alfa
    logical,             intent(in) :: lres
    integer :: i,j,k

    do j=1,nlat       
       do i=1,nlon
          omegaold(i,j,1)=boundaries(i,j,1)
          omegaold(i,j,nlev)=boundaries(i,j,nlev)
       enddo
    enddo
    do k=2,nlev-1
       do i=1,nlon
          omegaold(i,1,k)=boundaries(i,1,k)
          omegaold(i,nlat,k)=boundaries(i,nlat,k)
       enddo
    enddo

    omega=omegaold

!       write(*,*)'Boundary conditions given'

    do i=1,niter
       call updategen(omegaold,omega,sigma0,sigma,feta,corpar,d2zetadp,dudp,&
            dvdp,rhs,dx,dy,dlev,alfa)
    enddo
!
!      Calculate the residual = RHS - L(omega)

    if(lres)then
       call residgen(rhs,omega,resid,sigma,feta,corpar,d2zetadp,dudp,dvdp,&
            dx,dy,dlev)
    endif

  end subroutine solvegen
  
  subroutine updategen(omegaold,omega,sigma0,sigma,feta,f,d2zetadp,dudp,dvdp,&
       rhs,dx,dy,dlev,alfa)
!
!      Calculating new local values of omega, based on omega in the
!      surrounding points and the right-hand-side forcing (rhs).
!      
!      The left-hand-side of the omega equation as in (Räisänen 1995),
!      but terms reorganised following Pauley & Nieman (1992)
!         
!      Variables:
!
!      omegaold,omega = old and new omega
!      sigma0 = area means of sigma at each pressure level
!      sigma = local values of sigma (*after modifying for ellipticity*)
!      feta = f*eta (*after modifying for ellipticity*)
!      f = coriolis parameter
!      d2zetadp = second pressure derivative of relative vorticity 
!      dudp,dvdp = pressure derivatives of wind components
!      rhs = right-hand-side forcing
!
    implicit none

    integer i,j,k,nlon,nlat,nlev
    real,dimension(:,:,:),intent(inout) :: omegaold,omega 
    real,dimension(:,:,:),intent(in) :: sigma,feta,rhs,d2zetadp,dudp,dvdp
    real,dimension(:),    intent(in) :: sigma0,f
    real,                 intent(in) :: dx,dy,dlev,alfa
    
    real,dimension(:,:,:),allocatable :: lapl2,domedp2,coeff1,coeff2,coeff
    real,dimension(:,:,:),allocatable :: dum0,dum1,dum3,dum4,dum5,dum6

    nlon=size(rhs,1)
    nlat=size(rhs,2)
    nlev=size(rhs,3)
       
    allocate(lapl2(nlon,nlat,nlev),domedp2(nlon,nlat,nlev))
    allocate(coeff1(nlon,nlat,nlev),coeff2(nlon,nlat,nlev))
    allocate(coeff(nlon,nlat,nlev),dum0(nlon,nlat,nlev))
    allocate(dum1(nlon,nlat,nlev),dum3(nlon,nlat,nlev))
    allocate(dum4(nlon,nlat,nlev),dum5(nlon,nlat,nlev))
    allocate(dum6(nlon,nlat,nlev))
!
!   Top and bottom levels: omega directly from the boundary conditions,
!   does not need to be solved.
!
    call laplace2_cart(omegaold,dx,dy,lapl2,coeff1)
    call p2der2(omegaold,dlev,domedp2,coeff2) 
!
!   Calculate non-constant terms on the left-hand-side, based on 'omegaold'
!
!   a) Deviation of sigma from its normal value

    do k=2,nlev-1
       dum0(:,:,k)=omegaold(:,:,k)*(sigma(:,:,k)-sigma0(k))
    enddo
    call laplace_cart(dum0,dum1,dx,dy)         
!
!   b) f*omega*(d2zetadp): explicitly, later
!          
!   c) tilting
!
    call xder_cart(omegaold,dx,dum4) 
    call yder_cart(omegaold,dy,dum5) 
 
    do j=1,nlat
       dum6(:,j,:)=f(j)*(dudp(:,j,:)*dum5(:,j,:)-dvdp(:,j,:)*dum4(:,j,:))
    enddo
    call pder(dum6,dlev,dum3) 
!
!   Solving for omega 
!   Old values are retained at y and z boundaries.
!       

    do k=2,nlev-1
       do j=2,nlat-1
          coeff(:,j,k)=sigma0(k)*coeff1(:,j,k)+feta(:,j,k)*coeff2(:,j,k)-&
                          f(j)*d2zetadp(:,j,k)
          omega(:,j,k)=(rhs(:,j,k)-dum1(:,j,k)-dum3(:,j,k)-&
               sigma0(k)*lapl2(:,j,k)-feta(:,j,k)*domedp2(:,j,k)) &
               /coeff(:,j,k)
       enddo
    enddo
     
!    write(*,*)'updating omega'
    do k=2,nlev-1
       do j=2,nlat-1
          do i=1,nlon
             omegaold(i,j,k)=alfa*omega(i,j,k)+(1-alfa)*omegaold(i,j,k)
          enddo
       enddo
    enddo

  end subroutine updategen

  subroutine residgen(rhs,omega,resid,sigma,feta,f,d2zetadp,dudp,dvdp, &
       dx,dy,dlev)
!
!   Calculating the residual RHS - L(omega)
!      
!   Variables:
!
!   omega = approximation for omega
!   sigma = local values of sigma (*after modifying for ellipticity*)
!   feta = f*eta (*after modifying for ellipticity*)
!   f = coriolis parameter
!   d2zetadp = second pressure derivative of relative vorticity 
!   dudp,dvdp = pressure derivatives of wind components
!   rhs = right-hand-side forcing
!
    implicit none

    real,dimension(:,:,:),intent(in) :: rhs,omega,sigma,feta,d2zetadp
    real,dimension(:,:,:),intent(in) :: dudp,dvdp
    real,dimension(:,:,:),intent(out) :: resid
    real,dimension(:),intent(in) :: f 
    real,intent(in) :: dx,dy,dlev
    integer :: i,j,k,nlon,nlat,nlev
    real,dimension(:,:,:),allocatable :: dum0,dum1,dum2,dum3,dum4,dum5,dum6
    
    nlon=size(rhs,1)
    nlat=size(rhs,2)
    nlev=size(rhs,3)

    allocate(dum0(nlon,nlat,nlev),dum2(nlon,nlat,nlev))
    allocate(dum1(nlon,nlat,nlev),dum3(nlon,nlat,nlev))
    allocate(dum4(nlon,nlat,nlev),dum5(nlon,nlat,nlev))
    allocate(dum6(nlon,nlat,nlev))    
!
!   Calculate L(omega)

!    a) nabla^2(sigma*omega)

    dum0=omega*sigma

    call laplace_cart(dum0,dum1,dx,dy)         
!
!   f*eta*d2omegadp
!       
    call p2der(omega,dlev,dum2) 
!
    dum3=feta*dum2
!
!   c) -f*omega*(d2zetadp): explicitly, later
!                  
!   d) tilting
!
    call xder_cart(omega,dx,dum4) 
    call yder_cart(omega,dy,dum5) 
 
    do k=1,nlev
       do j=1,nlat
          do i=1,nlon
             dum6(i,j,k)=f(j)*(dudp(i,j,k)*dum5(i,j,k)-dvdp(i,j,k)*dum4(i,j,k))
          enddo
       enddo
    enddo
    call pder(dum6,dlev,dum2) 

    do k=1,nlev
       do j=1,nlat
          do i=1,nlon
             resid(i,j,k)=rhs(i,j,k)-(dum1(i,j,k)+dum2(i,j,k)+dum3(i,j,k)- &
                  f(j)*d2zetadp(i,j,k)*omega(i,j,k))
          enddo
       enddo
    enddo
    
  end subroutine residgen

  subroutine laplace2_cart(f,dx,dy,lapl2,coeff)
!
!      As laplace_cart but
!        - the contribution of the local value to the Laplacian is left out
!        - coeff is the coefficient for the local value
!
    real,dimension(:,:,:),intent(in) :: f
    real,dimension(:,:,:),intent(out) :: lapl2,coeff
    real,intent(in) :: dx,dy
    integer :: nlon,nlat,nlev,i,j,k,c
    
    nlon=size(f,1)
    nlat=size(f,2)
    nlev=size(f,3)
    c=1
    select case (c)
    case(1)
       
       ! x-direction
       lapl2 ( 2 : nlon - 1, :, : ) = f( 1: nlon - 2, :, : ) + f ( 3: nlon, :, : ) 
       lapl2 ( 1, :, : )    = f( nlon, :, : ) + f ( 2, :, : ) 
       lapl2 ( nlon, :, : ) = f( nlon - 1, :, : ) + f ( 1, :, : ) 
       lapl2 = lapl2 / ( dx * dx )
       
       ! y-direction
       lapl2 ( :, 2 : nlat -1, : ) = lapl2 ( :, 2 : nlat -1, : ) &
            + ( f ( :, 1 : nlat -2, : ) + f ( :, 3 : nlat, :) ) / ( dy * dy )
       
       coeff ( :, 2 : nlat -1, : ) = -2 / (dy**2.) -2 / (dx**2.)
       coeff ( :, 1, : ) = -2 / ( dx**2. )
       coeff ( :, nlat, : ) = -2 / ( dx**2. )
    case(2)
       do j=1,nlat
          do k=1,nlev
             do i=1,nlon
                ! x-direction
                if(i==1)then
                   lapl2(i,j,k)=(-(1./12.)*f(nlon-1,j,k)+(4./3.)*f(nlon,j,k) &
                        +(4./3.)*f(i+1,j,k)-(1./12.)*f(i+2,j,k))/(dx*dx)
                else if(i==2)then
                   lapl2(i,j,k)=(-(1./12.)*f(nlon,j,k)+(4./3.)*f(i-1,j,k) &
                        +(4./3.)*f(i+1,j,k)-(1./12.)*f(i+2,j,k))/(dx*dx)
                else if(i==nlon-1)then
                   lapl2(i,j,k)=(-(1./12.)*f(i-2,j,k)+(4./3.)*f(i-1,j,k) &
                        +(4./3.)*f(i+1,j,k)-(1./12.)*f(1,j,k))/(dx*dx)
                else if(i==nlon)then
                   lapl2(i,j,k)=(-(1./12.)*f(i-2,j,k)+(4./3.)*f(i-1,j,k) &
                        +(4./3.)*f(1,j,k)-(1./12.)*f(2,j,k))/(dx*dx)
                else
                   lapl2(i,j,k)=-(1./12.)*f(i-2,j,k)+(4./3.)*f(i-1,j,k) &
                        +(4./3.)*f(i+1,j,k)-(1./12.)*f(i+2,j,k)
                   lapl2(i,j,k)=lapl2(i,j,k)/(dx*dx)
                end if
             enddo
          enddo
       enddo
       
       do i=1,nlon
          do k=1,nlev
             do j=2,nlat-1
                ! y-direction
                if(j==2)then
                   lapl2(i,j,k)=lapl2(i,j,k)+(f(i,j-1,k)+f(i,j+1,k))/(dy*dy)
                else if(j==nlat-1)then
                   lapl2(i,j,k)=lapl2(i,j,k)+(f(i,j-1,k)+f(i,j+1,k))/(dy*dy)
                else
                   lapl2(i,j,k)=lapl2(i,j,k)+(-(1./12.)*f(i,j-2,k)+(4./3.)*f(i,j-1,k) &
                        +(4./3.)*f(i,j+1,k)-(1./12.)*f(i,j+2,k))/(dy*dy)
                end if
             enddo
          enddo
       enddo
       coeff(:,3:nlat-2,:)=-(5./2.)/(dx*dx)-(5./2.)/(dy*dy)
       coeff(:,2,:)=-(5./2.)/(dx*dx)-2./(dy*dy)
       coeff(:,nlat-1,:)=-(5./2.)/(dx*dx)-2./(dy*dy)
       coeff(:,1,:)=-(5./2.)/(dx*dx)
       coeff(:,nlat,:)=-(5./2.)/(dx*dx)
    end select

  end subroutine laplace2_cart

  subroutine p2der(f,dp,df2dp2) 
!
!      Estimation of second pressure derivatives.
!      At top and bottom levels, these are set to zero
!
    implicit none
    real,dimension(:,:,:),intent(in) :: f
    real,dimension(:,:,:),intent(out) :: df2dp2
    real,intent(in) :: dp
    integer :: nlev,c,k

    nlev=size(f,3)
    c=1
    select case(c)
       
    case(1)
       df2dp2(:,:,2:nlev-1)=(f(:,:,3:nlev)+f(:,:,1:nlev-2) &
            -2*f(:,:,2:nlev-1))/(dp*dp)
       df2dp2(:,:,1)=0.
       df2dp2(:,:,nlev)=0.
    case(2)
       do k=3,nlev-2
          df2dp2(:,:,k)=(-1./12.)*f(:,:,k-2)+(4./3.)*f(:,:,k-1)&
               -(5./2.)*f(:,:,k)+(4./3.)*f(:,:,k+1)-(1./12.)*f(:,:,k+2)
          df2dp2(:,:,k)=df2dp2(:,:,k)/(dp*dp)
       enddo
       df2dp2(:,:,2)=(f(:,:,3)+f(:,:,1)-2.*f(:,:,2))/(dp*dp)
       df2dp2(:,:,nlev-1)=(f(:,:,nlev)+f(:,:,nlev-2)-2.*f(:,:,nlev-1))/(dp*dp)
       df2dp2(:,:,1)=0.
       df2dp2(:,:,nlev)=0.
    end select
 end subroutine p2der

 subroutine p2der2(f,dp,df2dp22,coeff) 
!
!      As p2der, but
!        - the contribution of the local value is left out
!        - the coefficient 'coeff' of the local value is also calculated        
!
   implicit none
   real,dimension(:,:,:),intent(in) :: f
   real,dimension(:,:,:),intent(out) :: df2dp22,coeff
   real,intent(in) :: dp
   integer :: nlev,c,k

   c=1
   nlev=size(f,3)
 
   select case (c)
   case(1)
      df2dp22(:,:,2:nlev-1)=(f(:,:,3:nlev)+f(:,:,1:nlev-2)) &
           /(dp*dp)
      coeff(:,:,2:nlev-1)=-2./(dp*dp)
      df2dp22(:,:,1)=0.
      df2dp22(:,:,nlev)=0.
      coeff(:,:,1)=0.
      coeff(:,:,nlev)=0.
   case(2)
      do k=3,nlev-2
         df2dp22(:,:,k)=(-1./12.)*f(:,:,k-2)+(4./3.)*f(:,:,k-1) &
              +(4./3.)*f(:,:,k+1)-(1./12.)*f(:,:,k+2)
         df2dp22(:,:,k)=df2dp22(:,:,k)/(dp*dp)
      enddo
      coeff(:,:,3:nlev-2)=(-5./2.)/(dp*dp)
      df2dp22(:,:,2)=(f(:,:,3)+f(:,:,1))/(dp*dp)
      df2dp22(:,:,nlev-1)=(f(:,:,nlev)+f(:,:,nlev-2))/(dp*dp)
      coeff(:,:,2)=-2./(dp*dp)
      coeff(:,:,nlev-1)=-2./(dp*dp)
      coeff(:,:,1)=0
      coeff(:,:,nlev)=0
      df2dp22(:,:,1)=0.
      df2dp22(:,:,nlev)=0.
   end select
      
 end subroutine p2der2


!************ SUBROUTINES NOT CURRENTLY USED ****************************************

       subroutine coarsen3dbil(f,g,nlon1,nlat1,nlev1,nlon2,nlat2,nlev2)
!
!      DOES NOT WORK, DON'T KNOW WHY (JR 230316)
!
!      Bilinear interpolation from grid 1 (nlon1,nlat1,nlev1)
!      to grid 2 (nlon2,nlat2,nlev2). Assumptions:
!       - all grids have even spacing in all directions
!       - the boundaries are located at x/y/z=0.5 and x/y/z=nlon1-2/nlat1-2/nlev1-2+0.5
!       - the domain is periodic is in the x direction but not in y and z
!      With these assumptions, only the grid sizes are needed as input,
!      not the actual coordinates 
!
       implicit none
       integer nlon1,nlat1,nlon2,nlat2,nlev1,nlev2
       real f(nlon1,nlat1,nlev1),g(nlon2,nlat2,nlev2)
       integer i,j,k
       real rx,ry,rz,rlon,rlat,rlev
       integer maxsiz
       parameter (maxsiz=1000)
       integer ii1(maxsiz),jj1(maxsiz),kk1(maxsiz),ii2(maxsiz),jj2(maxsiz),kk2(maxsiz)       
       real ai(maxsiz),aj(maxsiz),ak(maxsiz),bi(maxsiz),bj(maxsiz),bk(maxsiz)       

       rlon=real(nlon1)/real(nlon2) 
       rlat=real(nlat1)/real(nlat2) 
       rlev=real(nlev1)/real(nlev2) 
!
!      Tabulate coordinates and coefficients to speed up
!
       do i=1,nlon2          
         rx=0.5+rlon*(i-0.5)
         ii1(i)=int(rx)
         bi(i)=rx-ii1(i)
         if(ii1(i).lt.1)ii1(i)=nlon1
         ii2(i)=ii1(i)+1
         if(ii2(i).gt.nlon1)ii2(i)=1
         ai(i)=1-bi(i)

!         write(*,*)'Coarsen3dbil: i,ii1,ii2,ai,bi',i,ii1(i),ii2(i),ai(i),bi(i)
       enddo 

       do j=1,nlat2          
         ry=0.5+rlat*(j-0.5)
         jj1(j)=int(ry)
         jj2(j)=jj1(j)+1
         bj(j)=ry-jj1(j)
         if(jj1(j).lt.1)then
           jj1(j)=1
           bj(j)=0.
         endif
         if(jj2(j).gt.nlat1)then
           jj2(j)=nlat1
           bj(j)=1.
         endif
         aj(j)=1-bj(j)
!         write(*,*)'Coarsen3dbil: j,jj1,jj2,aj,bj',j,jj1(j),jj2(j),aj(j),bj(j)
       enddo 

       do k=1,nlev2          
         rz=0.5+rlev*(k-0.5)
         kk1(k)=int(rz)
         kk2(k)=kk1(k)+1
         bk(k)=rz-kk1(k)
         if(kk1(k).lt.1)then
           kk1(k)=1
           bk(k)=0.
         endif
         if(kk2(k).gt.nlev1)then
           kk2(k)=nlev1
         endif
         ak(k)=1-bk(k)

!         write(*,*)'Coarsen3dbil: k,kk1,kk2,ak,bk',k,kk1(k),kk2(k),ak(k),bk(k)
       enddo 

       do k=1,nlev2
       do j=1,nlat2
       do i=1,nlon2
            g(i,j,k)=ai(i)*aj(j)*ak(k)*f(ii1(i),jj1(j),kk1(k))+&
                     bi(i)*aj(j)*ak(k)*f(ii2(i),jj1(j),kk1(k))+&
                     ai(i)*bj(j)*ak(k)*f(ii1(i),jj2(j),kk1(k))+&
                     bi(i)*bj(j)*ak(k)*f(ii2(i),jj2(j),kk1(k))+&
                     ai(i)*aj(j)*bk(k)*f(ii1(i),jj1(j),kk2(k))+&
                     bi(i)*aj(j)*bk(k)*f(ii2(i),jj1(j),kk2(k))+&
                     ai(i)*bj(j)*bk(k)*f(ii1(i),jj2(j),kk2(k))+&
                     bi(i)*bj(j)*bk(k)*f(ii2(i),jj2(j),kk2(k))
       enddo
       enddo
       enddo

       return
       end subroutine coarsen3dbil

       subroutine finen3Dbil(f,g,nlon2,nlat2,nlev2,nlon1,nlat1,nlev1)
!
!      Bilinear interpolation from grid 1 (nlon1,nlat1,nlev1)
!      to grid 2 (nlon2,nlat2,nlev2). Assumptions:
!       - all grids have even spacing in all directions
!       - the boundaries are located at x/y/z=0.5 and x/y/z=nlon1-2/nlat1-2/nlev1-2+0.5
!       - the domain is periodic is in the x direction but not in y and z
!      With these assumptions, only the grid sizes are needed as input,
!      not the actual coordinates 
!
!      SAME AS COARSEN3DBIL, except for the order of the arguments
!
       implicit none
       integer nlon1,nlat1,nlon2,nlat2,nlev1,nlev2
       real f(nlon1,nlat1,nlev1),g(nlon2,nlat2,nlev2)
       integer i,j,k
       real rx,ry,rz,rlon,rlat,rlev
       integer maxsiz
       parameter (maxsiz=1000)
       integer ii1(maxsiz),jj1(maxsiz),kk1(maxsiz),ii2(maxsiz),jj2(maxsiz),kk2(maxsiz)       
       real ai(maxsiz),aj(maxsiz),ak(maxsiz),bi(maxsiz),bj(maxsiz),bk(maxsiz)       

       rlon=real(nlon1)/real(nlon2) 
       rlat=real(nlat1)/real(nlat2) 
       rlev=real(nlev1)/real(nlev2) 
!
!      Tabulate coordinates and coefficients to speed up
!
       do i=1,nlon2          
         rx=0.5+rlon*(i-0.5)
         ii1(i)=int(rx)
         bi(i)=rx-ii1(i)
         ai(i)=1-bi(i)
         if(ii1(i).lt.1)ii1(i)=nlon1
         ii2(i)=ii1(i)+1
         if(ii2(i).gt.nlon1)ii2(i)=1
!         write(*,*)'Finen3D: i,ii1,ii2,ai,bi',i,ii1(i),ii2(i),ai(i),bi(i)
       enddo 

       do j=1,nlat2          
         ry=0.5+rlat*(j-0.5)
         jj1(j)=int(ry)
         jj2(j)=jj1(j)+1
         bj(j)=ry-jj1(j)
         if(jj1(j).lt.1)then
           jj1(j)=1
           bj(j)=0.
         endif
         if(jj2(j).gt.nlat1)then
           jj2(j)=nlat1
           bj(j)=1.
         endif
         aj(j)=1-bj(j)
!         write(*,*)'Finen3Dbil: j,jj1,jj2,aj,bj',j,jj1(j),jj2(j),aj(j),bj(j)
       enddo 

       do k=1,nlev2          
         rz=0.5+rlev*(k-0.5)
         kk1(k)=int(rz)
         kk2(k)=kk1(k)+1
         bk(k)=rz-kk1(k)
         if(kk1(k).lt.1)then
           kk1(k)=1
           bk(k)=0.
         endif
         if(kk2(k).gt.nlev1)then
           kk2(k)=nlev1
           bk(k)=1.
         endif
         ak(k)=1-bk(k)

!         write(*,*)'Finen3Dbil: k,kk1,kk2,ak,bk',k,kk1(k),kk2(k),ak(k),bk(k)
       enddo 
      
       do k=1,nlev2
       do j=1,nlat2
       do i=1,nlon2
            g(i,j,k)=ai(i)*aj(j)*ak(k)*f(ii1(i),jj1(j),kk1(k))+&
                     bi(i)*aj(j)*ak(k)*f(ii2(i),jj1(j),kk1(k))+&
                     ai(i)*bj(j)*ak(k)*f(ii1(i),jj2(j),kk1(k))+&
                     bi(i)*bj(j)*ak(k)*f(ii2(i),jj2(j),kk1(k))+&
                     ai(i)*aj(j)*bk(k)*f(ii1(i),jj1(j),kk2(k))+&
                     bi(i)*aj(j)*bk(k)*f(ii2(i),jj1(j),kk2(k))+&
                     ai(i)*bj(j)*bk(k)*f(ii1(i),jj2(j),kk2(k))+&
                     bi(i)*bj(j)*bk(k)*f(ii2(i),jj2(j),kk2(k))
       enddo
       enddo
       enddo

       return
       end subroutine finen3Dbil

     end module mod_omega_subrs

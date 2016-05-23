module mod_omega
  use mod_wrf_file
  use mod_common_subrs
  use mod_omega_subrs
  implicit none

contains

  subroutine calculate_omegas( t, u, v, w, z, lev, dx, dy, corpar, q, &
       xfrict, yfrict, utend, vtend, ttend, psfc, omegas)

    real,dimension(:,:,:,:),intent(inout) :: omegas
    real,dimension(:,:,:),intent(in) :: t,w,xfrict,yfrict
    real,dimension(:,:,:),intent(in) :: utend,vtend
    real,dimension(:,:),  intent(in) :: psfc
    real,dimension(:),intent(in) :: lev,corpar
    real,                 intent(in) :: dx,dy 
    real,dimension(:,:,:),intent(inout) :: z,q,u,v,ttend

    real,dimension(:,:,:,:),allocatable :: omega,omegaold,boundaries
    real,dimension(:,:,:),allocatable :: sigmaraw,omegaan,zetaraw,dum6
    real,dimension(:,:,:),allocatable :: dum4,dum5,dum2,dum1,lapl,dum3
    real,dimension(:,:,:),allocatable :: dudp,dvdp,zetatend,mulfact
    real,dimension(:,:,:),allocatable :: zero,feta
    real,dimension(:,:,:),allocatable :: sigma,domedp2
    real,dimension(:,:,:),allocatable :: df2dp2,zeta,d2zetadp
    real,dimension(:,:,:,:),allocatable :: fv,ft,ff,fq,fa,boundaries2
    real,dimension(:,:,:,:),allocatable :: zero2,sigma2,feta2,d2zetadp2
    real,dimension(:,:,:,:),allocatable :: dudp2,dvdp2,ftest
    real,dimension(:,:,:),allocatable :: dum0,resid,omega1
    real,dimension(:,:),allocatable :: corpar2,sigma02
    real,dimension(:),allocatable :: sigma0
    integer,dimension(:),allocatable :: nlonx,nlatx,nlevx
    real,dimension(:),allocatable :: dx2,dy2,dlev2

    integer :: nlev,nlon,nlat,nres,i,j,k
    real :: dlev
    logical :: lcensor
    character :: mode

    ! itermax = maximum number of multigrid cycles
    ! ny1 = number of iterations for each grid when going to coarser grids
    ! ny2 = number of iterations for each grid when returning to finer grids
    !
    ! Threshold values to keep the generalized omega equation elliptic.
    real,parameter :: sigmamin=2e-7,etamin=2e-6

!      For iubound, ilbound and iybound are 0, horizontal boundary
!      conditions are used at the upper, lower and north/south boundaries  
!      A value of 1 for any of these parameters means that the boundary
!      condition is taken directly from the "real" WRF omega. In practice,
!      only the lower boundary condition (ilbound) is important.
       integer iubound,ilbound,iybound

       iubound=1 ! 1 for "real" omega as upper-boundary condition
       ilbound=1 ! 1 for "real" omega as upper-boundary condition
       iybound=1 ! 1 for "real" omega as north/south boundary condtion
       mode='G'   ! Mode = generalized omega equation
       lcensor=.true. ! Forcing below surface is set to zero.

       nlon=size(t,1)
       nlat=size(t,2)
       nlev=size(t,3)

       if(mode.eq.'G')write(*,*)'Generalized omega equation'   
       if(mode.eq.'Q')write(*,*)'Quasi-geostrophic omega equation'   
       if(mode.eq.'T')write(*,*)'Generalized test version'   
       if(mode.eq.'t')write(*,*)'Quasigeostrophic test version'   
       if(mode.ne.'G'.and.mode.ne.'Q'.and.mode.ne.'T'.and.mode.ne.'t')then
         write(*,*)'Unknown mode of operation. Aborting'
         stop
       endif

       allocate(sigma0(nlev))
       allocate(sigmaraw(nlon,nlat,nlev),zetaraw(nlon,nlat,nlev))
       allocate(dum1(nlon,nlat,nlev),sigma(nlon,nlat,nlev))
       allocate(dum2(nlon,nlat,nlev))
       allocate(dum3(nlon,nlat,nlev),domedp2(nlon,nlat,nlev))
       allocate(dum4(nlon,nlat,nlev))
       allocate(dum5(nlon,nlat,nlev))
       allocate(lapl(nlon,nlat,nlev))
       allocate(dum6(nlon,nlat,nlev),df2dp2(nlon,nlat,nlev))
       allocate(dvdp(nlon,nlat,nlev),zeta(nlon,nlat,nlev))
       allocate(dudp(nlon,nlat,nlev),d2zetadp(nlon,nlat,nlev))
       allocate(zetatend(nlon,nlat,nlev),feta(nlon,nlat,nlev))
       allocate(mulfact(nlon,nlat,nlev),dum0(nlon,nlat,nlev))
       allocate(resid(nlon,nlat,nlev))
       allocate(omega1(nlon,nlat,nlev))
       allocate(zero(nlon,nlat,nlev),omegaan(nlon,nlat,nlev))
 
       omegaan=w
       zero=0.
!
!      Number of different resolutions in solving the equation = nres 
!      Choose so that the coarsest grid has at least 5 points
!
       nres=1+int(log(max(nlon,nlat,nlev)/5.)/log(2.))
       write(*,*)'nres=',nres
!
       allocate(omega(nlon,nlat,nlev,nres),omegaold(nlon,nlat,nlev,nres))
       allocate(boundaries(nlon,nlat,nlev,nres))
       allocate(fv(nlon,nlat,nlev,nres))
       allocate(ft(nlon,nlat,nlev,nres))
       allocate(ff(nlon,nlat,nlev,nres))
       allocate(fq(nlon,nlat,nlev,nres))
       allocate(fa(nlon,nlat,nlev,nres))
       allocate(boundaries2(nlon,nlat,nlev,nres))
       allocate(zero2(nlon,nlat,nlev,nres))
       allocate(sigma2(nlon,nlat,nlev,nres))
       allocate(feta2(nlon,nlat,nlev,nres))
       allocate(d2zetadp2(nlon,nlat,nlev,nres))
       allocate(dudp2(nlon,nlat,nlev,nres))
       allocate(dvdp2(nlon,nlat,nlev,nres))
       allocate(ftest(nlon,nlat,nlev,nres))
       allocate(corpar2(nlat,nres),sigma02(nlev,nres))
       allocate(nlonx(nres),nlatx(nres),nlevx(nres))
       allocate(dx2(nres),dy2(nres),dlev2(nres))
       dlev=lev(2)-lev(1)

!      Grid sizes for the different resolutions
!
       nlonx(1)=nlon
       nlatx(1)=nlat
       nlevx(1)=nlev
       dx2(1)=dx
       dy2(1)=dy
       dlev2(1)=dlev
       do i=2,nres
         nlonx(i)=max(nlonx(i-1)/2,5)
         nlatx(i)=max(nlatx(i-1)/2,5)
         nlevx(i)=max(nlevx(i-1)/2,5)
         dx2(i)=dx*real(nlon)/real(nlonx(i))
         dy2(i)=dy*real(nlat)/real(nlatx(i))
         dlev2(i)=dlev*real(nlev)/real(nlevx(i))
         write(*,*)'i,nlonx,nlatx,nlevx',i,nlonx(i),nlatx(i),nlevx(i)
         write(*,*)'i,dx2,dy2,dlev2',i,dx2(i),dy2(i),dlev2(i)
       enddo
!
!      For quasi-geostrophic equation: calculation of geostrophic winds
!
       if(mode.eq.'Q')then
         call gwinds(z,dx,dy,corpar,u,v)
       endif        
!
!      Multiplication factor for forcing: 
!      1 above the ground, smaller (or 0) below the ground
!
       if(lcensor)then         
         call calmul(psfc,lev,nlev,mulfact) 
       endif     
!
!      Calculation of vorticity and vorticity tendency 
!
       call curl_cart(u,v,dx,dy,zetaraw)
       if(mode.eq.'G')call curl_cart(utend,vtend,dx,dy,zetatend)
!
!      Calculation of forcing terms 
!
       if(mode.eq.'G'.or.mode.eq.'Q')then

       call fvort(u,v,zetaraw,corpar,dx,dy,dlev,mulfact,fv(:,:,:,1))
       call ftemp(u,v,t,lev,dx,dy,mulfact,ft(:,:,:,1))
       endif

       if(mode.eq.'G')then

       call ffrict(xfrict,yfrict,corpar,dx,dy,dlev,mulfact,ff(:,:,:,1))
       call fdiab(q,lev,dx,dy,mulfact,fq(:,:,:,1))
       call fimbal(zetatend,ttend,corpar,lev,dx,dy,dlev,mulfact,fa(:,:,:,1))
       endif 
!
!      Deriving quantities needed for the LHS of the 
!      QG and/or generalised omega equation.
!

!      1. Pressure derivatives of wind components

       call pder(u,dlev,dudp) 
       call pder(v,dlev,dvdp) 
!
!      2. Stability sigma
!
       call define_sigma(t,lev,dlev,sigmaraw)
!
!      3. Modifying stability and vorticity on the LHS to keep
!      the solution elliptic
!
       call modify(sigmaraw,sigmamin,etamin,zetaraw,&
                  corpar,dudp,dvdp,sigma,feta,zeta)      
!
!      4. Second pressure derivative of vorticity 
!
       call p2der(zeta,dlev,d2zetadp) 
!
!      5. Area mean of static stability over the whole grid
!
       do k=1,nlev
       call aave(sigmaraw(:,:,k),sigma0(k))                      
       enddo
       write(*,*)'Area mean static stability'
       do k=1,nlev
         write(*,*)lev(k),sigma0(k)
       enddo  
!
!      Left-hand side coefficients for the QG equation
!

       if(mode.eq.'Q'.or.mode.eq.'t')then
         do k=1,nlev
         do j=1,nlat       
         do i=1,nlon
           sigma(i,j,k)=sigma0(k)
           feta(i,j,k)=corpar(1)**2.
         enddo
         enddo
         enddo
       endif
!
!      Forcing for quasigeostrophic test case ('t')
!      In essence: calculating the LHS from the WRF omega (omegaan) 
!      and substituting it to the RHS
!
       if(mode.eq.'t')then

       call p2der(omegaan,dlev,df2dp2)
       call laplace_cart(omegaan,lapl,dx,dy)
!
       do k=1,nlev
       do j=1,nlat       
       do i=1,nlon
         ftest(i,j,k,1)=sigma(i,j,k)*lapl(i,j,k)+feta(i,j,k)*df2dp2(i,j,k) 
       enddo
       enddo
       enddo

       endif ! mode.eq.'t'   
!
!      Forcing for the general test case
!      In essence: calculating the LHS from the WRF omega (omegaan) 
!      and substituting it to the RHS
       
       if(mode.eq.'T')then

       do k=1,nlev
       do j=1,nlat
       do i=1,nlon
          dum2(i,j,k)=sigmaraw(i,j,k)*omegaan(i,j,k)
       enddo
       enddo
       enddo

       call laplace_cart(dum2,dum1,dx,dy)
       call p2der(omegaan,dlev,dum3)
       call p2der(zetaraw,dlev,dum5)

       do k=1,nlev
       do j=1,nlat       
       do i=1,nlon
         dum2(i,j,k)=(corpar(j)+zetaraw(i,j,k))*corpar(j)*dum3(i,j,k)
       enddo
       enddo
       enddo

       do k=1,nlev
       do j=1,nlat       
       do i=1,nlon
         dum3(i,j,k)=-corpar(j)*omegaan(i,j,k)*dum5(i,j,k)
       enddo
       enddo
       enddo

       call xder_cart(omegaan,dx,dum4) 
       call yder_cart(omegaan,dy,dum5) 
 
       do k=1,nlev
       do j=1,nlat
       do i=1,nlon
          dum6(i,j,k)=-corpar(j)*(dvdp(i,j,k)*dum4(i,j,k)-dudp(i,j,k)*dum5(i,j,k))
       enddo
       enddo
       enddo        
       call pder(dum6,dlev,dum4) 

       do k=1,nlev
       do j=1,nlat       
       do i=1,nlon
         ftest(i,j,k,1)=dum1(i,j,k)+dum2(i,j,k)+dum3(i,j,k)+dum4(i,j,k) 
       enddo
       enddo
       enddo

       endif ! (forcing for the general test case if mode.eq.'T')

!      Boundary conditions from WRF omega?  
!
       boundaries=0.
!
       do j=1,nlat       
       do i=1,nlon
         boundaries(i,j,1,1)=iubound*omegaan(i,j,1)
         boundaries(i,j,nlev,1)=ilbound*omegaan(i,j,nlev)
       enddo
       enddo      
       do k=2,nlev-1
       do i=1,nlon
         boundaries(i,1,k,1)=iybound*omegaan(i,1,k)
         boundaries(i,nlat,k,1)=iybound*omegaan(i,nlat,k)
       enddo
       enddo     

!
!      Regrid left-hand-side parameters and boundary conditions to 
!      coarser grids. Note that non-zero boundary conditions are only 
!      possibly given at the highest resolutions (As only the 
!      'residual omega' is solved at lower resolutions)

       do i=1,nres
         if(i.eq.1)then
           call coarsen3d(boundaries,boundaries2(:,:,:,i),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         else
           call coarsen3d(zero,boundaries2(:,:,:,i),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         endif 
         call coarsen3d(zero,zero2(:,:,:,i),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(sigma,sigma2(:,:,:,i),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(feta,feta2(:,:,:,i),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(d2zetadp,d2zetadp2(:,:,:,i),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(dudp,dudp2(:,:,:,i),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(dvdp,dvdp2(:,:,:,i),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(sigma0,sigma02(:,i),1,1,nlev,1,1,nlevx(i))
         call coarsen3d(corpar,corpar2(:,i),1,nlat,1,1,nlatx(i),1)
       enddo

! *************************************************************************
! ***** Solving for omega, using the forcing and the LHS coefficients *****    
! ***** and possibly boundary conditions **********************************
! *************************************************************************

!      1) Test cases ('T','t'): only one forcing ('ftest') is used, but
!         the results are written out for every resolution.
!      2) Other cases: the vertical motion associated with each individual
!         forcing term + boundary conditions is written out separately
!         (-> 2 + 1 = 3 terms for QG omega, 5 + 1 terms for generalized 
!         omega) 
!
!       iunit=1

       if(mode.eq.'T')then
         call callsolvegen(ftest,boundaries2,omega,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,&
              nres)
       endif 

       if(mode.eq.'t')then
         call callsolveQG(ftest,boundaries2,omega,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,feta2,nres)
       endif 

       if(mode.eq.'G')then            
         Write(*,*)'Vorticity advection'
         call callsolvegen(fv,zero2,omega,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,&
              nres)
         omegas(:,:,:,termV)=omega(:,:,:,1)

         Write(*,*)'Thermal advection'
         call callsolvegen(ft,zero2,omega,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,&
              nres)
         omegas(:,:,:,termT)=omega(:,:,:,1)

         Write(*,*)'Friction'
         call callsolvegen(ff,zero2,omega,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,&
              nres)
         omegas(:,:,:,termF)=omega(:,:,:,1)

         Write(*,*)'Diabatic heating'
         call callsolvegen(fq,zero2,omega,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,&
              nres)
         omegas(:,:,:,termQ)=omega(:,:,:,1)

         Write(*,*)'Imbalance term'
         call callsolvegen(fa,zero2,omega,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,&
              nres)
         omegas(:,:,:,termA)=omega(:,:,:,1)

         Write(*,*)'Boundary conditions'        
         call callsolvegen(zero2,boundaries,omega,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,&
              nres)
         omegas(:,:,:,termB)=omega(:,:,:,1)
                  
       endif

       if(mode.eq.'Q')then
         Write(*,*)'Vorticity advection'
         call callsolveQG(fv,zero2,omega,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,feta2,nres)
         Write(*,*)'Thermal advection'
         call callsolveQG(ft,zero2,omega,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,feta2,nres)
         Write(*,*)'Boundary conditions'        
         call callsolveQG(zero2,boundaries,omega,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,feta2,nres)
       endif 

     end subroutine calculate_omegas

     end module mod_omega

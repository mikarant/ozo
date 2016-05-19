module mod_omega
!  use mod_subrs
  use mod_const
  use mod_wrf_file
  use mod_common_subrs
  implicit none

contains

  subroutine calculate_omegas( t, u, v, w, z, lev, dx, dy, corpar, q, &
       xfrict, yfrict, utend, vtend, ttend, psfc, omegas)

    real,dimension(:,:,:,:),intent(inout) :: omegas
    real,dimension(:,:,:),intent(in) :: t,u,v,w,xfrict,yfrict
    real,dimension(:,:,:),intent(in) :: utend,vtend,ttend
    real,dimension(:,:),  intent(in) :: psfc
    real,dimension(:),intent(in) :: lev,corpar
    real,                 intent(in) :: dx,dy 
    real,dimension(:,:,:),intent(inout) :: z,q

    real,dimension(:,:,:,:),allocatable :: omega
    real,dimension(:,:,:),allocatable :: sigmaraw,omegaan,zetaraw,dum6
    real,dimension(:,:,:),allocatable :: dum4,dum5,dum2,dum1,lapl
    real,dimension(:,:,:),allocatable :: dudp,dvdp,zetatend
       integer nlev,nlon,nlat,nlevmax,nlonmax,nlatmax,nresmax,max3d,max4d    
!
!      Maximum size of the WRF output grid. 
!      nresmax = maximum number of resolutions used in the
!                multigrid algorithm

       parameter (nlonmax=160,nlatmax=320,nlevmax=40,nresmax=7)
       parameter (max3d=nlonmax*nlatmax*nlevmax)
       parameter (max4d=nresmax*max3d)

       real f1,f2 
       real forcing(max3d),ftest(max3d)
       real zero(max3d),boundaries(max3d)
       real omegaold(max3d)
       real sigma(max3d),sigma0(nlevmax)
       real laplome(max3d)
       real domedp2(max3d)
       real coeff(max3d),coeff1(max3d),coeff2(max3d)       
       real df2dp2(max3d)

       real zeta(max3d),d2zetadp(max3d)
       real feta(max3d)     ! coriolis * abs vorticity
       real fv(max3d),ft(max3d),ff(max3d),fq(max3d),fa(max3d)
       real dum0(max3d)
       real dum3(max3d)
       real resid(max3d),omega1(max3d)
       real mulfact(max3d)

       logical lfconst,lzeromean,lcensor        
       real fconst   ! constant value for coriolis parameter
       real alfa     ! relaxation coefficient
        
       integer i,j,k,ijk
!
!      itermax = maximum number of multigrid cycles
!      ny1 = number of iterations for each grid when going to coarser grids
!      ny2 = number of iterations for each grid when returning to finer grids
!
       integer itermax,ny1,ny2

!
!      Number of input fields (may need to be changed in the netCDF version) 
!
       integer, parameter :: nfieldin = 11
!
!      Some natural constants
!
       real, parameter :: a = 6371e3,pii=3.1415926536,r = 287.,cp = 1004.,&
             omg=7.292e-5,g=9.80665
!
!      Threshold values to keep the generalized omega equation elliptic.
!
       real,parameter :: sigmamin=2e-7,etamin=2e-6

!      For iubound, ilbound and iybound are 0, horizontal boundary
!      conditions are used at the upper, lower and north/south boundaries  
!      A value of 1 for any of these parameters means that the boundary
!      condition is taken directly from the "real" WRF omega. In practice,
!      only the lower boundary condition (ilbound) is important.
       integer iubound,ilbound,iybound

!      Horizontal and vertical resolution
       real dlev       
!
!      Mode of operation:
!
       character mode 
!
!      Tables for administrating the different resolutions 
!      (***CHANGED SUBSTANTIALLY FOR THE MULTIGRID VERSION ***)
!
       integer nres
       integer nlonx(nresmax),nlatx(nresmax),nlevx(nresmax)
       real dx2(nresmax),dy2(nresmax),dlev2(nresmax)
       real boundaries2(max4d)
       real zero2(max4d) ! (this is really stupid)
       real sigma2(max4d),feta2(max4d),corpar2(nlatmax*nresmax)
       real d2zetadp2(max4d),dudp2(max4d),dvdp2(max4d)
       real sigma02(nlevmax*nresmax)

       real toler                     

!       namelist/PARAM/infile,outfile,psfile,&
!               alfa,dx,dy,iubound,ilbound,iybound,itermax,&
!               f1,f2,nlon,nlat,nlev,lev1,lev2,t1,t2, &
!               fconst,lfconst,mode,toler,lzeromean,lcensor,ny1,ny2

       iubound=1 ! 1 for "real" omega as upper-boundary condition
       ilbound=1 ! 1 for "real" omega as upper-boundary condition
       iybound=1 ! 1 for "real" omega as north/south boundary condtion
       alfa=0.2  ! relaxation coefficient
       itermax=1000 ! maximum number of iterations
       mode='G'   ! Mode = generalized omega equation
       fconst=1e-4 ! default of coriolis parameter, used in lfconst=.true.
       f1=fconst ! coriolis parameter at southern boundary (if lfconst=.false.)
       f2=fconst ! coriolis parameter at northern boundary (if lfconst=.false.)
       toler=5e-5  ! threshold for stopping iterations
       lfconst=.true.   ! Coriolis parameter treated as constant
       lzeromean=.true. ! Area means of omega are set to zero
       lcensor=.true. ! Forcing below surface is set to zero.
!      For example, if the surface pressure is psfc = 970 hPa
!      and the pressure levels have 50 hPa spacing, the values of
!      vorticity advection, temperature advection etc. at 1000 hPa
!      will be multiplied by 0.4. If psfc <= 950 hPa,
!      the multiplication factor at 1000 hPa will be zero.
       ny1=2 ! number of iterations at each grid resolution when proceeding to coarser
       ny2=2 ! number of iterations at each grid resolution when returning to finer

       nlon=size(t,1)
       nlat=size(t,2)
       nlev=size(t,3)


       zero=0.

!       read(*,nml=PARAM)

       if(mode.eq.'G')write(*,*)'Generalized omega equation'   
       if(mode.eq.'Q')write(*,*)'Quasi-geostrophic omega equation'   
       if(mode.eq.'T')write(*,*)'Generalized test version'   
       if(mode.eq.'t')write(*,*)'Quasigeostrophic test version'   
       if(mode.ne.'G'.and.mode.ne.'Q'.and.mode.ne.'T'.and.mode.ne.'t')then
         write(*,*)'Unknown mode of operation. Aborting'
         stop
       endif

       allocate(sigmaraw(nlon,nlat,nlev),zetaraw(nlon,nlat,nlev))
       allocate(dum1(nlon,nlat,nlev))
       allocate(dum2(nlon,nlat,nlev))
!       allocate(dum3(nlon,nlat,nlev))
       allocate(dum4(nlon,nlat,nlev))
       allocate(dum5(nlon,nlat,nlev))
       allocate(lapl(nlon,nlat,nlev))
       allocate(dum6(nlon,nlat,nlev))
       allocate(dvdp(nlon,nlat,nlev))
       allocate(dudp(nlon,nlat,nlev))
       allocate(zetatend(nlon,nlat,nlev))
       omegaan=w
!
!      Number of different resolutions in solving the equation = nres 
!      Choose so that the coarsest grid has at least 5 points
!
       nres=1+int(log(max(nlon,nlat,nlev)/5.)/log(2.))
       write(*,*)'nres=',nres
!
       allocate(omega(nlon,nlat,nlev,nres))
!      Coriolis parameter as a function of latitude
!
!       do j=1,nlat
!         if(lfconst)then
!          corpar(j)=fconst
!         else
!          corpar(j)=f1+(j-1.)/(nlat-1.)*(f2-f1)
!         endif
!       enddo 
!
!      Pressure levels are converted from hPa to Pa  
!
!       do k=1,nlev
!         lev(k)=100*(lev1+(lev2-lev1)*(k-1.)/(nlev-1.))
!       enddo
!
!      dlev will be negative, if lev2 > lev1. This is OK.
!
!       dlev=100*(lev2-lev1)/(nlev-1.)
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
!      In the QG equation, coriolis parameter is treated
!      as constant, regardless of the namelist value of 'lfconst'

       if(mode.eq.'Q'.or.mode.eq.'t')then
         lfconst=.true.
       endif
!
!      Open output and input files.
!
!       call OPENFW(outfile,1,nlon*nlat*nlev)
!       call OPENFR(infile,11,nlon*nlat*nlev)
!       if(lcensor)call OPENFR(psfile,12,nlon*nlat)
!
!      Counter for writing output fields   
!
!       irec=0
!
!************ The main time loop begins here ************************
!
!       do time=t1,t2
!          write(*,*)'Time=',time
!
!      Read the input files (to be completely revised with netCDF!)
!
!       irecin=(time-1)*nfieldin

!       call reagra2(t,nlon*nlat*nlev,irecin+1,11)         
!       if(mode.eq.'G'.or.mode.eq.'T') call reagra2(u,nlon*nlat*nlev,irecin+2,11) 
!       if(mode.eq.'G'.or.mode.eq.'T') call reagra2(v,nlon*nlat*nlev,irecin+3,11)
!       call reagra2(omegaan,nlon*nlat*nlev,irecin+4,11)  
!       if(mode.eq.'Q')call reagra2(z,nlon*nlat*nlev,irecin+5,11)                  
!       if(mode.eq.'G')then
!         call reagra2(q,nlon*nlat*nlev,irecin+6,11)
!         call reagra2(xfrict,nlon*nlat*nlev,irecin+7,11)
!         call reagra2(yfrict,nlon*nlat*nlev,irecin+8,11)
!         call reagra2(ttend,nlon*nlat*nlev,irecin+9,11)
!         call reagra2(utend,nlon*nlat*nlev,irecin+10,11)
!         call reagra2(vtend,nlon*nlat*nlev,irecin+11,11)
!       endif
!
!      For quasi-geostrophic equation: calculation of geostrophic winds
!
       if(mode.eq.'Q')then
         call gwinds(z,u,v,nlon,nlat,nlev,dx,dy,corpar,g,dum1,dum2)
       endif        
!
!      Multiplication factor for forcing: 
!      1 above the ground, smaller (or 0) below the ground
!
       mulfact=1.
       if(lcensor)then         
!         call reagra2(psfc,nlon*nlat,time,12)
         call calmul(psfc,lev,nlon,nlat,nlev,mulfact) 
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

       call fvort(u,v,zetaraw,nlon,nlat,nlev,corpar,& 
                       dx,dy,dlev,dum1,dum2,dum3,dum4,fv,mulfact)     
       call ftemp(u,v,t,nlon,nlat,nlev,lev,r,& 
                       dx,dy,dum1,dum2,dum3,ft,mulfact)
       endif

       if(mode.eq.'G')then

       call ffrict(xfrict,yfrict,dum1,dum2,nlon,nlat,nlev,corpar,&
                   dx,dy,dlev,ff,mulfact)
       call fdiab(q,nlon,nlat,nlev,lev,r,dx,dy,fq,mulfact)

       call fimbal(zetatend,ttend,nlon,nlat,nlev,corpar,& 
                   r,lev,dx,dy,dlev,dum1,dum2,fa,mulfact)

       endif 
!
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
       call modify(sigmaraw,sigma,sigmamin,etamin,nlon,nlat,nlev,zetaraw,&
                  zeta,feta,corpar,dudp,dvdp)      
!
!      4. Second pressure derivative of vorticity 
!
       call p2der(zeta,d2zetadp,nlon,nlat,nlev,dlev) 
!
!      5. Area mean of static stability over the whole grid
!
       do k=1,nlev
       call aave(sigmaraw(:,:,k),nlon,nlat,sigma0(k))                      
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
           ijk=i+(j-1)*nlon+(k-1)*nlon*nlat
           sigma(ijk)=sigma0(k)
           feta(ijk)=fconst**2.
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

       do k=1,nlev
       do j=1,nlat
       do i=1,nlon
          dum2(i,j,k)=sigma0(k)*omegaan(i,j,k)
       enddo
       enddo
       enddo

       call p2der(omegaan,df2dp2,nlon,nlat,nlev,dlev)
       call laplace_cart(omegaan,lapl,dx,dy)
!
       do k=1,nlev
       do j=1,nlat       
       do i=1,nlon
         ijk=i+(j-1)*nlon+(k-1)*nlon*nlat
         ftest(ijk)=sigma(ijk)*lapl(i,j,k)+feta(ijk)*df2dp2(ijk) 
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
          ijk=i+(j-1)*nlon+(k-1)*nlon*nlat
          forcing(ijk)=omegaan(i,j,k)           
          dum2(i,j,k)=sigmaraw(i,j,k)*omegaan(i,j,k)
       enddo
       enddo
       enddo

       call laplace_cart(dum2,dum1,dx,dy)
       call p2der(omegaan,dum3,nlon,nlat,nlev,dlev)
       call p2der(zetaraw,dum5,nlon,nlat,nlev,dlev)

       do k=1,nlev
       do j=1,nlat       
       do i=1,nlon
         ijk=i+(j-1)*nlon+(k-1)*nlon*nlat
         dum2(i,j,k)=(corpar(j)+zetaraw(i,j,k))*corpar(j)*dum3(ijk)
       enddo
       enddo
       enddo

       do k=1,nlev
       do j=1,nlat       
       do i=1,nlon
         ijk=i+(j-1)*nlon+(k-1)*nlon*nlat
         dum3(ijk)=-corpar(j)*omegaan(i,j,k)*dum5(i,j,k)
       enddo
       enddo
       enddo

       call xder_cart(omegaan,dx,dum4) 
       call yder_cart(omegaan,dy,dum5) 
 
       do k=1,nlev
       do j=1,nlat
       do i=1,nlon
          ijk=i+(j-1)*nlon+(k-1)*nlon*nlat
          dum6(i,j,k)=-corpar(j)*(dvdp(i,j,k)*dum4(i,j,k)-dudp(i,j,k)*dum5(i,j,k))
       enddo
       enddo
       enddo        
       call pder(dum6,dlev,dum4) 

       do k=1,nlev
       do j=1,nlat       
       do i=1,nlon
         ijk=i+(j-1)*nlon+(k-1)*nlon*nlat
         ftest(ijk)=dum1(i,j,k)+dum2(i,j,k)+dum3(ijk)+dum4(i,j,k) 
       enddo
       enddo
       enddo

       endif ! (forcing for the general test case if mode.eq.'T')

!      Write the WRF omega field to the beginning of the output file   

!        irec=irec+1
!        call WRIGRA2(omegaan,nlon*nlat*nlev,irec,1)       
!
!      Boundary conditions from WRF omega?  
!
       boundaries=0.
!
       do j=1,nlat       
       do i=1,nlon
         ijk=i+(j-1)*nlon+(1-1)*nlon*nlat
         boundaries(ijk)=iubound*omegaan(i,j,1)
         ijk=i+(j-1)*nlon+(nlev-1)*nlon*nlat
         boundaries(ijk)=ilbound*omegaan(i,j,nlev)
       enddo
       enddo      
       do k=2,nlev-1
       do i=1,nlon
         ijk=i+(1-1)*nlon+(k-1)*nlon*nlat        
         boundaries(ijk)=iybound*omegaan(i,1,k)
         ijk=i+(nlat-1)*nlon+(k-1)*nlon*nlat        
         boundaries(ijk)=iybound*omegaan(i,nlat,k)
       enddo
       enddo     

!
!      Regrid left-hand-side parameters and boundary conditions to 
!      coarser grids. Note that non-zero boundary conditions are only 
!      possibly given at the highest resolutions (As only the 
!      'residual omega' is solved at lower resolutions)

       do i=1,nres
         ijk=1+nlon*nlat*nlev*(i-1)
         if(i.eq.1)then
           call coarsen3d(boundaries,boundaries2(ijk),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         else
           call coarsen3d(zero,boundaries2(ijk),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         endif 
         call coarsen3d(zero,zero2(ijk),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(sigma,sigma2(ijk),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(feta,feta2(ijk),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(d2zetadp,d2zetadp2(ijk),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(dudp,dudp2(ijk),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(dvdp,dvdp2(ijk),nlon,nlat,nlev,nlonx(i),nlatx(i),nlevx(i))
         call coarsen3d(sigma0,sigma02(1+(i-1)*nlev),1,1,nlev,1,1,nlevx(i))
         call coarsen3d(corpar,corpar2(1+(i-1)*nlat),1,nlat,1,1,nlatx(i),1)
       enddo   
!
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
         call callsolvegen(ftest,boundaries2,omega,omegaold,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,laplome,domedp2,&
              coeff1,coeff2,coeff,dum0,dum1,dum2,dum3,dum4,dum5,dum6,resid,omega1,&
              itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)
       endif 

       if(mode.eq.'t')then
         call callsolveQG(ftest,boundaries2,omega,omegaold,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,feta2,laplome,domedp2,dum1,&
              coeff1,coeff2,coeff,resid,omega1,itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)
       endif 

       if(mode.eq.'G')then            
         Write(*,*)'Vorticity advection'
         call callsolvegen(fv,zero2,omega,omegaold,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,laplome,domedp2,&
              coeff1,coeff2,coeff,dum0,dum1,dum2,dum3,dum4,dum5,dum6,resid,omega1,&
              itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)
         omegas(:,:,:,termV)=omega(:,:,:,1)

         Write(*,*)'Thermal advection'
         call callsolvegen(ft,zero2,omega,omegaold,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,laplome,domedp2,&
              coeff1,coeff2,coeff,dum0,dum1,dum2,dum3,dum4,dum5,dum6,resid,omega1,&
              itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)
         omegas(:,:,:,termT)=omega(:,:,:,1)

         Write(*,*)'Friction'
         call callsolvegen(ff,zero2,omega,omegaold,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,laplome,domedp2,&
              coeff1,coeff2,coeff,dum0,dum1,dum2,dum3,dum4,dum5,dum6,resid,omega1,&
              itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)
         omegas(:,:,:,termF)=omega(:,:,:,1)

         Write(*,*)'Diabatic heating'
         call callsolvegen(fq,zero2,omega,omegaold,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,laplome,domedp2,&
              coeff1,coeff2,coeff,dum0,dum1,dum2,dum3,dum4,dum5,dum6,resid,omega1,&
              itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)
         omegas(:,:,:,termQ)=omega(:,:,:,1)

         Write(*,*)'Imbalance term'
         call callsolvegen(fa,zero2,omega,omegaold,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,laplome,domedp2,&
              coeff1,coeff2,coeff,dum0,dum1,dum2,dum3,dum4,dum5,dum6,resid,omega1,&
              itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)
         omegas(:,:,:,termA)=omega(:,:,:,1)

         Write(*,*)'Boundary conditions'        
         call callsolvegen(zero2,boundaries,omega,omegaold,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,sigma2,feta2,corpar2,d2zetadp2,dudp2,dvdp2,laplome,domedp2,&
              coeff1,coeff2,coeff,dum0,dum1,dum2,dum3,dum4,dum5,dum6,resid,omega1,&
              itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)
         omegas(:,:,:,termB)=omega(:,:,:,1)
                  
       endif

       if(mode.eq.'Q')then
         Write(*,*)'Vorticity advection'
         call callsolveQG(fv,zero2,omega,omegaold,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,feta2,laplome,domedp2,dum1,&
              coeff1,coeff2,coeff,resid,omega1,itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)

         Write(*,*)'Thermal advection'
         call callsolveQG(ft,zero2,omega,omegaold,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,feta2,laplome,domedp2,dum1,&
              coeff1,coeff2,coeff,resid,omega1,itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)
         Write(*,*)'Boundary conditions'        
         call callsolveQG(zero2,boundaries,omega,omegaold,nlonx,nlatx,nlevx,dx2,dy2,dlev2,&
              sigma02,feta2,laplome,domedp2,dum1,&
              coeff1,coeff2,coeff,resid,omega1,itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)
       endif 
!       enddo 

!**************** The main time loop ends here ********************

       close(1)
       close(11)
       if(lcensor)close(12)

     end subroutine calculate_omegas
!
!****************** SUBROUTINES **************************************************
!

       subroutine calmul(psfc,lev,nlon,nlat,nlev,mulfact) 
!
!      Calculation of multiplication factors to attenuate below-surface forcing.
!      (mulfact = 1 if above surface)
!
      implicit none
       integer i,j,k,nlon,nlat,nlev 
       real psfc(nlon,nlat),lev(nlev),mulfact(nlon,nlat,nlev)

       mulfact=1.
       do i=1,nlon
       do j=1,nlat
       do k=2,nlev
         if(psfc(i,j).le.lev(k-1))then
             mulfact(i,j,k)=0.
         else
           if(psfc(i,j).le.lev(k))then
             mulfact(i,j,k)=(psfc(i,j)-lev(k-1))/(lev(k)-lev(k-1))             ! 
           endif
         endif   
       enddo
       enddo
       enddo 
       return
       end subroutine calmul


       subroutine coarsen3D(f,g,nlon1,nlat1,nlev1,nlon2,nlat2,nlev2)
!
!      Averages the values of field f (grid size nlon1 x nlat1 x nlev1) over larger
!      grid boxes (grid size nlon2 x nlat2 x nlev2), to field g

!      To keep the algorithm
!      simple, only 0/1 weights are used -> works only well if
!      div(nlon1/nlon2)=div(nlat1/nlat2)=div(nlev1/nlev2)
!
       implicit none
       integer nlon1,nlat1,nlon2,nlat2,nlev1,nlev2
       real f(nlon1,nlat1,nlev1),g(nlon2,nlat2,nlev2)
       integer i,i2,j,j2,k,k2,imin,imax,jmin,jmax,kmin,kmax
       real fsum
 
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

       return
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
       integer nlon1,nlat1,nlev1,nlon2,nlat2,nlev2
       real f(nlon2,nlat2,nlev2),g(nlon1,nlat1,nlev1)
       integer i,i2,j,j2,k,k2,imin,imax,jmin,jmax,kmin,kmax

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

       return
       end subroutine finen3D

       subroutine gwinds(z,u,v,nlon,nlat,nlev,dx,dy,corpar,g,dzdx,dzdy)
!
!      Calculation of geostrophic winds (u,v) from z. At the equator, mean of
!      the two neighbouring latitudes is used (should not be a relevant case).
!
       implicit none
       integer i,j,k,nlon,nlat,nlev
       real z(nlon,nlat,nlev),u(nlon,nlat,nlev),v(nlon,nlat,nlev)
       real dzdx(nlon,nlat,nlev),dzdy(nlon,nlat,nlev)
       real corpar(nlat),g,dx,dy

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

       return
       end subroutine gwinds


!      subroutine define_sigma(t,dum1,dum2,sigma,nlon,nlat,nlev,lev,r,cp,dlev)
!
!     Calculation of static stability (sigma) from t and the constants r and cp
!
!      implicit none
!      integer i,j,k,nlon,nlat,nlev
!      real t(nlon,nlat,nlev),dum1(nlon,nlat,nlev),dum2(nlon,nlat,nlev),&
!           sigma(nlon,nlat,nlev),lev(nlev)
!      real r,cp,dlev

!       do k=1,nlev
!       do j=1,nlat
!       do i=1,nlon
!         dum1(i,j,k)=log(t(i,j,k))-(r/cp)*log(lev(k)/1e5)
!       enddo
!       enddo
!       enddo
!       call pder(dum1,dum2,nlon,nlat,nlev,dlev) 
!       do k=1,nlev
!       do j=1,nlat
!       do i=1,nlon
!         sigma(i,j,k)=-R*t(i,j,k)/lev(k)*dum2(i,j,k)
!       enddo
!       enddo
!       enddo

!      return
!      end subroutine define_sigma

      subroutine modify(sigmaraw,sigma,sigmamin,etamin,nlon,nlat,nlev,zetaraw,&
                        zeta,feta,corpar,dudp,dvdp)      
!
!      Modifying stability and vorticity to keep the LHS of the genearlized
!      omega equation elliptic
!
      implicit none
      integer i,j,k,nlon,nlat,nlev
      real sigmaraw(nlon,nlat,nlev),sigma(nlon,nlat,nlev)
      real etamin,sigmamin,corpar(nlat)
      real zeta(nlon,nlat,nlev),zetaraw(nlon,nlat,nlev),feta(nlon,nlat,nlev)
      real dudp(nlon,nlat,nlev),dvdp(nlon,nlat,nlev)
     
       do k=1,nlev
       do j=1,nlat
       do i=1,nlon
         sigma(i,j,k)=max(sigmaraw(i,j,k),sigmamin)
         zeta(i,j,k)=0.
!
!         Northern Hemisphere
!
          if(corpar(j).gt.1e-7)then
            zeta(i,j,k)=max(zetaraw(i,j,k),etamin+corpar(j)/(4*sigma(i,j,k))*(dudp(i,j,k)**2.+dvdp(i,j,k)**2.)-corpar(j))
          endif
!
!         Southern Hemisphere
!
          if(corpar(j).lt.-1e-7)then
            zeta(i,j,k)=min(zetaraw(i,j,k),-etamin+corpar(j)/(4*sigma(i,j,k))*(dudp(i,j,k)**2.+dvdp(i,j,k)**2.)-corpar(j))
          endif
          feta(i,j,k)=(zeta(i,j,k)+corpar(j))*corpar(j)
       enddo
       enddo
       enddo

      return
      end subroutine modify

      subroutine aave(f,nlon,nlat,res)                      
!
!     Calculation of area mean (res) of field f in cartesian coordinates. Simplest possible way.
!
      implicit none
      integer i,j,nlon,nlat
      real f(nlon,nlat),res,sum,wsum
!      do k=1,nlev
      sum=0
      wsum=0
      do j=1,nlat
      do i=1,nlon
        sum=sum+f(i,j)
        wsum=wsum+1.
      enddo 
      enddo
      res=sum/wsum
!      enddo

      return
      end subroutine aave

!      subroutine advect_cart(u,v,f,dfdx,dfdy,adv,nlon,nlat,nlev,dx,dy)
!
!     Computing adv = u*dfdx + v*dfdy in cartesian coordinates
!
!      implicit none
!      integer i,j,k,nlon,nlat,nlev
!      real u(nlon,nlat,nlev),v(nlon,nlat,nlev),f(nlon,nlat,nlev)       
!      real dfdx(nlon,nlat,nlev),dfdy(nlon,nlat,nlev),adv(nlon,nlat,nlev)
!      real dx,dy
!      call xder_cart(f,dfdx,nlon,nlat,nlev,dx)          
!      call yder_cart(f,dfdy,nlon,nlat,nlev,dy) 
      
!      do k=1,nlev
!      do j=1,nlat
!      do i=1,nlon
!        adv(i,j,k)=u(i,j,k)*dfdx(i,j,k)+v(i,j,k)*dfdy(i,j,k)
!      enddo
!      enddo
!      enddo

!      return
!      end subroutine advect_cart

      subroutine fvort(u,v,zeta,nlon,nlat,nlev,f,& 
                       dx,dy,dp,eta,detadx,detady,dadvdp,adv,&
                       mulfact)
!
!     Calculation of vorticity advection forcing
!     Input: u,v,zeta
!     Output: stored in "adv" (bad style ...)
!
      implicit none
      integer i,j,k,nlon,nlat,nlev
      real dx,dy,dp
      real u(nlon,nlat,nlev),v(nlon,nlat,nlev),f(nlat)
      real zeta(nlon,nlat,nlev),eta(nlon,nlat,nlev),adv(nlon,nlat,nlev)
      real detadx(nlon,nlat,nlev),detady(nlon,nlat,nlev)
      real dadvdp(nlon,nlat,nlev)
      real mulfact(nlon,nlat,nlev)

      do k=1,nlev
      do j=1,nlat
      do i=1,nlon
        eta(i,j,k)=zeta(i,j,k)+f(j)
      enddo
      enddo
      enddo

      call advect_cart(u,v,eta,dx,dy,adv)
      adv=adv*mulfact
      call pder(adv,dp,dadvdp) 

      do k=1,nlev
      do j=1,nlat      
      do i=1,nlon
         adv(i,j,k)=f(j)*dadvdp(i,j,k)
      enddo
      enddo
      enddo

      return
      end subroutine fvort

      subroutine ftemp(u,v,t,nlon,nlat,nlev,lev,r,& 
                       dx,dy,dtdx,dtdy,lapladv,adv,mulfact)
!
!     Calculation of temperature advection forcing
!     Input: u,v,t
!     Output: stored in "adv" (bad style ...)
!
      implicit none
      integer i,j,k,nlon,nlat,nlev
      real dx,dy,r,lev(nlev)
      real u(nlon,nlat,nlev),v(nlon,nlat,nlev)
      real t(nlon,nlat,nlev),adv(nlon,nlat,nlev),lapladv(nlon,nlat,nlev)
      real dtdx(nlon,nlat,nlev),dtdy(nlon,nlat,nlev)
      real mulfact(nlon,nlat,nlev)

      call advect_cart(u,v,t,dx,dy,adv)
      adv=adv*mulfact
      call laplace_cart(adv,lapladv,dx,dy)         

      do k=1,nlev
      do j=1,nlat      
      do i=1,nlon
          adv(i,j,k)=lapladv(i,j,k)*r/lev(k)
      enddo
      enddo
      enddo

      return
      end subroutine ftemp

      subroutine ffrict(fx,fy,fcurl,dcurldp,nlon,nlat,nlev,f,& 
                       dx,dy,dp,res,mulfact)
!
!     Calculation of friction forcing
!     Input: fx,fy = x and y components of "friction force"
!     Output: res 
!
      implicit none
      integer i,j,k,nlon,nlat,nlev
      real dx,dy,dp
      real fx(nlon,nlat,nlev),fy(nlon,nlat,nlev),f(nlat)
      real fcurl(nlon,nlat,nlev),dcurldp(nlon,nlat,nlev),res(nlon,nlat,nlev)
      real mulfact(nlon,nlat,nlev)

      call curl_cart(fx,fy,dx,dy,fcurl)
      fcurl=fcurl*mulfact
      call pder(fcurl,dp,dcurldp) 

      do k=1,nlev
      do j=1,nlat      
      do i=1,nlon
          res(i,j,k)=-f(j)*dcurldp(i,j,k)
      enddo
      enddo
      enddo

      return
      end subroutine ffrict

      subroutine fdiab(q,nlon,nlat,nlev,lev,r,dx,dy,res,mulfact)
!
!     Calculation of diabatic heaging forcing
!     Input: q = diabatic temperature tendency (already normalized by cp)
!     Output: slightly illogically stored in "adv"
!
      implicit none
      integer i,j,k,nlon,nlat,nlev
      real dx,dy,r,lev(nlev)
      real q(nlon,nlat,nlev),res(nlon,nlat,nlev)
      real mulfact(nlon,nlat,nlev)

      q=q*mulfact
      call laplace_cart(q,res,dx,dy)         

      do k=1,nlev
      do j=1,nlat      
      do i=1,nlon
          res(i,j,k)=-r*res(i,j,k)/lev(k)
      enddo
      enddo
      enddo

      return
      end subroutine fdiab

      subroutine fimbal(dzetadt,dtdt,nlon,nlat,nlev,f,& 
                       r,lev,dx,dy,dp,ddpdzetadt,lapldtdt,res,mulfact)
!
!     Calculation of the FA ("imbalance") forcing term
!     Input: dzetadt, dtdt = vorticity & temperature tendencies
!     Output: res 
!
      implicit none
      integer i,j,k,nlon,nlat,nlev
      real dx,dy,dp,r,lev(nlev)
      real dzetadt(nlon,nlat,nlev),dtdt(nlon,nlat,nlev),f(nlat)
      real ddpdzetadt(nlon,nlat,nlev),lapldtdt(nlon,nlat,nlev),res(nlon,nlat,nlev)
      real mulfact(nlon,nlat,nlev)

      dzetadt=dzetadt*mulfact
      dtdt=dtdt*mulfact

      call pder(dzetadt,dp,ddpdzetadt) 
      call laplace_cart(dtdt,lapldtdt,dx,dy)         

      do k=1,nlev
      do j=1,nlat      
      do i=1,nlon
        ddpdzetadt(i,j,k)=f(j)*ddpdzetadt(i,j,k)
        res(i,j,k)=ddpdzetadt(i,j,k)+lapldtdt(i,j,k)*r/lev(k)
      enddo
      enddo
      enddo

      return
      end subroutine fimbal

       subroutine callsolveQG(rhs,boundaries,omega,omegaold,nlon,nlat,nlev,&
              dx,dy,dlev,sigma0,feta,laplome,domedp2,dum1,&
              coeff1,coeff2,coeff,resid,omega1,itermax,ny1,ny2,alfa, &
              nres,toler,lzeromean)
!
!      Calling solveQG + writing out omega. Multigrid algorithm.
!
       implicit none
       integer i,j,k,nres,nlon(nres),nlat(nres),nlev(nres)
       real dx(nres),dy(nres),dlev(nres)       
       real rhs(nlon(1),nlat(1),nlev(1),nres),boundaries(nlon(1),nlat(1),nlev(1),nres)
       real omega(nlon(1),nlat(1),nlev(1),nres),omegaold(nlon(1),nlat(1),nlev(1),nres) 
       real sigma0(nlev(1),nres)
       real feta(nlon(1),nlat(1),nlev(1),nres)
       real laplome(nlon(1),nlat(1),nlev(1)),domedp2(nlon(1),nlat(1),nlev(1)),dum1(nlon(1),nlat(1),nlev(1))              
       real coeff1(nlon(1),nlat(1),nlev(1)),coeff2(nlon(1),nlat(1),nlev(1)),coeff(nlon(1),nlat(1),nlev(1))
       real resid(nlon(1),nlat(1),nlev(1)),omega1(nlon(1),nlat(1),nlev(1))      
       real alfa,maxdiff,toler
       integer itermax,ny1,ny2,ires
       logical lzeromean
       real aomega
       integer iter

       omega=0.
       omegaold=boundaries
!
!      The whole multigrid cycle is written explicitly here. Better as a separate subroutine?  
!
!----------------------------------------------------------------------------------------
       do iter=1,itermax  ! Each iteration = one (fine->coarse->fine) multigrid cycle
!
!      Loop from finer to coarser resolutions
!
       do ires=1,nres
!         write(*,*)'fine-coarse:iter,ires',iter,ires
         call solveQG(rhs(1,1,1,ires),boundaries(1,1,1,ires),& 
             omega(1,1,1,ires),omegaold(1,1,1,ires),nlon(ires),nlat(ires),nlev(ires),&
             dx(ires),dy(ires),dlev(ires),sigma0(1,ires),feta(1,1,1,ires),&
             laplome,domedp2,coeff1,coeff2,coeff,&
             ny1,alfa,.true.,resid)
         if(ires.eq.1)omega1(:,:,:)=omega(:,:,:,1)
         if(ires.lt.nres)then
            call coarsen3d(resid,rhs(1,1,1,ires+1),nlon(ires),nlat(ires),nlev(ires),&
                 nlon(ires+1),nlat(ires+1),nlev(ires+1))          
         endif  
       enddo         
!
!      Loop from coarser to finer resolutions
!
       do ires=nres-1,1,-1
!        write(*,*)'coarse-fine:iter,ires',iter,ires
         call finen3D(omega(1,1,1,ires+1),dum1,nlon(ires),nlat(ires),nlev(ires),& 
              nlon(ires+1),nlat(ires+1),nlev(ires+1))          
!        Without the underrelaxation (coefficient alfa), the solution diverges
         omegaold(:,:,:,ires)=omega(:,:,:,ires)+alfa*dum1(:,:,:)
!         if(ires.eq.1)then
!         write(*,*)'omega',ires,omega(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2,ires)
!         write(*,*)'dum1',ires,dum1(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2)
!         write(*,*)'omegaold',ires,omegaold(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2,ires)
!         endif

         call solveQG(rhs(1,1,1,ires),boundaries(1,1,1,ires),& 
             omega(1,1,1,ires),omegaold(1,1,1,ires),nlon(ires),nlat(ires),nlev(ires),&
             dx(ires),dy(ires),dlev(ires),sigma0(1,ires),feta(1,1,1,ires),&
             laplome,domedp2,coeff1,coeff2,coeff,&
             ny2,alfa,.false.,resid)
       enddo  

        maxdiff=0.
        do k=1,nlev(1)  
        do j=1,nlat(1)  
        do i=1,nlon(1)  
           maxdiff=max(maxdiff,abs(omega(i,j,k,1)-omega1(i,j,k)))
        enddo
        enddo
        enddo
        if(maxdiff.lt.toler.or.iter.eq.itermax)then
          write(*,*)'iter,maxdiff',iter,maxdiff
          goto 10
        endif

        omegaold=omega
          
        enddo ! iter=1,itermax
 10     continue         
!----------------------------------------------------------------------------------------
!
!       Subtract the area mean of omega
!
         if(lzeromean)then
           do k=1,nlev(1) 
           call aave(omega(1,1,k,1),nlon(1),nlat(1),aomega)
           do j=1,nlat(1)
           do i=1,nlon(1)
             omega(i,j,k,1)=omega(i,j,k,1)-aomega
           enddo 
           enddo 
           enddo 
         endif

!        irec=irec+1
!        call WRIGRA2(omega,nlon(1)*nlat(1)*nlev(1),irec,iunit)

       return
       end subroutine callsolveQG

       subroutine solveQG(rhs,boundaries,omega,omegaold,nlon,nlat,nlev,&
              dx,dy,dlev,sigma0,feta,&
              laplome,domedp2,coeff1,coeff2,coeff,&
              niter,alfa,lres,resid)
!
!      Solving the QG omega equation using 'niter' iterations.
!
       implicit none

       integer i,j,k,nlon,nlat,nlev
       real omegaold(nlon,nlat,nlev),omega(nlon,nlat,nlev)
       real boundaries(nlon,nlat,nlev) 
       real sigma0(nlev),feta(nlon,nlat,nlev),rhs(nlon,nlat,nlev) 
       real laplome(nlon,nlat,nlev),domedp2(nlon,nlat,nlev)        
       real coeff1(nlon,nlat,nlev),coeff2(nlon,nlat,nlev),coeff(nlon,nlat,nlev)
       real dx,dy,dlev,maxdiff,alfa
       integer niter
       logical lres
       real resid(nlon,nlat,nlev)       

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
       call updateQG(omegaold,omega,sigma0,feta,&
            rhs,nlon,nlat,nlev, &
            dx,dy,dlev,maxdiff,laplome,domedp2,&
            coeff1,coeff2,coeff,alfa)
       enddo

       if(lres)then
             call residQG(rhs,omega,resid,sigma0,feta,&
             nlon,nlat,nlev,dx,dy,dlev,laplome,domedp2)
       endif

       return
       end subroutine solveQG

       subroutine updateQG(omegaold,omega,sigma,etasq,rhs,nlon,nlat,nlev, &
            dx,dy,dlev,maxdiff,lapl2,domedp2,coeff1,coeff2,coeff,alfa)
!
!      New estimate for the local value of omega, using omega in the 
!      surrounding points and the right-hand-side forcing (rhs)
!
!      QG version: for 'sigma' and 'etasq', constant values from the QG theory
!      are used.
! 
       implicit none

       integer i,j,k,nlon,nlat,nlev
       real omegaold(nlon,nlat,nlev),omega(nlon,nlat,nlev) 
       real sigma(nlev),etasq(nlon,nlat,nlev),rhs(nlon,nlat,nlev) 
       real lapl2(nlon,nlat,nlev),domedp2(nlon,nlat,nlev)        
       real coeff1(nlon,nlat,nlev),coeff2(nlon,nlat,nlev),coeff(nlon,nlat,nlev)
       real dx,dy,dlev,maxdiff,alfa
!
!      Top and bottom levels: omega directly from the boundary conditions,
!      does not need to be solved.
!
       call laplace2_cart(omegaold,lapl2,coeff1,nlon,nlat,nlev,dx,dy)
       call p2der2(omegaold,domedp2,coeff2,nlon,nlat,nlev,dlev) 

!       write(*,*)'Calculate the coefficients'
!       write(*,*)nlon,nlat,nlev
!       write(*,*)'coeff(nlon,nlat,nlev-1)',coeff(nlon,nlat,nlev-1)
!       write(*,*)'coeff1(nlon,nlat,nlev-1)',coeff1(nlon,nlat,nlev-1)
!       write(*,*)'coeff2(nlon,nlat,nlev-1)',coeff2(nlon,nlat,nlev-1)
!       write(*,*)'sigma(nlon,nlat,nlev-1)',sigma(nlon,nlat,nlev-1)
!       write(*,*)'domedp2(nlon,nlat,nlev-1)',domedp2(nlon,nlat,nlev-1)
!       write(*,*)'lapl2(nlon,nlat,nlev-1)',lapl2(nlon,nlat,nlev-1)
!       write(*,*)'etasq(nlon,nlat,nlev-1)',etasq(nlon,nlat,nlev-1)

       do k=2,nlev-1
       do j=2,nlat-1
       do i=1,nlon
         coeff(i,j,k)=sigma(k)*coeff1(i,j,k)+etasq(i,j,k)*coeff2(i,j,k)
         omega(i,j,k)=(rhs(i,j,k)-sigma(k)*lapl2(i,j,k)-etasq(i,j,k)*domedp2(i,j,k)) &
         /coeff(i,j,k) 
       enddo
       enddo
       enddo

!       write(*,*)'Updating omega'
       maxdiff=0.
       do k=2,nlev-1
       do j=2,nlat-1
       do i=1,nlon
         maxdiff=max(maxdiff,abs(omega(i,j,k)-omegaold(i,j,k)))
         omegaold(i,j,k)=alfa*omega(i,j,k)+(1-alfa)*omegaold(i,j,k)
       enddo
       enddo
       enddo

       return
       end subroutine updateQG


       subroutine residQG(rhs,omega,resid,sigma,etasq,nlon,nlat,nlev,&
            dx,dy,dlev,laplome,domedp2)
!
!      Calculating the residual RHS - LQG(omega)
!      
!      Variables:
!
!      omega = approximation for omega
!      sigma = local values of sigma (*after modifying for ellipticity*)
!      feta = f*eta (*after modifying for ellipticity*)
!      f = coriolis parameter
!      d2zetadp = second pressure derivative of relative vorticity 
!      dudp,dvdp = pressure derivatives of wind components
!      rhs = right-hand-side forcing
!
       implicit none

       integer i,j,k,nlon,nlat,nlev
       real rhs(nlon,nlat,nlev),omega(nlon,nlat,nlev),resid(nlon,nlat,nlev) 
       real sigma(nlev),etasq(nlon,nlat,nlev)
       real laplome(nlon,nlat,nlev),domedp2(nlon,nlat,nlev)        
       real dx,dy,dlev

       call laplace_cart(omega,laplome,dx,dy)         
       call p2der(omega,domedp2,nlon,nlat,nlev,dlev) 
 
       do k=1,nlev
       do j=1,nlat
       do i=1,nlon
         resid(i,j,k)=rhs(i,j,k)-(sigma(k)*laplome(i,j,k)+etasq(i,j,k)*domedp2(i,j,k))
!         if(i.eq.nlon/2.and.j.eq.nlat/2.and.k.eq.nlev/2)write(*,*)'rhs,resid',rhs(i,j,k),resid(i,j,k)
       enddo
       enddo
       enddo

       return
       end subroutine residQG


       subroutine callsolvegen(rhs,boundaries,omega,omegaold,nlon,nlat,nlev,&
              dx,dy,dlev,sigma0,sigma,feta,corpar,d2zetadp,dudp,dvdp,&
              laplome,domedp2,coeff1,coeff2,coeff,dum0,dum1,dum2,dum3,&
              dum4,dum5,dum6,resid,omega1,itermax,ny1,ny2,alfa,&
              nres,toler,lzeromean)
!
!      Calling solvegen + writing out omega. Multigrid algorithm
!            
       implicit none
       integer i,j,k,nres,nlon(nres),nlat(nres),nlev(nres)
       real dx(nres),dy(nres),dlev(nres)       
       real rhs(nlon(1),nlat(1),nlev(1),nres),boundaries(nlon(1),nlat(1),nlev(1),nres)
       real omega(nlon(1),nlat(1),nlev(1),nres),omegaold(nlon(1),nlat(1),nlev(1),nres) 
       real sigma0(nlev(1),nres),sigma(nlon(1),nlat(1),nlev(1),nres)
       real feta(nlon(1),nlat(1),nlev(1),nres),corpar(nlat(1),nres)
       real d2zetadp(nlon(1),nlat(1),nlev(1),nres),dudp(nlon(1),nlat(1),nlev(1),nres),dvdp(nlon(1),nlat(1),nlev(1),nres)
       real laplome(nlon(1),nlat(1),nlev(1)),domedp2(nlon(1),nlat(1),nlev(1))        
       real coeff1(nlon(1),nlat(1),nlev(1)),coeff2(nlon(1),nlat(1),nlev(1)),coeff(nlon(1),nlat(1),nlev(1))
       real dum0(nlon(1),nlat(1),nlev(1)),dum1(nlon(1),nlat(1),nlev(1)),dum2(nlon(1),nlat(1),nlev(1))
       real dum3(nlon(1),nlat(1),nlev(1)),dum4(nlon(1),nlat(1),nlev(1)),dum5(nlon(1),nlat(1),nlev(1))
       real dum6(nlon(1),nlat(1),nlev(1))      
       real resid(nlon(1),nlat(1),nlev(1)),omega1(nlon(1),nlat(1),nlev(1))      
       real alfa,maxdiff,toler
       integer itermax,ny1,ny2,ires,iter
       logical lzeromean
       real aomega
 
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
         call solvegen(rhs(1,1,1,ires),boundaries(1,1,1,ires),& 
             omega(1,1,1,ires),omegaold(1,1,1,ires),nlon(ires),nlat(ires),nlev(ires),&
             dx(ires),dy(ires),dlev(ires),sigma0(1,ires),sigma(1,1,1,ires),&
             feta(1,1,1,ires),corpar(1,ires),&
             d2zetadp(1,1,1,ires),dudp(1,1,1,ires),dvdp(1,1,1,ires),&
             laplome,domedp2,coeff1,coeff2,coeff,dum0,dum1,dum2,dum3,&
             dum4,dum5,dum6,ny1,alfa,.true.,resid)
!             write(*,*)'ires,omega',ires,omega(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2,1)
!             write(*,*)'ires,resid',ires,resid(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2)
         if(ires.eq.1)omega1(:,:,:)=omega(:,:,:,1)
         if(ires.lt.nres)then
            call coarsen3d(resid,rhs(1,1,1,ires+1),nlon(ires),nlat(ires),nlev(ires),&
                 nlon(ires+1),nlat(ires+1),nlev(ires+1))          
!             write(*,*)'ires,rhs',ires,rhs(nlon(ires+1)/2,nlat(ires+1)/2,nlev(ires+1)/2,ires)
         endif  
       enddo         
!
!      Loop from coarser to finer resolutions
!
       do ires=nres-1,1,-1
 !        write(*,*)'coarse-fine:iter,ires',iter,ires
         call finen3D(omega(1,1,1,ires+1),dum1,nlon(ires),nlat(ires),nlev(ires),& 
              nlon(ires+1),nlat(ires+1),nlev(ires+1))          
!      Without the underrelaxation (coeffient alfa) the soultion diverges
          omegaold(:,:,:,ires)=omega(:,:,:,ires)+alfa*dum1(:,:,:)
!         if(ires.eq.1)then
!         write(*,*)'omega',ires,omega(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2,ires)
!         write(*,*)'dum1',ires,dum1(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2)
!         write(*,*)'omegaold',ires,omegaold(nlon(ires)/2,nlat(ires)/2,nlev(ires)/2,ires)
!         endif

         call solvegen(rhs(1,1,1,ires),boundaries(1,1,1,ires),& 
             omega(1,1,1,ires),omegaold(1,1,1,ires),nlon(ires),nlat(ires),nlev(ires),&
             dx(ires),dy(ires),dlev(ires),sigma0(1,ires),sigma(1,1,1,ires),&
             feta(1,1,1,ires),corpar(1,ires),&
             d2zetadp(1,1,1,ires),dudp(1,1,1,ires),dvdp(1,1,1,ires),&
             laplome,domedp2,coeff1,coeff2,coeff,dum0,dum1,dum2,dum3,&
             dum4,dum5,dum6,ny2,alfa,.false.,resid)
       enddo  

        maxdiff=0.
        do k=1,nlev(1)  
        do j=1,nlat(1)  
        do i=1,nlon(1)  
          maxdiff=max(maxdiff,abs(omega(i,j,k,1)-omega1(i,j,k)))
        enddo
        enddo
        enddo
        if(maxdiff.lt.toler.or.iter.eq.itermax)then
          write(*,*)'iter,maxdiff',iter,maxdiff
          goto 10
        endif

        omegaold=omega
          
        enddo ! iter=1,itermax
 10     continue
!----------------------------------------------------------------------------------

!       Subtracting area mean of omega
         if(lzeromean)then
           do k=1,nlev(1) 
           call aave(omega(1,1,k,1),nlon(1),nlat(1),aomega)
           do j=1,nlat(1)
           do i=1,nlon(1)
             omega(i,j,k,1)=omega(i,j,k,1)-aomega
           enddo 
           enddo 
           enddo 
         endif

!        irec=irec+1
!         call WRIGRA2(omega,nlon(1)*nlat(1)*nlev(1),irec,iunit)

       return
       end subroutine callsolvegen
 

       subroutine solvegen(rhs,boundaries,omega,omegaold,nlon,nlat,nlev,&
              dx,dy,dlev,sigma0,sigma,feta,corpar,d2zetadp,dudp,dvdp,&
              laplome,domedp2,coeff1,coeff2,coeff,dum0,dum1,dum2,dum3,&
              dum4,dum5,dum6,niter,alfa,lres,resid)
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

       integer i,j,k,nlon,nlat,nlev
       real omegaold(nlon,nlat,nlev),omega(nlon,nlat,nlev) 
       real sigma(nlon,nlat,nlev),feta(nlon,nlat,nlev),rhs(nlon,nlat,nlev) 
       real boundaries(nlon,nlat,nlev) 
       real sigma0(nlev),corpar(nlat)
       real d2zetadp(nlon,nlat,nlev),dudp(nlon,nlat,nlev),dvdp(nlon,nlat,nlev)
       real laplome(nlon,nlat,nlev),domedp2(nlon,nlat,nlev)        
       real coeff1(nlon,nlat,nlev),coeff2(nlon,nlat,nlev),coeff(nlon,nlat,nlev)
       real dx,dy,dlev,maxdiff
       real dum0(nlon,nlat,nlev),dum1(nlon,nlat,nlev),dum2(nlon,nlat,nlev)
       real dum3(nlon,nlat,nlev),dum4(nlon,nlat,nlev),dum5(nlon,nlat,nlev)
       real dum6(nlon,nlat,nlev)    
       real resid(nlon,nlat,nlev)  
       real alfa
       integer niter
       logical lres

!       omegaold=boundaries

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
       call updategen(omegaold,omega,sigma0,sigma,feta,corpar,&
            d2zetadp,dudp,dvdp,rhs,nlon,nlat,nlev, &
            dx,dy,dlev,maxdiff,laplome,domedp2,&
            coeff1,coeff2,coeff,dum0,dum1,dum3,dum4,dum5,dum6,alfa)
       enddo
!
!      Calculate the residual = RHS - L(omega)

       if(lres)then
         call residgen(rhs,omega,resid,&
            sigma,feta,corpar,d2zetadp,dudp,dvdp,nlon,nlat,nlev, &
            dx,dy,dlev,dum0,dum1,dum2,dum3,dum4,dum5,dum6)
       endif

       return
       end subroutine solvegen

       subroutine updategen(omegaold,omega, &
            sigma0,sigma,feta,f,d2zetadp,dudp,dvdp,rhs,nlon,nlat,nlev, &
            dx,dy,dlev,maxdiff,lapl2,domedp2,coeff1,coeff2,coeff, &
            dum0,dum1,dum3,dum4,dum5,dum6,alfa)
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
       real omegaold(nlon,nlat,nlev),omega(nlon,nlat,nlev) 
       real sigma(nlon,nlat,nlev),feta(nlon,nlat,nlev),rhs(nlon,nlat,nlev) 
       real sigma0(nlev),f(nlat)
       real d2zetadp(nlon,nlat,nlev),dudp(nlon,nlat,nlev),dvdp(nlon,nlat,nlev)
       real lapl2(nlon,nlat,nlev),domedp2(nlon,nlat,nlev)        
       real coeff1(nlon,nlat,nlev),coeff2(nlon,nlat,nlev),coeff(nlon,nlat,nlev)
       real dx,dy,dlev,maxdiff
       real dum0(nlon,nlat,nlev),dum1(nlon,nlat,nlev)
       real dum3(nlon,nlat,nlev),dum4(nlon,nlat,nlev),dum5(nlon,nlat,nlev)
       real dum6(nlon,nlat,nlev)      
       real alfa
!
!      Top and bottom levels: omega directly from the boundary conditions,
!      does not need to be solved.
!
       call laplace2_cart(omegaold,lapl2,coeff1,nlon,nlat,nlev,dx,dy)
       call p2der2(omegaold,domedp2,coeff2,nlon,nlat,nlev,dlev) 
!
!      Calculate non-constant terms on the left-hand-side, based on 'omegaold'
!
!       a) Deviation of sigma from its normal value

       do k=2,nlev-1
       do j=1,nlat
       do i=1,nlon
        dum0(i,j,k)=omegaold(i,j,k)*(sigma(i,j,k)-sigma0(k))
       enddo
       enddo
       enddo        
       call laplace_cart(dum0,dum1,dx,dy)         
!
!      b) f*omega*(d2zetadp): explicitly, later
!          
!      c) tilting
!
       call xder_cart(omegaold,dx,dum4) 
       call yder_cart(omegaold,dy,dum5) 
 
       do k=1,nlev
       do j=1,nlat
       do i=1,nlon
         dum6(i,j,k)=f(j)*(dudp(i,j,k)*dum5(i,j,k)-dvdp(i,j,k)*dum4(i,j,k))
       enddo
       enddo
       enddo        
       call pder(dum6,dlev,dum3) 
!
!      Solving for omega 
!      Old values are retained at y and z boundaries.
!       
       do k=2,nlev-1
       do j=2,nlat-1
       do i=1,nlon
         coeff(i,j,k)=sigma0(k)*coeff1(i,j,k)+feta(i,j,k)*coeff2(i,j,k)-f(j)*d2zetadp(i,j,k)
         omega(i,j,k)=(rhs(i,j,k)-dum1(i,j,k)-dum3(i,j,k)-sigma0(k)*lapl2(i,j,k)-feta(i,j,k)*domedp2(i,j,k)) &
         /coeff(i,j,k)
       enddo
       enddo
       enddo

!       write(*,*)'updating omega'
       maxdiff=0.
       do k=2,nlev-1
       do j=2,nlat-1
       do i=1,nlon
         maxdiff=max(maxdiff,abs(omega(i,j,k)-omegaold(i,j,k)))
         omegaold(i,j,k)=alfa*omega(i,j,k)+(1-alfa)*omegaold(i,j,k)
       enddo
       enddo
       enddo

       return
       end subroutine updategen

       subroutine residgen(rhs,omega,resid, &
            sigma,feta,f,d2zetadp,dudp,dvdp,nlon,nlat,nlev, &
            dx,dy,dlev,dum0,dum1,dum2,dum3,dum4,dum5,dum6)
!
!      Calculating the residual RHS - L(omega)
!      
!      Variables:
!
!      omega = approximation for omega
!      sigma = local values of sigma (*after modifying for ellipticity*)
!      feta = f*eta (*after modifying for ellipticity*)
!      f = coriolis parameter
!      d2zetadp = second pressure derivative of relative vorticity 
!      dudp,dvdp = pressure derivatives of wind components
!      rhs = right-hand-side forcing
!
       implicit none

       integer i,j,k,nlon,nlat,nlev
       real rhs(nlon,nlat,nlev),omega(nlon,nlat,nlev),resid(nlon,nlat,nlev) 
       real f(nlat)
       real sigma(nlon,nlat,nlev),feta(nlon,nlat,nlev)
       real d2zetadp(nlon,nlat,nlev),dudp(nlon,nlat,nlev),dvdp(nlon,nlat,nlev)
       real dx,dy,dlev
       real dum0(nlon,nlat,nlev),dum1(nlon,nlat,nlev),dum2(nlon,nlat,nlev)
       real dum3(nlon,nlat,nlev),dum4(nlon,nlat,nlev),dum5(nlon,nlat,nlev)
       real dum6(nlon,nlat,nlev)      
!
!      Calculate L(omega)

!       a) nabla^2(sigma*omega)

       dum0=omega*sigma

       call laplace_cart(dum0,dum1,dx,dy)         
!
!      f*eta*d2omegadp
!       
       call p2der(omega,dum2,nlon,nlat,nlev,dlev) 
!
       dum3=feta*dum2
!
!      c) -f*omega*(d2zetadp): explicitly, later
!                  
!      d) tilting
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
         resid(i,j,k)=rhs(i,j,k)-(dum1(i,j,k)+dum2(i,j,k)+dum3(i,j,k)-f(j)*d2zetadp(i,j,k)*omega(i,j,k))
       enddo
       enddo
       enddo

       return
       end subroutine residgen


!       subroutine curl_cart(u,v,zeta,nlon,nlat,nlev,dx,dy)
!
!      Relative vorticity in cartesian coordinates.
!      
!      The domain is assumed to be periodic in east-west-direction
!      ** At the northern and southern boundaries, one-sided y derivatives are used.**
!
!       integer i,i1,i2,j,k,nlon,nlat,nlev
!       real u(nlon,nlat,nlev),v(nlon,nlat,nlev),zeta(nlon,nlat,nlev)
!       real dudy,dvdx
!       real dx,dy

!       do k=1,nlev
!       do i=1,nlon 
!         i1=i-1
!         i2=i+1
!         if(i1.lt.1)i1=i1+nlon
!         if(i2.gt.nlon)i2=i2-nlon
!         do j=1,nlat
!            if(j.gt.1.and.j.lt.nlat)then
!              dudy=(u(i,j+1,k)-u(i,j-1,k))/(2*dy)
!            else
!              if(j.eq.1)dudy=(u(i,j+1,k)-u(i,j,k))/dy
!              if(j.eq.nlat)dudy=(u(i,j,k)-u(i,j-1,k))/dy
!            endif
!            dvdx=(v(i2,j,k)-v(i1,j,k))/(2*dx)
!            zeta(i,j,k)=dvdx-dudy
!         enddo
!       enddo
!       enddo

!       return
!       end subroutine curl_cart

!       subroutine laplace_cart(f,lapl,nlon,nlat,nlev,dx,dy)         
!
!      Laplace operator in cartesian coordinates
!
!      The domain is assumed to be periodic in east-west-direction
!      ** At the northern and southern boundaries, second y derivative is assumed to be zero      
!      
!       integer i,i1,i2,j,k,nlon,nlat,nlev
!       real f(nlon,nlat,nlev),lapl(nlon,nlat,nlev)
!       double precision d2fdy,d2fdx
!       real dx,dy

!       do k=1,nlev
!       do i=1,nlon 
!         i1=i-1
!         i2=i+1
!         if(i1.lt.1)i1=i1+nlon
!         if(i2.gt.nlon)i2=i2-nlon
!         do j=1,nlat
!            d2fdx=(f(i2,j,k)+f(i1,j,k)-2*f(i,j,k))/(dx**2.)
!            if(j.gt.1.and.j.lt.nlat)then
!              d2fdy=(f(i,j+1,k)+f(i,j-1,k)-2*f(i,j,k))/(dy**2.)
!              lapl(i,j,k)=d2fdx+d2fdy
!            else
!              lapl(i,j,k)=d2fdx
!            endif
!         enddo
!       enddo
!       enddo

!       return
!       end subroutine laplace_cart

       subroutine laplace2_cart(f,lapl2,coeff,nlon,nlat,nlev,dx,dy)
!
!      As laplace_cart but
!        - the contribution of the local value to the Laplacian is left out
!        - coeff is the coefficient for the local value
!
       integer i,i1,i2,j,k,nlon,nlat,nlev
       real f(nlon,nlat,nlev),lapl2(nlon,nlat,nlev),coeff(nlon,nlat,nlev)
       double precision d2fdy,d2fdx
       real dx,dy

       do k=1,nlev
       do i=1,nlon 
         i1=i-1
         i2=i+1
         if(i1.lt.1)i1=i1+nlon
         if(i2.gt.nlon)i2=i2-nlon
         do j=2,nlat-1
            d2fdx=(f(i2,j,k)+f(i1,j,k))/(dx**2.)
            if(j.gt.1.and.j.lt.nlat)then
               d2fdy=(f(i,j+1,k)+f(i,j-1,k))/(dy**2.)
               lapl2(i,j,k)=d2fdx+d2fdy
               coeff(i,j,k)=-2/(dy**2.)-2/(dx**2.)
            else
               lapl2(i,j,k)=d2fdx
               coeff(i,j,k)=-2/(dx**2.)
            endif
         enddo
       enddo 
       enddo
 
       return
       end subroutine laplace2_cart

!       subroutine pder(f,dfdp,nlon,nlat,nlev,dp) 
!
!      Estimation of pressure derivatives.
!      One-sided derivatives are used at top and bottom levels
!
!       implicit none
!       integer i,j,k,nlon,nlat,nlev
!       real f(nlon,nlat,nlev),dfdp(nlon,nlat,nlev),dp 
!       do i=1,nlon
!       do j=1,nlat
!       do k=2,nlev-1
!         dfdp(i,j,k)=(f(i,j,k+1)-f(i,j,k-1))/(2.*dp)
!       enddo
!       dfdp(i,j,1)=(f(i,j,2)-f(i,j,1))/dp
!       dfdp(i,j,nlev)=(f(i,j,nlev)-f(i,j,nlev-1))/dp
 
!       enddo
!       enddo
!       return
!       end subroutine pder      


       subroutine p2der(f,df2dp2,nlon,nlat,nlev,dp) 
!
!      Estimation of second pressure derivatives.
!      At top and bottom levels, these are set to zero
!
       implicit none
       integer i,j,k,nlon,nlat,nlev
       real f(nlon,nlat,nlev),df2dp2(nlon,nlat,nlev),dp 
       do i=1,nlon
       do j=1,nlat
       do k=2,nlev-1
         df2dp2(i,j,k)=(f(i,j,k+1)+f(i,j,k-1)-2*f(i,j,k))/(dp*dp)
       enddo
       df2dp2(i,j,1)=0.
       df2dp2(i,j,nlev)=0.
 
       enddo
       enddo
       return
       end subroutine p2der      

       subroutine p2der2(f,df2dp22,coeff,nlon,nlat,nlev,dp) 
!
!      As p2der, but
!        - the contribution of the local value is left out
!        - the coefficient 'coeff' of the local value is also calculated        
!
       implicit none
       integer i,j,k,nlon,nlat,nlev
       real f(nlon,nlat,nlev),df2dp22(nlon,nlat,nlev),coeff(nlon,nlat,nlev),dp 
       do i=1,nlon
       do j=1,nlat
       do k=2,nlev-1
         df2dp22(i,j,k)=(f(i,j,k+1)+f(i,j,k-1))/(dp*dp)
         coeff(i,j,k)=-2/(dp*dp)
       enddo
       df2dp22(i,j,1)=0.
       df2dp22(i,j,nlev)=0.
       coeff(i,j,1)=0.
       coeff(i,j,nlev)=0.
 
       enddo
       enddo
       return
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

     end module mod_omega

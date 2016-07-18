module mod_subrs
  implicit none
contains
!
! This module contains all the Zwack-Okossi-specified subroutines.
!
!------------SUBROUTINES--------------
!
  subroutine calculate_tendencies(omegas,t,u,v,w,z,lev,dx,dy,corpar,q,xfrict,&
       yfrict,ztend,ttend,zeta,zetatend,uKhi,vKhi,mulfact,calc_b,hTends)

!   This is the main subroutine of solving the Zwack-Okossi equation. Input 
!   arguments are variables from WRF and omegas, and output of this subroutine
!   is height tendencies, stored in hTends.
    use mod_const
    use mod_wrf_file
    use mod_common_subrs
    implicit none

    real,dimension(:,:,:,:),intent(in) :: omegas
    real,dimension(:,:,:),  intent(in) :: t,u,v,w,xfrict,yfrict,zeta
    real,dimension(:,:,:),  intent(in) :: ttend,mulfact,uKhi,vKhi
    real,dimension(:),      intent(in) :: lev,corpar
    real,                   intent(in) :: dx,dy 
    logical,                intent(in) :: calc_b
    real,dimension(:,:,:,:),intent(inout) :: hTends
    real,dimension(:,:,:),  intent(inout) :: z,q,ztend,zetatend

    real,dimension(:,:,:,:),allocatable :: vortTends
    real,dimension(:,:,:),allocatable :: sigma,sp,tadv,tadvs
    real,dimension(:,:,:),allocatable :: vorTend_omegaWRF
    real,dimension(:,:,:,:),allocatable :: temptend
    integer :: nlon,nlat,nlev
    real :: dlev

    nlon=size(t,1); nlat=size(t,2); nlev=size(t,3)
    
    allocate(vortTends(nlon,nlat,nlev,n_terms))
    allocate(temptend(nlon,nlat,nlev,n_terms))
    allocate(vorTend_omegaWRF(nlon,nlat,nlev))

    vorTend_omegaWRF=0.
    vortTends=0.
   
    dlev=lev(2)-lev(1)

!   Interpolation of 1000mb geopotential height and height tendency
!   p_interp does not do that.
    call interp1000(z,ztend,t,ttend)

    call vorticity_tendencies(omegas,u,v,w,uKhi,vKhi,zeta,zetatend,dx,dy,&
                        corpar,dlev,xfrict,yfrict,ztend,vortTends,mulfact,&
                        vorTend_omegaWRF)

!   Calculation of stability sigma and Sp (stability parameter)
    sigma = define_sigma(t,lev,dlev)
    sp = define_sp(sigma,lev)

!   Calculation of thermal advection
    tadv = advect_cart(u,v,t,dx,dy)
    tadv=-tadv*mulfact ! Minus sign because it's usually considered as negative

!   Calculation of thermal advection by nondivergent/irrotational wind
    tadvs = advect_cart(uKhi,vKhi,t,dx,dy)
    tadvs=-tadvs*mulfact

!   Vorticity tendencies were multiplied by 'mulfact' in  'vorticity_tendencies'
    q=q*mulfact

!   Calculate height tendencies
    call zwack_okossi(vortTends,vorTend_omegaWRF,w,tadv,tadvs,q,omegas,sp,&
                      corpar,dx,dy,calc_b,hTends,ztend)

!   Area mean correction
    call ht_correction(hTends,temptend,lev,omegas,sp,tadv,tadvs,q,calc_b)

  end subroutine calculate_tendencies

  subroutine ht_correction(hTends,temptend,lev,omegas,sp,tadv,tadvs,q,calc_b)

!   This subroutine does area mean correction for height tendencies
!   by using mean temperature tendencies and hypsometric equation.
!    use mod_type_vars
    use mod_wrf_file
    use mod_const
    implicit none
    real,dimension(:,:,:,:),intent(in) :: omegas
    real,dimension(:,:,:),intent(in) :: sp,tadv,tadvs,q
    real,dimension(:),    intent(in) :: lev
    real,dimension(:,:,:,:),intent(inout) :: hTends
    real,dimension(:,:,:,:),intent(inout) :: temptend
    logical,intent(in) :: calc_b

    integer :: i,j,k,l,m,nlon,nlat,nlev
    real,dimension(:,:),    allocatable :: mtt,mht,mzo,diff
    real,dimension(:),      allocatable :: pres,help1,help2,diffsum
    integer,dimension(:),   allocatable :: terms

    nlon=size(q,1);nlat=size(q,2);nlev=size(q,3)

    allocate(pres(nlev),terms(n_terms-1),diffsum(nlev))
    allocate(help1(n_terms),help2(n_terms))
    allocate(mtt(n_terms,nlev),mht(n_terms,nlev))
    allocate(mzo(n_terms,nlev),diff(n_terms,nlev))
    
    terms=(/1,2,3,4,5,6,7/)
    if (calc_b) terms=(/1,2,3,4,5,6,7,8/)


!   Pressure levels       
    pres=lev(1:size(lev))/100.

!   Initializate variables    
    help1=0.
    help2=0.
    temptend=0.

!   Temperature tendencies of all terms
    do m=1,size(terms)
       l=terms(m)
       temptend(:,:,:,l)=sp*omegas(:,:,:,l)
    enddo
    temptend(:,:,:,termT)=temptend(:,:,:,termT)+tadv
    temptend(:,:,:,termQ)=temptend(:,:,:,termQ)+q
    temptend(:,:,:,termTkhi)=temptend(:,:,:,termTkhi)+tadvs

!   Area mean temperature and calculated area mean height tendencies
!   for each pressure level separately
    do m=1,size(terms)
       l=terms(m)
       do k=1,nlev
          do i=1,nlon
             do j=1,nlat
                help1(l)=help1(l)+temptend(i,j,k,l)
                help2(l)=help2(l)+hTends(i,j,k,l)
             enddo
          enddo
          mtt(l,k)=help1(l)/(nlon*nlat)
          help1=0.
          mzo(l,k)=help2(l)/(nlon*nlat)
          help2=0.   
       enddo
    enddo

!   Area mean temperature and calculated height tendency
!   for layers (1000-950,1000-900,1000-850...)
    do m=1,size(terms)
       l=terms(m)
       do k=2,nlev
          mtt(l,k-1)=sum(mtt(l,1:k))/k
          mzo(l,k-1)=sum(mzo(l,1:k))/k
       enddo
    enddo

!   Area mean "thickness" tendencies, calculated from temperature tendencies
!   by using hypsometric equation
    do m=1,size(terms)
       l=terms(m)
       do k=2,nlev
          mht(l,k-1)=(r/g)*mtt(l,k-1)*log(1000./pres(k))
       enddo
    enddo

!   Difference between hypsometric thickness tendencies and zwack-okossi-
!   calculated thickness tendencies
    do m=1,size(terms)
       l=terms(m)
       diff(l,:)=mht(l,:)-mzo(l,:)
    enddo

!   Adding the difference to calculated height tendencies
    diffsum=0.
    do m=1,size(terms)
       l=terms(m)
       do k=2,nlev
          hTends(:,:,k,l)=hTends(:,:,k,l)+diff(l,k-1)
          if(l<=5)then
             diffsum(k)=diffsum(k)+diff(l,k-1)
          endif
       enddo
    enddo
!    do k=1,nlev
!       hTends(:,:,k,termB)=hTends(:,:,k,termB)-diffsum(k)
!    enddo
  end subroutine ht_correction
      
  subroutine vorticity_tendencies(omegas,u,v,w,uKhi,vKhi,zeta,zetatend,&
                                  dx,dy,corpar,dlev,xfrict,yfrict,ztend,&
                                  vortTends,mulfact,vorTend_omegaWRF)
!   This function calculates both direct and indirect (associated with vertical
!   motions) vorticity tendencies of all five forcings (vorticity and thermal
!   advection, friction, diabatic heating and ageostrophic vorticity tendency). 
!   Input: omegas,u,v and vorticity
!   output: vorticity tendencies of six forcings 

    use mod_const
    use mod_wrf_file
    use mod_common_subrs
    implicit none

    real,dimension(:,:,:,:),intent(in) :: omegas
    real,dimension(:,:,:),intent(in) :: u,v,w,zeta,zetatend,uKhi,vKhi,ztend
    real,dimension(:,:,:),intent(in) :: xfrict,yfrict,mulfact
    real,dimension(:),intent(in) :: corpar
    real,dimension(:,:,:),intent(inout) :: vorTend_omegaWRF
    real,dimension(:,:,:,:),intent(inout) :: vortTends
    
    real,dimension(:,:,:,:,:),allocatable :: vTend
    real,dimension(:,:,:,:),allocatable :: vTend_omegaWRF
    real,dimension(:,:,:),allocatable :: eta,vadv,vadvs,avortt,fvort
    real :: dlev,dx,dy
    integer :: i,j,nlon,nlat,nlev

    nlon=size(u,1); nlat=size(u,2); nlev=size(u,3)

    allocate(eta(nlon,nlat,nlev))
    allocate(avortt(nlon,nlat,nlev))
    allocate(fvort(nlon,nlat,nlev))
    allocate(vTend(nlon,nlat,nlev,3,n_terms))
    allocate(vTend_omegaWRF(nlon,nlat,nlev,3))

    vTend_omegaWRF=0.
    vTend=0.

    do j=1,nlat
       eta(:,j,:)=zeta(:,j,:)+corpar(j)
    enddo

!   Omega-related vorticity equation terms
    do i=1,n_terms
       call vorterms(omegas(:,:,:,i),dx,dy,eta,u,v,zeta,dlev,vTend(:,:,:,:,i))
    enddo

    call vorterms(w,dx,dy,eta,u,v,zeta,dlev,vTend_omegaWRF)

!   Vorticity advection term
    vadv = advect_cart(u,v,eta,dx,dy)
    vadv=-vadv*mulfact ! Minus sign because it's considered negative

!   Irrotational/nondivergent vorticity advection term
    vadvs = advect_cart(uKhi,vKhi,eta,dx,dy)
    vadvs=-vadvs*mulfact

!   Friction-induced vorticity tendency
    fvort = curl_cart(xfrict,yfrict,dx,dy)
    fvort=fvort*mulfact

!   Ageostrophic vorticity tendency
    call ageo_tend(zetatend,ztend,dx,dy,corpar,avortt)
    avortt=avortt*mulfact

!   Totalling all terms, both direct and indirect effects:
    do i=1,n_terms
       do j=1,3
          vortTends(:,:,:,i)=vortTends(:,:,:,i)+&
                                 vTend(:,:,:,j,i) 
       enddo
    enddo
    vortTends(:,:,:,termV)=vadv+vortTends(:,:,:,termV)
    vortTends(:,:,:,termf)=fvort+vortTends(:,:,:,termF)
    vortTends(:,:,:,termA)=avortt+vortTends(:,:,:,termA)
    vortTends(:,:,:,termVKhi)=vadvs+vortTends(:,:,:,termVKhi)
    
    ! Vorticity tendency with WRF omega
    do j=1,3
       vorTend_omegaWRF=vorTend_omegaWRF+vTend_omegaWRF(:,:,:,j)
    enddo
    vorTend_omegaWRF=vorTend_omegaWRF+vadv+fvort+avortt

  end subroutine vorticity_tendencies

  subroutine vorterms(omega,dx,dy,eta,u,v,zeta,dlev,vortt)

!   This function calculates omega-related terms of vorticity equation
    implicit none
    real,dimension(:,:,:),intent(in) :: omega,eta,u,v,zeta
    real,                 intent(in) :: dx,dy,dlev  
    real,dimension(:,:,:,:),intent(inout) :: vortt

!   Vertical advection of vorticity
    call f2(omega,zeta,dlev,vortt(:,:,:,1))
       
!   Divergence term
    call f3(omega,dlev,eta,vortt(:,:,:,2))

!   Tilting/twisting term
    call f4(omega,u,v,dx,dy,dlev,vortt(:,:,:,3))
    
  end subroutine vorterms

  subroutine f2(omega,zeta,dlev,vortadv)
    use mod_common_subrs
!   This function calculates the vertical vorticity advection term,
!   stored in vortadv.

    implicit none
    real,dimension(:,:,:),intent(in) :: omega,zeta
    real,                 intent(in) :: dlev
    real,dimension(:,:,:),intent(inout) :: vortadv
    
    real,dimension(:,:,:),allocatable :: dvortdp

!   Pressure derivative of vorticity
    dvortdp = pder(zeta,dlev)

!   Vertical advection
    vortadv=-omega*dvortdp

  end subroutine f2

  subroutine f3(omega,dlev,eta,vortdiv)
    use mod_common_subrs
    implicit none
!   This function calculates divergence term of vorticity equation,
!   stored in vortdiv.
    real,dimension(:,:,:),intent(in) :: omega,eta
    real,                 intent(in) :: dlev
    real,dimension(:,:,:),intent(inout) :: vortdiv
    
    real,dimension(:,:,:),allocatable :: domegadp

!   Pressure derivative of omega
    domegadp = pder(omega,dlev)

!   Product of absolute vorticity and pressure derivative of omega
    vortdiv=eta*domegadp
       
  end subroutine f3

  subroutine f4(omega,u,v,dx,dy,dlev,vortt4)
    use mod_common_subrs
!   This function calculates tilting/twisting term of vorticity equation,
!   stored in vortt4.      
    implicit none
    real,dimension(:,:,:),intent(in) :: omega,u,v
    real,                 intent(in) :: dx,dy,dlev
    real,dimension(:,:,:),intent(inout) :: vortt4  
    real,dimension(:,:,:),allocatable :: domegadx,domegady,dudp,dvdp
       
!   Gradient of omega
    domegadx = xder_cart(omega,dx)
    domegady = yder_cart(omega,dy)

!   Pressure derivative of wind vectors
    dudp = pder(u,dlev)
    dvdp = pder(v,dlev)

    vortt4=(dudp*domegady-dvdp*domegadx)

  end subroutine f4
  
  subroutine ageo_tend(zetatend,ztend,dx,dy,corpar,avortt)
    use mod_const
    use mod_common_subrs
    implicit none
!   This function calculates ageostrophic vorticity tendency (needed in 
!   ageostrophic vorticity tendency forcing).
!   Input: Real vorticity tendency (zetatend), real height tendency (ztend)
!   Output: Ageostrophic vorticity tendency (avortt)

    real,dimension(:,:,:),intent(in) :: zetatend,ztend
    real,dimension(:),intent(in) :: corpar
    real,                 intent(in) :: dx,dy
    real,dimension(:,:,:),intent(inout) :: avortt

    real,dimension(:,:,:),allocatable :: lapl,gvort
    integer :: nlon,nlat,nlev,j
    nlon=size(ztend,1); nlat=size(ztend,2); nlev=size(ztend,3)
    allocate(gvort(nlon,nlat,nlev),lapl(nlon,nlat,nlev))

!   Laplacian of height tendency
    lapl=0.
    call laplace_cart(ztend,lapl,dx,dy)

!   Geostrophic vorticity tendency
    do j=1,nlat
       gvort(:,j,:)=(g/corpar(j))*lapl(:,j,:)
    enddo

!   Ageostrophic vorticity tendency is equal to real vorticity tendency 
!   minus geostrophic vorticity tendency.
    avortt=-(zetatend-gvort) 

  end subroutine ageo_tend

  subroutine zwack_okossi(vortTends,vorTend_omegaWRF,w,tadv,tadvs,q,omegas,&
                          sp,corpar,dx,dy,calc_b,hTends,ztend)

!   This function calculates zwack-okossi equation for all forcings.
!   Input: vorticity tendencies of forcings, omegas, thermal advection and
!   diabatic heating.
!   Output: Height tendencies of different forcingterms.  
    use mod_const
    use mod_wrf_file
    use mod_common_subrs
    use mod_poisson_DFT
    implicit none
    
    real,dimension(:,:,:,:),intent(in) :: omegas
    real,dimension(:,:,:,:),intent(inout) :: vortTends
    real,dimension(:,:,:),intent(in) :: tadv,tadvs,q,sp
    real,                 intent(in) :: dx,dy
    real,dimension(:,:,:),intent(in) :: w,ztend
    real,dimension(:),intent(in) :: corpar
    logical,intent(in) :: calc_b
    real,dimension(:,:,:),intent(inout) :: vorTend_omegaWRF
    real,dimension(:,:,:,:),intent(inout) :: hTends

    real,dimension(:,:,:),allocatable :: corf,laplz
    real,dimension(:,:,:),allocatable :: ttend_omegaWRF,gvtend_omegaWRF
    real,dimension(:,:,:,:),allocatable :: tempTends,gvortTends
    integer :: nlon,nlat,nlev,i,j,k
    double precision, dimension ( : , : ), allocatable :: &
         bd_ay, bd_by, bd_0
        
    nlon=size(q,1); nlat=size(q,2); nlev=size(q,3)

    allocate(tempTends(nlon,nlat,nlev,n_terms))
    allocate(gvortTends(nlon,nlat,nlev,n_terms))
    allocate ( bd_ay ( nlon + 1, nlev ), bd_by ( nlon + 1, nlev ) )
    allocate ( bd_0 ( nlon + 1, nlev ) )
    allocate(laplz(nlon,nlat,nlev),ttend_omegaWRF(nlon,nlat,nlev))
    allocate(corf(nlon,nlat,nlev))
    allocate(gvtend_omegaWRF(nlon,nlat,nlev))

    ttend_omegaWRF=0.
    gvtend_omegaWRF=0.

    do j=1,nlat
       corf(:,j,:)=corpar(j)
    enddo

!   Temperature tendencies
    do i=1,n_terms
       tempTends(:,:,:,i)=sp*omegas(:,:,:,i)       
    enddo
    tempTends(:,:,:,termT)=tadv+tempTends(:,:,:,termT)
    tempTends(:,:,:,termQ)=q+tempTends(:,:,:,termQ)
    tempTends(:,:,:,termTKhi)=tadvs+tempTends(:,:,:,termTKhi)

!   Calculation of geostrophic vorticity tendencies of all forcings using 
!   zwack-okossi vorticity equation.
    
    do k = 1, nlev
       do i = 2, nlon
          bd_ay ( i, k ) = ( ztend ( i, 1, k ) + &
               ztend ( i - 1, 1, k ) ) / 2
          bd_by ( i, k ) = ( ztend ( i, nlat, k ) + &
               ztend ( i - 1, nlat, k ) ) / 2
       enddo
    enddo
    bd_ay ( 1, : ) = ( ztend ( 1, 1, : ) + &
         ztend ( nlon, 1, : ) ) / 2
    bd_by ( 1, : ) = ( ztend ( 1, nlat, : ) + &
         ztend ( nlon, nlat, : ) ) / 2
    bd_ay ( nlon + 1, : ) = bd_ay ( 1, : )
    bd_by ( nlon + 1, : )  =bd_by ( 1, : )

    bd_0 = 0.0e0
 
! Integration
    do i=1,n_terms
       call zo_integral(vortTends(:,:,:,i),tempTends(:,:,:,i),dx,dy,corf,&
            gvortTends(:,:,:,i))
    enddo

! Height tendency with WRF omega
    ttend_omegaWRF=sp*w+tadv+q
    call zo_integral(vorTend_omegaWRF,ttend_omegaWRF,dx,dy,corf,&
         gvtend_omegaWRF)

    do k=1,nlev
       do i=1,5
          ! Five first terms with zero y-boundaries
          call poisson_solver_2D( gvortTends( :, :, k, i ), &
               dx, dy, hTends(:,:,k,i), bd_0 ( :, k ), bd_0 ( :, k ) )
       enddo
       ! B-term
       if (calc_b) then
          call poisson_solver_2D( gvortTends ( :, :, k, 8 ), &
               dx, dy, hTends(:,:,k,8), bd_0 ( :, k ), bd_0 ( :, k ) )
       end if
       ! "Pseudo" height tendency
       call poisson_solver_2D( gvortTends ( :, :, k, termVKhi ), & 
            dx, dy, hTends(:,:,k,termVKhi), bd_ay ( :, k ), bd_by ( :, k ) )
       ! WRF omega height tendency
       call poisson_solver_2D( gvortTends ( :, :, k, termTKhi ), & 
            dx, dy, hTends(:,:,k,termTKhi), bd_ay ( :, k ), bd_by ( :, k ) )
       ! WRF omega height tendency
       call poisson_solver_2D( gvtend_omegaWRF ( :, :, k ), & 
            dx, dy, hTends(:,:,k,termTKhi), bd_ay ( :, k ), bd_by ( :, k ) )
    enddo
    
  end subroutine zwack_okossi

  subroutine zo_integral(vorttend,temptend,dx,dy,corpar,geo_vort)
!
!   This function calculates integrals of zwack-okossi equation. 
!   It is done by slightly undocumented way.
!   
    use mod_const
    use mod_common_subrs
    implicit none
    real,dimension(:,:,:),intent(in) :: vorttend,temptend,corpar
    real,                 intent(in) :: dx,dy
    real,dimension(:,:,:),intent(inout) :: geo_vort 

    real,dimension(:,:,:),allocatable ::lapltemp,inttemp,int_tot
    real,dimension(:,:),allocatable :: temp_mean,vort_mean
    integer :: nlon,nlat,nlev,k

    nlon=size(vorttend,1)
    nlat=size(vorttend,2)
    nlev=size(vorttend,3)
    allocate(temp_mean(nlon,nlat),vort_mean(nlon,nlat))
    allocate(lapltemp(nlon,nlat,nlev),inttemp(nlon,nlat,nlev))
    allocate(int_tot(nlon,nlat,nlev))

!   Vertical mean of vorticity tendency. It's multiplied by coriolisparameter 
!   so that temperature tendency doesn't have to be divided by f 
!   (problems in equator). 
!   vort_mean=sum(vorttend,dim=3)/nlev
!
!   Changed 9.12.2015: half-weight for lowest and highest level 
!   (hope this is correct ...) 
!
    call vertave(vorttend,vort_mean,nlon,nlat,nlev)
    vort_mean(:,:)=vort_mean(:,:)*corpar(:,:,1)

!   Laplacian of temperature tendency
    lapltemp=0.
    call laplace_cart(temptend,lapltemp,dx,dy)

!   Divide laplacian of temperature tendency by "pressure" (Actually it is
!   divided only by k (index), but that has been taken into account in the
!   integral by not multiplying it by dp) and multiplying it by R.
    do k=1,nlev
       lapltemp(:,:,k)=r*lapltemp(:,:,k)/(nlev-k+1)
    enddo

!   Integration from p to pl. Pl goes from surface to the top of the atmosphere.
    inttemp(:,:,1)=0
    do k=2,nlev
       inttemp(:,:,k)=inttemp(:,:,k-1)+(lapltemp(:,:,k)+lapltemp(:,:,k-1))/2.
    enddo

!   Mean of that integration (outer integral of the equation)      
!    temp_mean=sum(inttemp,dim=3)/nlev
!
!   Changed 9.12.2015: half-weight for lowest and highest level 
!   (hope this is correct ...) 
!
    call vertave(inttemp,temp_mean,nlon,nlat,nlev)

!   Actually not sure what's going on in this line (copied from matlab-scripts)
    do k=1,nlev
       int_tot(:,:,k)=-temp_mean(:,:)+inttemp(:,:,k)
    enddo

!   Sum of vorticity tendency and integrated temperature tendency
    do k=1,nlev
       geo_vort(:,:,k)=vort_mean(:,:)+int_tot(:,:,k)
    enddo

!   Dividing by gravitational acceleration
    geo_vort=geo_vort/g
       
  end subroutine zo_integral

  subroutine vertave(f,fmean,nlon,nlat,nlev)
!
!   Vertical average, half-weight to top and bottom levels
!
    implicit none
    real,dimension(:,:,:),intent(in) :: f
    real,dimension(:,:),intent(inout) :: fmean

    integer :: i,j,k,nlon,nlat,nlev

    do i=1,nlon
       do j=1,nlat
          fmean(i,j)=(f(i,j,1)+f(i,j,nlev))/2.
          do k=2,nlev-1
             fmean(i,j)=fmean(i,j)+f(i,j,k) 
          enddo
          fmean(i,j)=fmean(i,j)/(nlev-1.)
       enddo
    enddo

  end subroutine vertave

  subroutine interp1000(z,ztend,t,ttend)
!
!   Interpolation of 1000 hPa geopotential height from 950 hPa height by using
!   mean temperature of 1000-950 hPa layer.
    use mod_const
    implicit none
    real,dimension(:,:,:),intent(inout) :: z,ztend
    real,dimension(:,:,:),intent(in) :: t,ttend
    real,dimension(:,:),allocatable :: meanTemp,meanTempTend
    integer :: nlon,nlat,i,j

    nlon=size(z,1)
    nlat=size(z,2)
    allocate(meanTemp(nlon,nlat),meanTempTend(nlon,nlat))
          
    do i=1,nlon
       do j=1,nlat
          meanTemp(i,j)=(t(i,j,1)+t(i,j,2))/2.
          meanTempTend(i,j)=(ttend(i,j,1)+ttend(i,j,2))/2.
          z(i,j,1)=z(i,j,2)-(r/g)*meanTemp(i,j)*log(100000./95000.)
          ztend(i,j,1)=ztend(i,j,2)-(r/g)*meanTempTend(i,j)*log(100000./95000.)
       enddo
    enddo
    
  end subroutine interp1000

end module mod_subrs

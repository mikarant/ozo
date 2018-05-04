module mod_wrf_file
  use netcdf
  implicit none

  integer :: n_terms = 7
  integer,parameter :: termV=1,& ! names of omega and height tendency terms
       termT=2,&
       termF=3,&
       termQ=4,&
       termA=5,&
       termVKhi=6,&
       termTKhi=7
  character ( 6 ), dimension ( 7 ), parameter :: &
       omega_term_names &
       = [ 'ome_v ', 'ome_t ', 'ome_f ', 'ome_q ', 'ome_a ', &
       'ome_vs', 'ome_ts' ], &
       htend_term_names &
       = [ 'htv   ', 'htt   ', 'htf   ', 'htq   ', 'hta   ', &
       'htvKhi', 'httKhi' ]
  character ( 8 ), dimension ( 2 ), parameter :: &
       QG_omega_term_names &
       = [ 'ome_v_QG', 'ome_t_QG' ]
  character ( 5 ), parameter :: ome_b_name='ome_b'
  character ( 3 ), parameter :: htend_b_name='htb'
  character ( 9 ), parameter :: ztend_name='ztend_WRF'
  character ( 7 ), parameter :: ome_name='ome_WRF'
  character ( 11 ), dimension ( 4 ), parameter :: rname = &
       [ 'west_east  ', 'south_north', 'vlevs      ', &
       'time       ' ]
  character ( 37 ), dimension ( 2 ), parameter :: QG_omega_long_names = &
       [ 'QG omega due to vorticity advection  ', &
       'QG omega due to temperature advection' ]
  character ( 49 ), dimension ( 7 ), parameter :: omega_long_names = &
       [ 'omega due to vorticity advection                 ', &
       'omega due to temperature advection               ', &
       'omega due to friction                            ', &
       'omega due to diabatic heating                    ', &
       'omega due to ageostrophic vorticity tendency     ', &
       'omega due to vorticity advection by irrot. wind  ',&
       'omega due to temperature advection by irrot. wind' ]
  character ( 59 ), dimension ( 7 ), parameter :: htend_long_names = &
       [ 'height tendency due to vorticity advection                 ', &
       'height tendency due to temperature advection               ', &
       'height tendency due to friction                            ', &
       'height tendency due to diabatic heating                    ', &
       'height tendency due to ageostrophic vorticity tendency     ', &
       'height tendency due to vorticity advection by irrot. wind  ',&
       'height tendency due to temperature advection by irrot. wind' ]
  character ( 49 ), parameter :: ome_b_long_name = &
       'omega due to boundary conditions                 '
  character ( 59 ), parameter :: htend_b_long_name = &
       'height tendency due to boundary conditions                 '

  type wrf_file
     integer :: ncid, dims ( 4 )
     integer :: wrf_cu_phys
     !     real :: dx, dy
     real, dimension ( : ), allocatable :: times, pressure_levels, corpar
     integer, dimension ( : ), allocatable :: xdim,ydim
     integer, dimension ( : ), allocatable :: nlon, nlat, nlev
     real, dimension ( : ), allocatable :: dx, dy, dlev
  end type wrf_file

contains

  function create_out_file ( fname, wrf_infile, mode, calc_b, forc ) result ( f )
    character ( * ),   intent ( in ) :: fname ! file name
    type ( wrf_file ), intent ( in ) :: wrf_infile ! wrf inputfile
    logical,           intent ( in ) :: calc_b, forc
    type ( wrf_file ) :: f
    integer :: dimids(4),i,status,varid,varids(4)
    character :: mode

    ! Create a new netcdf file
    call check( nf90_create ( fname, NF90_CLOBBER, f % ncid ) )

    ! Copy dimension information from wrf input file
    f % dims = wrf_infile % dims
    f % pressure_levels = wrf_infile % pressure_levels

    ! Create dimensions and their variables to the new file
    do i=1,3
       call def_dim(f%ncid,rname(i),dimids(i),f%dims(i),.FALSE.)
    enddo
    call def_dim(f%ncid,rname(4),dimids(4),f%dims(4),.TRUE.)

    ! Add axis attributes to dimensions
    do i = 1, 4
       call check ( nf90_inq_varid ( &
            f % ncid, trim ( rname ( i ) ), varid ) )
       varids ( i ) = varid
    end do
    call check( nf90_put_att(f % ncid, varids (1),&
         trim('standard_name'),trim('longitude') ) )
    call check( nf90_put_att(f % ncid, varids (1),&
         trim('units'),trim('degrees_east') ) )
    call check( nf90_put_att(f % ncid, varids (1),&
         trim('axis'),trim('X') ) )
    call check( nf90_put_att(f % ncid, varids (2),&
         trim('standard_name'),trim('latitude') ) )
    call check( nf90_put_att(f % ncid, varids (2),&
         trim('units'),trim('degrees_north') ) )
    call check( nf90_put_att(f % ncid, varids (2),&
         trim('axis'),trim('Y') ) )
    call check( nf90_put_att(f % ncid, varids (3),&
         trim('standard_name'),trim('air_pressure') ) )
    call check( nf90_put_att(f % ncid, varids (3),&
         trim('units'),trim('hPa') ) )
    call check( nf90_put_att(f % ncid, varids (3),&
         trim('positive'),trim('down') ) )
    call check( nf90_put_att(f % ncid, varids (3),&
         trim('axis'),trim('Z') ) )
    call check( nf90_put_att(f % ncid, varids (4),&
         trim('standard_name'),trim('time') ) )
    call check( nf90_put_att(f % ncid, varids (4),&
         trim('units'),trim('hours since 1997-01-01 00:00:00') ) )
    call check( nf90_put_att(f % ncid, varids (4),&
         trim('calendar'),trim('standard') ) )

    ! Create quasi-geostrophic omega variables
    if (mode.eq.'Q')then
       do i = 1, size ( QG_omega_term_names )
          status = nf90_def_var ( f % ncid, trim ( QG_omega_term_names ( i ) ),&
               NF90_FLOAT, dimids, varid )
          if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
               call check ( status )
          call check( nf90_put_att(f % ncid, varid,trim('description'), &
               trim(QG_omega_long_names(i)) ) )
          call check( nf90_put_att(f % ncid, varid,trim('units'),&
               trim('Pa s-1') ) )
       end do

    else

       ! Create generalized omega variables
       do i = 1, size ( omega_term_names )
          status = nf90_def_var ( f % ncid, trim ( omega_term_names ( i ) ), &
               NF90_FLOAT, dimids, varid )
          if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
               call check ( status )
          call check( nf90_put_att(f % ncid, varid,trim('description'),&
               trim(omega_long_names(i)) ) )
          call check( nf90_put_att(f % ncid, varid,trim('units'),&
               trim('Pa s-1') ) )
       end do

       if (calc_b) then ! create omega b-variable if wanted
          status = nf90_def_var ( f % ncid, trim ( ome_b_name ), &
               NF90_FLOAT, dimids, varid )
          if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
               call check ( status )
          call check( nf90_put_att(f % ncid, varid,trim('description'),&
               trim(ome_b_long_name) ) )
          call check( nf90_put_att(f % ncid, varid,trim('units'),&
               trim('Pa s-1') ) )
       end if

       ! Create ome_WRF variable
       status = nf90_def_var ( f % ncid, trim ( ome_name ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('omega from WRF') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('Pa s-1') ) )

       ! Create height tendency variables
       do i = 1, size ( htend_term_names )
          status = nf90_def_var ( f % ncid, trim ( htend_term_names ( i ) ), &
               NF90_FLOAT, dimids, varid )
          if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
               call check ( status )
          call check( nf90_put_att(f % ncid, varid,trim('description'),&
               trim(htend_long_names(i)) ) )
          call check( nf90_put_att(f % ncid, varid,trim('units'),&
               trim('m s-1') ) )
       end do

       if (calc_b) then ! create htend b-variable if wanted
          status = nf90_def_var ( f % ncid, trim ( htend_b_name ), &
               NF90_FLOAT, dimids, varid )
          if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
               call check ( status )
          call check( nf90_put_att(f % ncid, varid,trim('description'),&
               trim(htend_b_long_name) ) )
          call check( nf90_put_att(f % ncid, varid,trim('units'),&
               trim('m s-1') ) )
       end if

       ! Create ztend_wrf variable
       status = nf90_def_var ( f % ncid, trim ( ztend_name ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('height tendency from WRF') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('m s-1') ) )
    endif

    ! Create z variable
    status = nf90_def_var ( f % ncid, trim ( 'GHT' ), NF90_FLOAT, &
         dimids, varid )
    if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
         call check ( status )
    call check( nf90_put_att(f % ncid, varid,trim('description'),&
         trim('geopotential height') ) )
    call check( nf90_put_att(f % ncid, varid,trim('units'),&
         trim('gpm') ) )

    ! Create psfc variable
    status = nf90_def_var ( f % ncid, trim ( 'PSFC' ), NF90_FLOAT, &
         (/ 1, 2, 4/), varid )
    if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
         call check ( status )
    call check( nf90_put_att(f % ncid, varid,trim('description'),&
         trim('sfc pressure') ) )
    call check( nf90_put_att(f % ncid, varid,trim('units'),&
         trim('Pa') ) )

    if (forc) then

       ! Create vadv variable
       status = nf90_def_var ( f % ncid, trim ( 'vadv' ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('vorticity advection') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('1/s^2') ) )
       ! Create tadv variable
       status = nf90_def_var ( f % ncid, trim ( 'tadv' ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('temperature advection') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('K/s') ) )
       ! Create friction variable
       status = nf90_def_var ( f % ncid, trim ( 'fvort' ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('vorticity tendency due to friction') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('1/s^2') ) )
       ! Create diabatic heating variable
       status = nf90_def_var ( f % ncid, trim ( 'diab' ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('Diabatic heating') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('K/s') ) )
       ! Create ageostrophic vorticity tendency variable
       status = nf90_def_var ( f % ncid, trim ( 'ageo' ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('Ageostrophic vorticity tendency') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('1/s^2') ) )
    end if

    ! Stop defining mode
    call check( nf90_enddef ( f % ncid ) )

    ! Copy dimension data to the new file type-variable
    allocate(f % xdim (f%dims(1)))
    allocate(f % ydim (f%dims(2)))

    f % xdim = (/(i, i=1,f%dims(1), 1)/)
    f % ydim = (/(i, i=1,f%dims(2), 1)/)

    print*,"Outputfile created!"
  end function create_out_file

  function open_wrf_file ( fname ) result ( f )
    character ( * ), intent ( in ) :: fname
    type ( wrf_file ) :: f
    integer :: i, dimid, varid, nres

    print*,"Opening file: ",fname
    call check( nf90_open ( fname, NF90_WRITE, f % ncid ) )

    print*,"Inquiring dimensions from the input file..."
    do i = 1, 4
       call check ( nf90_inq_dimid ( &
            f % ncid, trim ( rname ( i ) ), dimid ) )
       call check ( nf90_inquire_dimension ( &
            f % ncid, dimid, len = f % dims ( i ) ) )
    end do

    ! Number of resolutions in the solving of omega equation
    nres=1+int(log(max(f%dims(1),f%dims(2),f%dims(3))/5.)/log(2.))
    allocate(f % nlon ( nres ), f % nlat ( nres ), f % nlev ( nres ) )
    allocate(f % dx ( nres ), f % dy ( nres ), f % dlev ( nres ) )

    print*,"Getting attributes from the input file..."
    call check ( nf90_get_att ( &
         f % ncid, NF90_GLOBAL, 'DX', f % dx(1) ) )
    call check ( nf90_get_att ( &
         f % ncid, NF90_GLOBAL, 'DY', f % dy(1) ) )
    call check ( nf90_get_att ( &
         f % ncid, NF90_GLOBAL, 'CU_PHYSICS', f % wrf_cu_phys ) )

    print*,"Getting pressure level information from the input file..."
    allocate ( f % pressure_levels ( f % dims ( 3 ) ) )
    call check ( nf90_inq_varid ( f % ncid, 'LEV', varid ) )
    call check ( nf90_get_var ( f % ncid, varid, f % pressure_levels, &
         start = [ 1 ], count = [ size ( f % pressure_levels ) ] ) )

    f % nlon(1) = f % dims(1)
    f % nlat(1) = f % dims(2)
    f % nlev(1) = f % dims(3)
    f % dlev(1) = f % pressure_levels(2) - f % pressure_levels(1)

    ! Number of different resolutions in solving the equation = nres
    ! Choose so that the coarsest grid has at least 5 points
    do i=2,nres
       f % nlon(i) = max(f % nlon(i-1)/2,5)
       f % nlat(i) = max(f % nlat(i-1)/2,5)
       f % nlev(i) = max(f % nlev(i-1)/2,5)
       f % dx(i) = f % dx(1)*real(f % nlon(1))/real(f % nlon(i))
       f % dy(i) = f % dy(1)*real(f % nlat(1))/real(f % nlat(i))
       f % dlev(i) = f % dlev(1)*real(f % nlev(1))/real(f % nlev(i))
    enddo

    print*,"Getting time information from the input file..."
    allocate ( f % times ( f % dims ( 4 ) ) )
    !    f % times = (/ ( i, i = 0, f % dims ( 4 ) - 1, 1 ) /)
    call check ( nf90_inq_varid ( f % ncid, 'XTIME', varid ) )
    call check ( nf90_get_var ( f % ncid, varid, f % times, &
         start = [ 1 ], &
         count = [ size ( f % times ) ] ) )
    f % times = f % times * 60

    print*,"Getting coriolisparameter from the input file..."
    allocate ( f % corpar ( f % dims ( 2 ) ) )
    call check ( nf90_inq_varid ( f % ncid, 'F', varid ) )
    call check ( nf90_get_var ( f % ncid, varid, f % corpar, &
         start = [ 1, 1, 2 ], &
         count = [ 1, size ( f % corpar ), 1 ] ) )

    print*,"Input file opened succesfully!"
  end function open_wrf_file

  subroutine close_wrf_file ( file )
    type ( wrf_file ) :: file
    call check ( nf90_close ( file % ncid ) )
  end subroutine close_wrf_file

  function open_out_file ( fname ) result ( f )
    character ( * ), intent ( in ) :: fname
    type ( wrf_file ) :: f
    integer :: i,dimid

    call check( nf90_open ( fname, NF90_WRITE, f % ncid ) )

    do i = 1, 4
       call check ( nf90_inq_dimid ( &
            f % ncid, trim ( rname ( i ) ), dimid ) )
       !       dimids ( i ) = dimid
       call check ( nf90_inquire_dimension ( &
            f % ncid, dimid, len = f % dims ( i ) ) )
    end do

  end function open_out_file

  function real2d ( file, time, names )
    type ( wrf_file ), intent ( in ) :: file
    integer, intent ( in ) :: time
    character ( * ), dimension ( : ), intent ( in ) :: names
    real, dimension ( :, : ), allocatable :: real2d
    integer :: i

    real2d = data ( names ( 1 ) )
    do i = 2, size ( names )
       real2d = real2d + data ( names ( i ) )
    end do

  contains

    function data ( name )
      real, dimension ( :, : ), allocatable :: data
      character ( * ) :: name
      integer :: varid

      allocate ( data ( file % dims ( 1 ), file % dims ( 2 ) ) )
      call check ( nf90_inq_varid ( file % ncid, trim ( name ), varid ) )
      call check ( nf90_get_var ( file % ncid, varid, data, &
           start = [ 1, 1, time ], count = [ shape ( data ), 1 ] ) )

    end function data

  end function real2d

  function real3d ( file, time, names )
    type ( wrf_file ), intent ( in ) :: file
    integer, intent ( in ) :: time
    character ( * ), dimension ( : ), intent ( in ) :: names
    real, dimension ( :, :, : ), allocatable :: real3d
    integer :: i

    !    print*,'Reading ',trim(names(1))
    real3d = data ( names ( 1 ) )

    do i = 2, size ( names )
       !       print*,'Reading ',trim(names(i))
       real3d = real3d + data ( names ( i ) )
    end do

  contains

    function data ( name )
      real, dimension ( :, :, : ), allocatable :: data
      character ( * ) :: name
      integer :: varid

      allocate ( data ( &
           file % dims ( 1 ), file % dims ( 2 ), file % dims ( 3 ) ) )

      call check ( nf90_inq_varid ( file % ncid, trim ( name ), varid ) )
      call check ( nf90_get_var ( file % ncid, varid, data, &
           start = [ 1, 1, 1, time ], count = [ shape ( data ), 1 ] ) )
    end function data

  end function real3d

  subroutine check ( status )
    integer, intent ( in ) :: status
    if ( status /= NF90_NOERR ) then
       write(*,*) 'Error in ', trim ( nf90_strerror ( status ) )
       stop 1
    end if
  end subroutine check

  subroutine def_dim(ncid,dim_name,dimid,len,unlimit)
    implicit none
    character(*) :: dim_name
    integer :: ncid,dimid,varid,len
    logical :: unlimit

    if(unlimit)then
       call check( nf90_def_dim(ncid, trim(dim_name), NF90_UNLIMITED, dimid) )
    else
       call check( nf90_def_dim(ncid, trim(dim_name), len, dimid) )
    endif
    call check( nf90_def_var(ncid, trim(dim_name), NF90_REAL, dimid, varid) )

  end subroutine def_dim

end module mod_wrf_file

module mod_wrf_file
  use netcdf
  implicit none

  integer,parameter :: n_terms = 8
  integer,parameter :: termV=1,& ! names of omega and height tendency terms
                       termT=2,&
                       termF=3,&
                       termQ=4,&
                       termA=5,&
                       termB=6,&
                       termVKhi=7,&
                       termTKhi=8
  character ( 6 ), dimension ( n_terms ), parameter :: &
       omega_term_names &
       = [ 'ome_v ', 'ome_t ', 'ome_f ', 'ome_q ', 'ome_a ', 'ome_b ', &
       'ome_vs', 'ome_ts' ], &
       htend_term_names &
       = [ 'htv   ', 'htt   ', 'htf   ', 'htq   ', 'hta   ', 'htb   ', &
       'htvKhi', 'httKhi' ]

  type wrf_file
     integer :: ncid, dims ( 4 )
     integer :: wrf_cu_phys
     real :: dx, dy
     real, dimension ( : ), allocatable :: times, pressure_levels, corpar
  end type wrf_file

contains

  function open_wrf_file ( fname ) result ( f )
    character ( * ), intent ( in ) :: fname
    type ( wrf_file ) :: f
    character ( 4 ), dimension ( 4 ), parameter :: rname = &
         [ 'lon ', 'lat ', 'lev ', 'Time' ]
    integer :: i, dimid, varid, dimids ( 4 ), status

    call check( nf90_open ( fname, NF90_WRITE, f % ncid ) )
    
    do i = 1, 4
       call check ( nf90_inq_dimid ( &
            f % ncid, trim ( rname ( i ) ), dimid ) )
       dimids ( i ) = dimid
       call check ( nf90_inquire_dimension ( &
            f % ncid, dimid, len = f % dims ( i ) ) )
    end do
    
    call check ( nf90_get_att ( &
         f % ncid, NF90_GLOBAL, 'DX', f % dx ) )
    call check ( nf90_get_att ( &
         f % ncid, NF90_GLOBAL, 'DY', f % dy ) )
    call check ( nf90_get_att ( &
         f % ncid, NF90_GLOBAL, 'CU_PHYSICS', f % wrf_cu_phys ) )
    
    allocate ( f % times ( f % dims ( 4 ) ) )
    call check ( nf90_inq_varid ( f % ncid, 'Time', varid ) )
    call check ( nf90_get_var ( f % ncid, varid, f % times, &
         start = [ 1 ], &
         count = [ size ( f % times ) ] ) )
    f % times = f % times * 3600

    allocate ( f % pressure_levels ( f % dims ( 3 ) ) )
    call check ( nf90_inq_varid ( f % ncid, 'PRES', varid ) )
    call check ( nf90_get_var ( f % ncid, varid, f % pressure_levels, &
         start = [ 1, 1, 1, 2 ], &
         count = [ 1, 1, size ( f % pressure_levels ), 1 ] ) )

    allocate ( f % corpar ( f % dims ( 2 ) ) )
    call check ( nf90_inq_varid ( f % ncid, 'F', varid ) )
    call check ( nf90_get_var ( f % ncid, varid, f % corpar, &
         start = [ 1, 1, 2 ], &
         count = [ 1, size ( f % corpar ), 1 ] ) )

    call check( nf90_redef ( f % ncid ) )
    do i = 1, size ( htend_term_names )
       status = nf90_def_var ( f % ncid, trim ( htend_term_names ( i ) ), &
            NF90_FLOAT, dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
    end do
    
    call check( nf90_enddef ( f % ncid ) )

  end function open_wrf_file

  subroutine close_wrf_file ( file )
    type ( wrf_file ) :: file
    call check ( nf90_close ( file % ncid ) )
  end subroutine close_wrf_file

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

    print*,'Reading ',trim(names(1))
    real3d = data ( names ( 1 ) )

    do i = 2, size ( names )
       print*,'Reading ',trim(names(i))
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

end module mod_wrf_file

module mod_type_vars
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
  integer,parameter :: temp=1,& ! names of dynamic variables
                       uwind=2,&
                       vwind=3,&
                       wwind=4,&
                       ght=5,&
                       pressure=6
  integer,parameter :: rthshten=1,& ! names of physics variables
                       rthcuten=2,&
                       rthraten=3,&
                       rthblten=4,&
                       h_diabatic=5,&
                       rushten=6,&
                       rucuten=7,&
                       rublten=8,&
                       rvshten=9,&
                       rvcuten=10,&
                       rvblten=11
  integer,parameter :: surfacePres=1,& ! names of 2d variables
                       mub=2,&
                       mu=3

  integer, parameter :: max_len=100
  private :: max_len

  type dimension
     integer :: varid, &  !Variable ID in the input file
                dimid, &  !Dimension ID in the input file
                length    !Length of dimension
     
     character(len=max_len) :: dim_name, &       !Name in the dimension block
                               var_name, &       !Name in the variable block
                               units, &          !Units
                               unit_name="units" !Attribute name for units

     real,allocatable,dimension(:) :: data       !Dimension data
  end type dimension
  
  type var4d
     integer :: varid                            !ID for a variable 
     character(len=max_len) :: name              !Name of a variable
     real,allocatable,dimension(:,:,:,:) :: data !data for a variable (x,y,z,t)
  end type var4d
  
  type var3d
     integer :: varid                            !ID for a variable 
     character(len=max_len) :: name              !Name of a variable
     real,allocatable,dimension(:,:,:) :: data !data for a variable (x,y,z)
  end type var3d
  
  type terms3d
     type ( var3d ), dimension ( n_terms ) :: term
  end type terms3d

  type vterm
     type(var3d),dimension(3) :: term
  end type vterm

  type vorterms3d
     type(vterm),dimension(n_terms) :: term
  end type vorterms3d

  type var2d
     integer :: varid
     character(len=max_len) :: name
     real,allocatable,dimension(:,:) :: data !(x,y)
  end type var2d
  
end module mod_type_vars

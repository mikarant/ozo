# README

This is a README-file for OZO software. 


## Source code
Source code files of the OZO are located in the _src_ directory. 

1. mod\_common\_subrs.f90  
	Includes calculation subroutines which are common for both omega and zwack-okossi equation, such as approximation of horizontal and vertical derivatives, laplacian etc.
	
2. mod\_const.f90  
	This short module contains only definition of some natural constants.

3. mod\_omega.f90  
	This module contains subroutines related to the solving of the omega equation.
	
4. mod\_poisson\_DFT.f90  
	This module contains subroutines related to the solving of the Poisson's equation in the Zwack-Okossi equation. 
	
6. mod\_zo.f90  
	This module contains the solving of the Zwack-Okossi equation and all other relevant subroutines.
	
7. mod\_time\_step\_loop.f90  
	Contains the main time stepping loop + some input/output routines. 

8. mod\_wrf\_file.f90  
	This module includes netcdf IO routines to create and handle input and output files.
	
9. ozo.f90  
	The main program of the OZO.
	
## Running WRF and OZO
_run\_wrf_ includes instructions to compile and run WRF model to produce input data for OZO.  

_run\_ozo_ includes instructions to compile and run OZO itself. 

## Test data
Test data for testing the functionality of OZO can be found from the Downloads. It contains output of the WRF baroclinic wave simulation from timesteps 117-120.  
Thus, with that data, OZO can be run for two timesteps (h=118-119).

## Namelist

        &PARAM
        infile='/ozo/test/test_WRF.nc',
        outfile='/ozo/test/ozo_output_test.nc',
        time_1=2,time_n=3,
        alfa=0.2,toler=5e-5,
        ny1=4,ny2=2,
        mode='G',
        calc_omegas=.true.
        /

`infile`: Complete path to the input file. 

`outfile`: Complete path to the output file.  

`time_1`: Starting timestep. Note that due to numerical derivatives, you cannot choose starting timestep to one.  

`time_n`: Ending timestep. This has to be N-1, where N is the total number of timesteps.  

`alfa`: Relaxation coefficient in solving of the omega equation. 0.2 for 100 km resolution and 0.1 for 50 km resolution.  

`toler`: Threshold for testing the convergence in the solving of omega equation. By default, 5e-5 is given.  

`ny1, ny2`: Numbers of sub-cycle iterations in the descending and ascending phases of the multigrid cycle, respectively. By default, 4 and 2 for 100 km resolution are given.  

`mode`: Mode of the omega equation. Choose "G" for generalized one and "Q" for quasi-geostrophic one. "T" and "t" are test modes for both equations.  

`calc_omegas`: True, if you want to calculate vertical motion fields. False, if you have already calculated them, and want now recalculate only height tendencies.  
In the latter case, omegas are read from the output file.

`calc_div`: Choose if you want to calculate omega and height tendency fields due to divergent vorticity and thermal advection.

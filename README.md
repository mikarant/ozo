# README

This is a README-file for OZO software. 


## Source code
Source code files of the OZO are located in the directory called _src_. 

1. mod\_common\_subrs.f90  
	Includes calculation subroutines which are common for solving both omega and zwack-okossi equation, such as horizontal and vertical derivatives, laplacian and relative vorticity.
	
2. mod\_const.f90  
	This short module contains only definition of some natural constants.

3. mod\_omega.f90  
	This file contains the main routine of the solving of omega equation. All the other subroutines called by this routine are located in a file mod\_omega\_subrs.f90.

4. mod\_omega\_subrs.f90  
	All of the omega-related subroutines are located here.
	
5. mod\_poisson\_DFT.f90  
	This module contains subroutines related to solving of Poisson's equation in the Zwack-Okossi equation. 
	
6. mod\_subrs.f90  
	In this module, all the subroutines related to Zwack-Okossi equation are located.
	
7. mod\_time\_step\_loop.f90  
	Time stepping loop + input/output routines are described here. 

8. mod\_wrf\_file.f90  
	This module includes routines to create and handle input and output netcdf files.
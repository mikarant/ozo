# OZO program README

This is a manual for downloading, compiling and running the OZO program in your own computer.

## 1. Technical requirements:

1. Standard NETCDF-library

2. Intel Math Kernel Library. The program uses Intel's MKL library for solving Poisson's equation. 
It can be downloaded for free, but registration is required. 
   See more: [https://software.intel.com/en-us/articles/free_mkl] [MKL]
   
3. GNU's gfortran compiler

## 2. Downloading the source code

1. Launch a terminal window

2. Go to the local directory where you want to put the program

3. Write to the command line:

        git clone git@bitbucket.org:mikarant/ozo.git`

If the clone was successful, you should now have ozo-directory appeared on your local drive.

## 3. Downloading test data

For running the test case, you need to download test data. Datafile is in nc-format and contains WRF-output variables and calculated vertical motion fields from two timesteps.

1. Go to test-directory:

        cd ozo/test

2. Download the data by command (provided that you are added as an user to the private repository) (note this is all one line):

        wget --user=<email> --password=<password>  
        https://bitbucket.org/mikarant/ozo/downloads/wrf_4.nc


## 4. Compiling the program

Downloaded directory contains a makefile for compiling and running the program. At first, you should change paths for Netcdf- and MKL libraries.

1. Go to ozo-directory

2. Open makefile, for example with emacs:

        emacs makefile

3. Change following paths according to where netcdf and mkl libraries are located locally in your computer:

        NETCDF_INCLUDES = -I/usr/include  
        NETCDF_LIBS     = -L/usr/lib -lnetcdff  
        MKLROOT         = /home/mikarant/intel/compilers_and_libraries_2016.2.181/linux/mkl  

4. Save your changes and close the makefile

If your changed paths are correct, you should be now able to compile the program. Compiling can be done by just writing command
```bash
make
```
in the ozo-directory. If the compiling was successful, you should have executable called ``` ozo ``` appeared in your directory.

## 5. Running the test case

Once you have compiled the program, you can test whether it is working by running the test case. You can do it by writing command
```bash
make test
```
in the ozo-directory.


[//]: # (Reference links)

[MKL]: <https://software.intel.com/en-us/articles/free_mkl>


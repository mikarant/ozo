# OZO program README

This is a manual for downloading, compiling and running the OZO program in your own computer.

## 1. Technical requirements:

1. Standard NETCDF-library

2. Intel Math Kernel Library.  
The program uses Intel's MKL library for solving Poisson's equation. It can be downloaded for free, but registration is required. See more: [https://software.intel.com/en-us/intel-mkl] [MKL]. The package includes an installation script which is recommended to use.
   
3. GNU's gfortran compiler

## 2. Downloading the source code

1. Launch a terminal window

2. Go to the local directory where you want to put the program

3. Write to the command line:

        git clone https://mikarant@bitbucket.org/mikarant/ozo.git

 That command clones the repository for your own computer.

4. Extract the tarball by command

        tar -zxvf master.tar.gz


If the download was succesful, you should now have folder called `mikarant-ozo-xxxx` appeared to your local drive. In the folder name, xxxx refers the last commit ID.

## 3. Downloading test data

For running the test case, you need to download a test data. Data file is in netcdf-format and contains WRF-output at hours 117-120 of the idealized baroclinic wave simulation.

1. Go to test-directory:

        cd ozo/test


2. Download the data by command:

        wget https://bitbucket.org/mikarant/ozo/downloads/test_WRF.nc


## 4. Compiling the program

Downloaded directory contains a makefile for compiling and running the program. At first, you should change paths for Netcdf- and MKL libraries.

1. Go back to the ozo-directory

2. Open makefile, for example with emacs:

        emacs makefile

3. Change following paths according to where netcdf and mkl libraries are located locally in your computer. If you are using newer version of the MKL, change also the version number from the path.

        NETCDF_INCLUDES = -I/usr/include  
        NETCDF_LIBS     = -L/usr/lib -lnetcdff  
        MKLROOT         = /home/mikarant/intel/compilers_and_libraries_2016.2.181/linux/mkl  

4. Open _deps.mk_ and change also the paths of the two MKL objects used in the program:

        emacs deps.mk
        mkl_poisson.o : /home/mikarant/intel/compilers_and_libraries_2016.2.181/linux/mkl/include/mkl_poisson.f90 mkl_dfti.o 
        mkl_dfti.o : /home/mikarant/intel/compilers_and_libraries_2016.2.181/linux/mkl/include/mkl_dfti.f90 

5. Save your changes and close the makefile and deps.mk

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

[MKL]: <https://software.intel.com/en-us/intel-mkl>


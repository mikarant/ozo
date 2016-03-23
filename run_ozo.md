# Running Omega - Zwack-Okossi program

This is a manual for downloading, compiling and running the OZO-program.

Technical requirements:

1. Standard NETCDF-library

2. Intel Math Kernel Library. The program uses Intel's MKL library for solving poisson's equation. It can be downloaded for free, but registration is required. 

See more: [https://software.intel.com/en-us/articles/free_mkl] [MKL]

## 1. Downloading the source code


1. Launch a terminal window
2. Go to the local directory where you want to put the program
3. Write to the command line:
```sh
git clone git@bitbucket.org:mikarant/ozo.git
```
If the clone was successful, you should now have ozo-directory appeared on your local drive.

## 2. Compiling the program

Downloaded directory contains a makefile for compiling and running the program. At first, you should change paths for netcdf- and mkl-libraries.


[//]: # Reference links

[MKL]: <https://software.intel.com/en-us/articles/free_mkl>
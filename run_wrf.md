# WRF

## Download

Here you can find different versions of the WRF model:

[http://www2.mmm.ucar.edu/wrf/users/download/get_sources.html] [WRF]

## Tutorial

And here is information about downloading and compiling WRF on your local computer

[http://www2.mmm.ucar.edu/wrf/OnLineTutorial/index.htm] [tutorial]

## Instructions to produce input data to OZO by WRF idealized simulation

### 1. Downloading and compiling WRF on your local computer
```sh
# Build
wget http://www2.mmm.ucar.edu/wrf/src/WRFV3.7.1.TAR.gz
tar xvf WRFV3.7.1.TAR.gz
cd WRFV3/
export NETCDF=/usr WRFIO_NCD_LARGE_FILE_SUPPORT=1
# GCC, serial, no nesting
./configure <<<'32

'
export WRF_EM_CORE=1
./compile em_b_wave >& compile.log
```
Check the compile.log file for any errors.

### 2. Running the idealized baroclinic wave simulation in WRF
First, you have to edit the namelist to set correct values of some parameters.
```sh
cd test/em_b_wave/
```
Open _namelist.input_ with some text editor and change at least following values:  
Section _&time\_control_:
`run_days = 10` Sets simulation time to 10 days.  
`end_day = 10` Ending day of the simulation.  
`history_interval = 60` Output interval in seconds.  
`iofields_filename = = "iofield_list.txt"` This is optional. By default, WRF outputs huge number of unnecessary variables. With that file you can change the number of output variables. By default, there is no place for this in the namelist, but it can be added there to the _&time\_control_ section.

Section &physics:
`mp_physics = 2`  
`sf_sfclay_physics = 1`  
`sf_surface_physics = 1`  
`bl_pbl_physics = 1`  
`cu_physics = 1`  


```sh

./run_me_first.csh
ln -s ../../run/LANDUSE.TBL


./ideal.exe
./wrf.exe
```




[//]: # (Reference links)

[WRF]: <https://software.intel.com/en-us/articles/free_mkl>
[tutorial]: <http://www2.mmm.ucar.edu/wrf/OnLineTutorial/index.htm>


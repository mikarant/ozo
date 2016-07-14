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

### 2. Editing the namelist 
Go to the baroclinic wave directory:  
```sh
cd test/em_b_wave/
```
Open _namelist.input_ with some text editor and change at least following values:  

Section _&time\_control_:  


`run_days = 10` Sets simulation time to 10 days. 

`end_day = 10` Ending day of the simulation.  

`history_interval = 60` Output interval in seconds. 

`iofields_filename = = "iofield_list.txt"` This is optional. By default, WRF outputs huge number of unnecessary variables.  
With that file you can change the number of output variables. There is no place for this in the namelist by default, but it can be added there to the _&time\_control_ section.


Section &physics:  


`mp_physics = 2`  

`sf_sfclay_physics = 1`  

`sf_surface_physics = 1`  

`bl_pbl_physics = 1`  

`cu_physics = 1`  


### 3. Running the model
Once you have set correct values in the namelist, you have to link some files to running directory.  
This can be done by executing _run\_me\_first.csh_ script and linking one other file:

```sh
./run_me_first.csh
ln -s ../../run/LANDUSE.TBL
```

Now you are ready to create initial state of the simulation by running _ideal.exe_:

```sh
./ideal.exe
```

If the file _wrfinput\_d01_ appears to directory, your initial state is created.  
After that, start running the model:

```sh
./wrf.exe
```




[//]: # (Reference links)

[WRF]: <https://software.intel.com/en-us/articles/free_mkl>
[tutorial]: <http://www2.mmm.ucar.edu/wrf/OnLineTutorial/index.htm>


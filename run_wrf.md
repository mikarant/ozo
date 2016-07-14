# WRF

## Download

Here you can find different versions of the WRF model:

[http://www2.mmm.ucar.edu/wrf/users/download/get_sources.html] [WRF]

## Tutorial

And here is information about downloading and compiling WRF on your local computer

[http://www2.mmm.ucar.edu/wrf/OnLineTutorial/index.htm] [tutorial]

## Downloading and compiling WRF on your local computer
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

## Running the idealized baroclinic wave simulation in WRF
```sh
cd test/em_b_wave/
./run_me_first.csh
ln -s ../../run/LANDUSE.TBL
./ideal.exe
./wrf.exe
```




[//]: # (Reference links)

[WRF]: <https://software.intel.com/en-us/articles/free_mkl>
[tutorial]: <http://www2.mmm.ucar.edu/wrf/OnLineTutorial/index.htm>


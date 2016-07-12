# WRF

## Download

[http://www2.mmm.ucar.edu/wrf/users/download/get_sources.html] [WRF]

## Tutorial

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

# Test
cd test/em_b_wave/
./run_me_first.csh
ln -s ../../run/LANDUSE.TBL
./ideal.exe
./wrf.exe
```
[//]: # (Reference links)

[WRF]: <https://software.intel.com/en-us/articles/free_mkl>


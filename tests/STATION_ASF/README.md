# *Station Air-Sea Fluxes* demonstration case

Last successful test done with NEMOGCM trunk: `r13263`

Author: Laurent Brodeau, 2020

NOTE: if working with the trunk of NEMO, you are strongly advised to use the same test-case but on the `NEMO-examples` GitHub depo:
https://github.com/NEMO-ocean/NEMO-examples/tree/master/STATION_ASF

## Objectives

```STATION_ASF``` is a demonstration test-case that mimics a (static) in-situ station (buoy, platform) dedicated to the estimation of surface air-sea fluxes by means of *widely-measured* (bulk) meteorological surface parameters.

```STATION_ASF``` has been constructed by merging the *single column* and the *standalone surface module* configurations of NEMO. In short, it can be defined as "SAS meets C1D". As such, the spatial domain of ```STATION_ASF``` is punctual (1D, well actually 3 x 3 as in C1D).

```STATION_ASF``` is therefore a versatile tool, and extremely lightweight in terms of computing requirements, to test the different bulk algorithms and cool-skin/warm-layer parameterization options included in NEMO.

As input ```STATION_ASF``` will require the traditional *bulk* sea surface parameters:

- Bulk sea surface temperature (SST) at _z<sub>SST</sub>_ meters below the surface
- Surface current vector
- Sea surface salinity

as well as the usual surface atmospheric state:

- air temperature at _z<sub>t</sub>_ meters above the surface
- air humidity  at _z<sub>t</sub>_ meters above the surface (specific humidity or relative humidity or dew-point temperature)
- wind speed vector at _z<sub>u</sub>_ meters above the surface
- Sea level atmospheric pressure (SLP)
- Downwelling solar radiation
- Downwelling longwave radiation

### Example of diagnostics from `STATION_ASF`

(Generated with script `./EXPREF/plot_station_asf_simple.py`)

![plot](https://github.com/NEMO-ocean/NEMO-examples/blob/master/STATION_ASF/figs/01_temperatures_ECMWF.svg)

![plot](https://github.com/NEMO-ocean/NEMO-examples/blob/master/STATION_ASF/figs/Cd.svg)

![plot](https://github.com/NEMO-ocean/NEMO-examples/blob/master/STATION_ASF/figs/dT_skin.svg)

![plot](https://github.com/NEMO-ocean/NEMO-examples/blob/master/STATION_ASF/figs/Qlat.svg)


## Physical description

### Important namelist parameters specific to STATION_ASF

* ```rn_dept1@namusr_def:``` depth (m) at which the prescribed SST is taken (*i.e.* depth of first T-point); important due to impact on warm-layer estimate, the deeper, the more pronounced!

* ```rn_lat1d,rn_lon1d@namc1d:``` fixed coordinates of the location of the station (buoy, platform, etc).

* ```namsbc_blk:``` to be filled carefully, just as for "C1D", the prescribed surface ATMOSPHERIC state (files) are time series of shape 3x3 in space

* ```namsbc_sas:``` to be filled carefully, just as for "C1D", the prescribed surface OCEAN state (files) are time series of shape 3x3 in space



## Input files to test STATION ASF

One full year (2018) of processed hourly data from the PAPA station (buoy) is found into the `input_data` directory.
These three files are everything you need to play with the set of *namelists* provided for this test-case.

- ```Station_PAPA_50N-145W_atm_hourly_y2018.nc```  → contains hourly surface atmospheric state
- ```Station_PAPA_50N-145W_precip_daily_y2018.nc``` → contains daily precipitation
- ```Station_PAPA_50N-145W_oce_hourly_y2018.nc``` → contains hourly sea surface state

For station PAPA (50.1 N, 144.9 W), air temperature and humidity are measured at 2.5 m, the wind speed at 4 m, and the SST at 1 m below the surface, hence the following namelist parameters are given:

- `&namusr_def`
  - ```rn_dept1 =    1.  ```
- `&namc1d`
  - ```rn_lat1d =  50.1 ```
  - ```rn_lon1d = 215.1```
- `&namsbc_blk`
  - ```rn_zqt   =   2.5```
  - ```rn_zu    =    4.```



## Playing with STATION_ASF

First compile the test-case as follows (compile with xios-2.5 support → check your ARCH file):

```./makenemo -a STATION_ASF -m <your_arch> -n STATION_ASF2 -j 4```

Then you can use the script ``launch_sasf.sh`` found in  ```EXPREF/``` to launch 3 simulations (one for each bulk parameterization available). You need to adapt the following variable to your environment in the script:

- ```NEMO_ROOT_DIR``` : NEMO root directory where to fetch compiled STATION_ASF ```nemo.exe``` + setup (such as ```${NEMO_ROOT_DIR}/tests/STATION_ASF```)

- ```PROD_DIR``` :  Directory where to run the simulation

- ```DATA_IN_DIR``` : Directory containing sea-surface + atmospheric forcings (found here in ```input_data/```)

If everything goes according to plan, ``launch_sasf.sh`` should have generated the 3 following sets of output files into `${PROD_DIR}/output`:

    STATION_ASF-COARE3p6_1h_20180101_20181231_gridT.nc
    STATION_ASF-COARE3p6_1h_20180101_20181231_gridU.nc 
    STATION_ASF-COARE3p6_1h_20180101_20181231_gridV.nc 
    STATION_ASF-ECMWF_1h_20180101_20181231_gridT.nc 
    STATION_ASF-ECMWF_1h_20180101_20181231_gridU.nc 
    STATION_ASF-ECMWF_1h_20180101_20181231_gridV.nc 
    STATION_ASF-NCAR_1h_20180101_20181231_gridT.nc 
    STATION_ASF-NCAR_1h_20180101_20181231_gridU.nc 
    STATION_ASF-NCAR_1h_20180101_20181231_gridV.nc

---

*/Laurent, July 2020.*


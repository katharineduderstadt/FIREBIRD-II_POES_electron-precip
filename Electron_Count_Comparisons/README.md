# FIREBIRD-II_POES_electron-precip


This repository contains the code used to estimate FIREBIRD count equivalents 
based on POES/MetOP MEPED instruments. 

Sample code is for a conjunction between FIRBIRD-II FU3 and MetOP2a on 2018-12-31.
Modify dates and parameters in the shell script  
      sh FB-POES_counts.sh

## Description

This study compares energetic electron precipitation measurements
by the Focused Investigations of Relativistic Electron Burst Intensity
Range, and Dynamics (FIREBIRD-II) CubeSats with NOAA Polar-orbiting
Operational Environmental Satellite (POES) and ESA Meteorological 
Operational satellite (MetOp) satellites, which are equipped with the
Medium-Energy Proton Electron Detector (MEPED).

See "Comparisons of Energetic Electron Observations between 
FIREBIRD-II CubeSats and POES/MetOp Satellites from 2018-2020"
by Householder et al. (submitted 2024)

## Getting Started

sh FB-POES_counts.sh

Modify dates and parameters in the shell script FB-POES_counts.sh
    
Input files for G-factors - see GEANT folder (files created with FB_GEANT_all.py) 

      FU3-col.txt
      FU3_Gfactors.txt

Input files in DATA folders FIREBIRD and POES

      FIREBIRD data downloaded from http://solar.physics.montana.edu/FIREBIRD_II/Data/ 
      FU3_Hires_2018-12-31_L2.txt

      POES data downloaded from https://www.ngdc.noaa.gov/stp/satellite/poes/dataaccess.html 
      poes_m02_20181231_proc.nc
      poes_m02_20181231_raw.nc

      POES corrected fluxes and bounce loss cone estimates provided by Josh Pettit 
      POES_combinedSpectrum_2sec_m02_00_20181231.nc
      POES_combinedSpectrum_2sec_m02_90_20181231.nc
      POES_combinedSpectrum_m02_BLC_highres_20181231.nc

      L-shells for POES data come from 
      poes_m02_20181231_raw_magephem.txt

### Dependencies

python imports include: netcdf4, numpy, scipy, matplotlib, sys, datetime, ast, math, calendar

### Installing

git init

git clone https://github.com/katharineduderstadt/FIREBIRD-II_POES_electron-precip/


### Executing program

Modify dates and parameters in the shell script FB-POES_counts.sh

Ensure all input files are available for dates.  

Execute shell script

    sh FB-POES_counts.sh

The following python codes will run
      FB_select_times.py
      FB_datetime_average.py
      FB_flux_estimate.py
      POES_counts_corrected.py


## Authors

Isabella Householder   isabella.householder@unh.edu 
Katharine Duderstadt   katharine.duderstadt@unh.edu


## Acknowledgments

This research was supported by NSF CEDAR (1650738), NSF 1035642 and NASA (135260, NNX15AF66G). 
FIREBIRD-II data was made possible by the NSF (0838034, 1339414). 







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

## Directories 

Electron_Count_Comparisons
      Contains Python code and associated files to calculate equivalent FIREBIRD counts 
      for comparison with POES/MetOp MEPED. 

Data_Summary
      Contains .csv file with summary data as well as summary plot for L=5.

## DATA

FIREBIRD-II data 
https://solar.physics.montana.edu/FIREBIRD_II/

POES/MetOp data available from 
https://www.ngdc.noaa.gov/stp/satellite/poes/dataaccess.html

Corrected POES data available directly from Josh Pettit
joshua.m.pettit@nasa.gov 

## Authors

Isabella Householder   isabella.householder@unh.edu 

Katharine Duderstadt   katharine.duderstadt@unh.edu


## Acknowledgments

This research was supported by NSF CEDAR (1650738), NSF 1035642 and NASA (135260, NNX15AF66G). 
FIREBIRD-II data was made possible by the NSF (0838034, 1339414).  






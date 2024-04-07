#!/bin/sh
#
# sh FB-POES_counts.sh
#
# Edit lines 13-56 and 81-83 for specific FIREBIRD-POES data
# Note: avoid extra spaces...
#
# This version makes plots from estimated counts files
#
# Modified May 2023 to include new Pettit POES data
# ----------------------------------------------------

FB_dir='DATA/FIREBIRD/'
POES_dir='DATA/POES/'
output_dir='OUTPUT/'
GEANT_dir='G-factors/'

FBdatetime='2018-12-31_0720-0724'
POEStxtdate='2018-12-31_0720-0724'
POESdatetime='20181231'
FBsc='FU3'
POESsc='m02'

FBhires=$FB_dir$FBsc'_Hires_2018-12-31_L2.txt'

startFB='2018,12,31,7,20,00'
endFB='2018,12,31,7,24,00'

#FU4
#energy_FB_width='[63.7,100.2,136.7,200.4,264.3]'
#energy_FB='[251.5,333.5,452.0,620.5,852.8]'
#energy_legend='[252_keV,334_keV,452_keV,621_keV,853_keV]'

#FU3
energy_FB_width='[68.7,107.9,147.2,215.9,284.5]'
#energy_FB='[265.4,353.7,481.2,662.7,913.0]'
#energy_legend='[265_keV,354_keV,481_keV,663_keV,913_keV]'

energy_FB='[200.,300.,450.,600.,1000.]'
energy_legend='[200_keV,300_keV,450_keV,600_keV,1000_keV]'


startPOES='2018,12,31,7,21,00'
endPOES='2018,12,31,7,28,00'

timePLOTstart='2018,12,31,7,20,00'
timePLOTend='2018,12,31,7,28,00'

conj_time='2018,12,31,7,22,30'  # L 6.5

FBtxtlat='lat_-57_to_-72,_lon_-24_to_-56'
POEStxtlat='xxx'
FBtxtlat_conj='xxx'
POEStxtlat_conj='xxx'


# -----------------------------------------------------

POEScount=$POES_dir'poes_'$POESsc'_'$POESdatetime'_raw.nc'
#POEScorrect=$POES_dir'POES_combinedSpectrum_'$POESsc'_00_'$POESdatetime'.nc'
POES_Lshell=$POES_dir'poes_'$POESsc'_'$POESdatetime'_raw_magephem.txt'

FBslice=$output_dir'TimeSlice_'$FBsc'_'$FBdatetime'.txt'
FBavg=$output_dir'sec-avg_'$FBsc'_'$FBdatetime'.txt'
Gfactors=$GEANT_dir$FBsc'_Gfactors.txt'
GEANT=$GEANT_dir$FBsc'-col.txt'
fluxpar=$output_dir'flux-params_'$FBsc'_'$FBdatetime'.txt'

outvalues=$output_dir'Values_'$FBsc'_'$FBdatetime'_'$POESsc'.txt'
outtxt=$output_dir'FB_equiv_counts_'$FBsc'_'$FBdatetime'.txt'
outtxt_0tel=$output_dir'POES_counts_0tel'_''$POESsc'_'$FBdatetime'.txt'
outtxt_90tel=$output_dir'POES_counts_90tel'_''$POESsc'_'$FBdatetime'.txt'

outplot_L=$output_dir'FU-POES_counts_'$FBsc'_'$FBdatetime'_'$POESsc'_L.ps'
outplot_time=$output_dir'FU-POES_counts_'$FBsc'_'$FBdatetime'_'$POESsc'_time.ps'
title_plot=$FBsc'_and_'$POESsc'_conjunction'

FBtxt=$FBsc'_'$FBdatetime
POEStxt=$POESsc'_'$POEStxtdate

POES_Pettit_file_0=$POES_dir'POES_combinedSpectrum_2sec_m02_00_20181231.nc'  #Josh Pettit POES May 2023
POES_Pettit_file_90=$POES_dir'POES_combinedSpectrum_2sec_m02_90_20181231.nc'  #Josh Pettit POES May 2023
POES_Pettit_file_BLC=$POES_dir'POES_combinedSpectrum_m02_BLC_highres_20181231.nc'  #Josh Pettit POES May 2023
# ------------------------------------------------------
#To re-run this code to only plot (and not repeat other steps), comment out the first three python calls
#and only keep python POES_counts_corrected.py line

echo "Run FB_select_times_shell.py"
python FB_select_times.py $FBhires $FBslice $startFB $endFB
# ---------------------------------------------------
echo "Run FB_datetime_average_shell.py"
python FB_datetime_average.py $FBslice $FBavg
# ----------------------------------------------------
echo "Run FB_flux_estimate_shell.py"
echo "This one takes a few minutes"
python FB_flux_estimate.py $Gfactors $FBavg $GEANT $fluxpar $energy_FB_width $energy_FB
# ----------------------------------------------------
echo "POES_counts_corrected.py"
python POES_counts_corrected.py $POESsc $FBsc $POEScount $POES_Lshell $fluxpar $FBavg $outvalues $outtxt $outtxt_0tel $outtxt_90tel $outplot_L $outplot_time $title_plot $conj_time $timePLOTstart $timePLOTend $startPOES $endPOES $FBtxt $POEStxt $FBtxtlat $FBtxtlat_conj $POEStxtlat $POEStxtlat_conj $energy_FB $energy_legend $POES_Pettit_file_0 $POES_Pettit_file_90 $POES_Pettit_file_BLC  #Josh Pettit POES May 2023
#-----------------------------------------------------

exit

'''
Filename: 
    POES_counts_corrected.py
    
Author: 
    Isabella M. Householder
    Katharine A. Duderstadt
    
Synopsis:
    Use electron flux calculated from FIREBIRD counts to calculate counts the POES instruments
    should measure using Yando GEANT tables
    
Date:
    20 Nov 2019
    last modified 1 Apr 2023 to prepare for archive
    
Description:
    Read in Yando POES electron GEANT table
          Yando_TableB2.csv - lowest energy in bin, channel width (keV), E1 Gfactor, E2 Gfactor, E3 Gfactor
          Yando_TableB4.csv - lowest energy in bin, channel width (keV), P6 Gfactor
          values have been divided be 100 to arrive at correct units  cm^2*sr
    Read in FIREBIRD environmental flux values
          in electrons keV-1 cm-2 s-1 sr-1
    Read in POES observed counts
    Calculated POES counts 
         Interpolate exponential fit to flux to yando energies
         counts_bin(E1, E2, E3) = (flux in bin)*G(E1, E2, E3)
         counts_channel(E1, E2, E3) - Integrate all energy bins for each channel
    Plot calculated vs observed POES counts for a selected time step
    
    Note:
    The parameters J0 and E0 are from curve fit in "FB_flux_estimate.py" that fits an
    exponential function to Energy (x-values) in MeV and Flux (y-values) in 
    electrons cm^-2 s^-1 st^-1 keV^-1. This means that when using the parameters 
    to calculate flux for a given energy, the energy must be specified in MeV 
    while the return flux is uses keV. (We find this approach necessary in order
    for the function "curve_fit" to estimate parameters without crashing.)

    Updated August 2022 to read in Josh Pettit's files
    
    
Input files:
    Yando_TableB2.csv
    Yando_TableB4.csv
    flux-params_FU3_2018-12-31_0720-0724.txt
    sec-avg_FU3_2018-12-31_0720-0724.txt
    poes_m01_20180922_raw.nc
    poes_m01_20180922_proc.nc
    POES_combinedSpectrum_2sec_m02_00_20181231.nc
    POES_combinedSpectrum_2sec_m02_90_20181231.nc
    POES_combinedSpectrum_m02_BLC_highres_20181231.nc
    
Output files: 
    FU-POES_counts_FU3_2018-12-31_0720-0724_m02_L_ALL.ps
    FU-POES_counts_FU3_2018-12-31_0720-0724_m02_L.ps
    FU-POES_counts_FU3_2018-12-31_0720-0724_m02_time.ps
    FB_equiv_counts_FU3_2018-12-31_0720-0724.txt
    POES_counts_0tel_FU3_2018-12-31_0720-0724.txt
    POES_counts_90tel_FU3_2018-12-31_0720-0724.txt
    POES_counts_geometric_mean.txt
    Values_FU3_2018-12-31_0720-0724_m02.txt

To Run:
    python POES_counts_corrected.py
    
Modify for new POES data from Josh Pettit on May 2023

-------------------------------------------------------------------------------
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.dates as mdates
import datetime as dt
from datetime import datetime, timedelta, timezone
import calendar
import netCDF4 as nc
from netCDF4 import Dataset
import sys
#import pytz
import scipy
from scipy.interpolate import interp1d

# ------------------------------------------------------------------------------
# UPDATE THIS SECTION FOR EACH CONJUNCTION
# You'll need to use https://sscweb.gsfc.nasa.gov/ to find
#         lat and lon for the range of FB and POES data (locator)
#         conjunctions times and lat and lon of conjunction(query)
#
# To find the L-shells associated with the conjunction time:
# FB L-shells are in the OUTPUT/sec-avg_FU4_2018-09-28_0031-0036.txt (the last value in each line)
# POES L-shells are in  /DATA/POES/poes_n18_20180928_raw_magephem.txt - the second value (time, L, MLT)
# ------------------------------------------------------------------------------

POESsc=sys.argv[1]                   #$POESsc
FBsc=sys.argv[2]                     #$FBsc

filename_POES_counts=sys.argv[3]     #$POEScount
filename_POES_L = sys.argv[4]        #$POES_Lshell
#filename_correct=sys.argv[6]               #$POEScorrect

filename_FB_flux=sys.argv[5]         #$fluxpar
filename_FB_Lshell=sys.argv[6]       #$FBavg
filename_values=sys.argv[7]          #$outvalues

filename_output=sys.argv[8]          #$outtxt
filename_output_POES_0tel=sys.argv[9]         #$outtxt_0tel
filename_output_POES_90tel=sys.argv[10]       #$outtxt_90tel

filename_plot_L=sys.argv[11]         #$outplot_L
filename_plot_L_all=filename_plot_L[:-3]+'_ALL.ps'
filename_plot_time=sys.argv[12]      #$outplot_time

title_plot = sys.argv[13]            #$title_plot
conj = sys.argv[14]                  #$conj_time

startplot=sys.argv[15]               #$timePLOTstart
endplot=sys.argv[16]                 #$timePLOTend
startarg=sys.argv[17]                #$startPOES
endarg=sys.argv[18]                  #$endPOES

FBfig_text=sys.argv[19]              #$FBtxt
POESfig_text=sys.argv[20]            #$POEStxt
FBfig_text_latlon=sys.argv[21]       #$FBtxtlat
FBfig_txt_conj = sys.argv[22]        #$FBtxtlat_conj
POESfig_text_latlon=sys.argv[23]     #$POEStxtlat
POESfig_txt_conj = sys.argv[24]      #$POEStxtlat_conj

FB_energy_middles=sys.argv[25]       #$energy_FB
legend_string=sys.argv[26]           #$energy_legend

filename_Pettit_0_nc=sys.argv[27]      #$POES_Pettit_file -  May 2023
filename_Pettit_90_nc=sys.argv[28]      #$POES_Pettit_file -  May 2023
filename_Pettit_BLC_nc=sys.argv[29]      #$POES_Pettit_file -  May 2023


sd=[int(s) for s in startarg.split(',')]
ed=[int(s) for s in endarg.split(',')]
sd2=[int(s) for s in startplot.split(',')]
ed2=[int(s) for s in endplot.split(',')]
sd3=[int(s) for s in conj.split(',')]
ed3=[int(s) for s in conj.split(',')]
startPOES = dt.datetime(sd[0],sd[1],sd[2],sd[3],sd[4],sd[5])
endPOES = dt.datetime(ed[0],ed[1],ed[2],ed[3],ed[4],ed[5])
starttime_plot = dt.datetime(sd2[0],sd2[1],sd2[2],sd2[3],sd2[4],sd2[5])
endtime_plot = dt.datetime(ed2[0],ed2[1],ed2[2],ed2[3],ed2[4],ed2[5])
conj_time = dt.datetime(ed3[0],ed3[1],ed3[2],ed3[3],ed3[4],ed3[5])

energy_FB_list = [float(s) for s in FB_energy_middles.strip('[]').split(',')]
energy_FB = np.array(energy_FB_list, dtype=float)
energy_legend = legend_string.strip('[]').split(',')

'''
POESsc ='noaa18'
FBsc = 'FU4'

filename_POES_counts = 'DATA/POES/poes_n18_20180927_raw.nc'
filename_POES_L = 'DATA/POES/poes_n18_20180927_raw_magephem.txt'
#filename_correct = 'DATA/POES/POES_combinedSpectrum_n18_00_20180928.nc'  #NOT USED YET

filename_FB_flux = 'OUTPUT/flux-params_FU4_2018-09-27_0135-0142.txt'  #exiting
filename_FB_Lshell = 'OUTPUT/sec-avg_FU4_2018-09-27_0135-0142.txt'    #existing
filename_values = 'OUTPUT/Values_n18_20180927.txt'                 #exiting

filename_output = 'OUTPUT/FB_equiv_counts_FU4_2018-09-27_0135-0142.txt'     #new
filename_output_POES_0tel = 'OUTPUT/POES_counts_0tel_FU4_2018-09-27_0135-0142.txt'  #new
filename_output_POES_90tel = 'OUTPUT/POES_counts_90tel_FU4_2018-09-27_0135-0142.txt' #new

filename_plot_L = 'OUTPUT/Plot_FU4_noaa18_2018-09-27_0135-0142_Lshell.ps'  #new
filename_plot_time = 'OUTPUT/Plot_FU4_noaa18_2018-09-27_0135-0142_time.ps'  #new

title_plot = 'FU4 and NOAA-18 conjunction'
conj_time = dt.datetime(2018,9,27,1,38,00)   # if a range, choose something in the middle

starttime_plot=dt.datetime(2018,9,27,1,35,00)  #Use full extent of FB adn POES times
endtime_plot=dt.datetime(2018,9,27,1,44,00)    #Use full extent of FB adn POES times

startPOES=dt.datetime(2018,9,27,1,36,00)       #POES times to plot 3-7
endPOES=dt.datetime(2018,9,27,1,44,00)         #POES times to plot 3-7

FBfig_text = 'FU4 2018-09-27 00:35-00:39'         # FIREBIRD hi-res download range
POESfig_text = 'noaa18 2018-09-27 00:36-00:42'     # POES time range
FBfig_text_latlon = 'lat x, lon x'  # latitude for FB data
FBfig_txt_conj = 'conj lat x, lon x, L x'      #FB lat, lon, and L of actual conjunction (or range)
POESfig_text_latlon = 'lat x to x1, lon x to x'
POESfig_txt_conj = 'conj lat -x lon x, L x'    #POES lat, lon, and L of actual conjunction (or range)

energy_FB = [252., 333., 452., 621., 853.]     # for FIREBIRD Unit and campaign
energy_legend = ['252 keV', '333 keV', '452 keV', '621 keV','853 keV']  # for FIREBIRD Unit & campaign
'''

horz_size = 8.5       # horizontal size of plots in inches (8-10 is good)    #MLT
#vert_size_L = 11.       # vertical size of plots in inches (6-8 is good)     #MLT
vert_size_L = 13.       # vertical size of plots in inches (6-8 is good)     #MLT
#vert_size_time = 11.   #MLT
vert_size_time = 13.   #MLT

# -----------------------------------------------------------------------------
# END UPDATE
# -----------------------------------------------------------------------------

filename_POES_GfactorsE = 'G-factors/Yando_TableB4.csv'
filename_POES_GfactorsP = 'G-factors/Yando_TableB2.csv'


def read_POES_Gfactors(filename_GE, filename_GP):
    poesE_GEANT = np.genfromtxt(filename_GE, delimiter=',')
    poesP_GEANT = np.genfromtxt(filename_GP, delimiter=',')
    return poesE_GEANT, poesP_GEANT

def read_FB_flux(filename_FB, filename_FB_L):
    timestamp_list = []
    J0_value = []
    E0_value = []
    L_value = []
    MLT_value = []   #MLT
    file = open(filename_FB, 'r')
    nline = 0
    for line in file:
        line = line.strip().split(',')
        datetime_str = line[0]
        timestamp_list.append(dt.datetime.strptime(datetime_str, '%Y-%m-%dT%H:%M:%S.%fZ'))
        J0_value.append(float(line[1]))
        E0_value.append(float(line[2]))
        nline = nline+1
    J0_value = np.array(J0_value)
    E0_value = np.array(E0_value)
    file.close()
    file = open(filename_FB_L, 'r')
    for line in file:
        line = line.strip().split(',')
        L_value.append(float(line[7]))
        MLT_value.append(float(line[8]))   #MLT
    L_value = np.array(L_value)
    MLT_value = np.array(MLT_value)    #MLT
    #print('FB L-values',L_value)
    file.close()
    return timestamp_list, J0_value, E0_value, L_value, MLT_value, nline   #MLT

def read_POES_counts(filename_raw_nc,nlines):
    nc_fid = Dataset(filename_raw_nc, "r")
    e1_cps = nc_fid.variables['mep_ele_tel0_cps_e1']
    e2_cps = nc_fid.variables['mep_ele_tel0_cps_e2']
    e3_cps = nc_fid.variables['mep_ele_tel0_cps_e3']
    p6_cps = nc_fid.variables['mep_pro_tel0_cps_p6']
    e1_90_cps = nc_fid.variables['mep_ele_tel90_cps_e1']
    e2_90_cps = nc_fid.variables['mep_ele_tel90_cps_e2']
    e3_90_cps = nc_fid.variables['mep_ele_tel90_cps_e3']
    p6_90_cps = nc_fid.variables['mep_pro_tel90_cps_p6']
    timeindex = len(e1_90_cps)
    year = nc_fid.variables['year']
    day = nc_fid.variables['day']
    print("formatting time for POES...this takes a couple minutes")
    timeclock = [dt.datetime.strptime(f'{nc_fid["year"][i]}, {nc_fid["day"][i]}', '%Y, %j') +
                 timedelta(milliseconds=int(nc_fid['msec'][i]))
                 for i in range(len(nc_fid['time'][:]))]
    e1_counts = np.zeros((timeindex),'f')
    e2_counts = np.zeros((timeindex),'f')
    e3_counts = np.zeros((timeindex),'f')
    p6_counts = np.zeros((timeindex),'f')
    e1_90_counts = np.zeros((timeindex),'f')
    e2_90_counts = np.zeros((timeindex),'f')
    e3_90_counts = np.zeros((timeindex),'f')
    p6_90_counts = np.zeros((timeindex),'f')
    e1_Gmean_counts = np.zeros((timeindex),'f')
    e2_Gmean_counts = np.zeros((timeindex),'f')
    e3_Gmean_counts = np.zeros((timeindex),'f')
    p6_Gmean_counts = np.zeros((timeindex),'f')
    e1_counts[:] = e1_cps[:]
    e2_counts[:] = e2_cps[:]
    e3_counts[:] = e3_cps[:]
    p6_counts[:] = p6_cps[:]
    e1_90_counts[:] = e1_90_cps[:]
    e2_90_counts[:] = e2_90_cps[:]
    e3_90_counts[:] = e3_90_cps[:]
    p6_90_counts[:] = p6_90_cps[:]
    e1_Gmean_counts[:] = np.sqrt(e1_cps[:]* e1_90_cps[:])
    e2_Gmean_counts[:] = np.sqrt(e2_cps[:]* e2_90_cps[:])
    e3_Gmean_counts[:] = np.sqrt(e3_cps[:]* e3_90_cps[:])
    p6_Gmean_counts[:] = np.sqrt(p6_cps[:]* p6_90_cps[:])

    POES_obs = np.zeros((timeindex,4),'f')
    POES_obs_90 = np.zeros((timeindex,4),'f')
    POES_obs_Gmean = np.zeros((timeindex,4),'f')
    POES_obs[:,0] = e1_counts
    POES_obs[:,1] = e2_counts
    POES_obs[:,2] = e3_counts
    POES_obs[:,3] = p6_counts
    POES_obs_90[:,0] = e1_90_counts
    POES_obs_90[:,1] = e2_90_counts
    POES_obs_90[:,2] = e3_90_counts
    POES_obs_90[:,3] = p6_90_counts
    POES_obs_Gmean[:,0] = e1_Gmean_counts
    POES_obs_Gmean[:,1] = e2_Gmean_counts
    POES_obs_Gmean[:,2] = e3_Gmean_counts
    POES_obs_Gmean[:,3] = p6_Gmean_counts
    POES_obs_Gmean[POES_obs_Gmean==0]=['nan']
    '''
     for i in range(3):
        for n in range(nlines):
            if POES_obs_Gmean[n,i] <= 0.:
                POES_obs_Gmean[n,i] = 'nan'
     '''
    
    nc_fid.close()
    return timeclock, POES_obs, POES_obs_90, POES_obs_Gmean, timeindex


# -----------------------------------------------------------------------------
# Read Josh Pettit's corrected POES and BLC flux files.  Set all missing values to zero
# to be consistent with the POES raw datafile.
#
# Modified May 2023 for new POES data from Josh Pettit
# Zeros reset to 0.01 for geometric mean calculations
# -----------------------------------------------------------------------------

def read_Pettit_POES_counts(filename_Pettit_0, filename_Pettit_90, filename_Pettit_BLC, POEStime_orig, FB_en, POES_indexstart, POES_indexend):
    timeindex = len(POEStime_orig)
    nc_fid_0 = Dataset(filename_Pettit_0, "r")
    nc_fid_90 = Dataset(filename_Pettit_90, "r")
    nc_fid_BLC = Dataset(filename_Pettit_BLC, "r")
    time_correct = nc_fid_0.variables['time']
    len_time_correct = len(time_correct)
    energy_27 = nc_fid_0.variables['energy']
    E_correct_0_orig = nc_fid_0.variables['EOcounts_corrected']
    E_correct_90_orig = nc_fid_90.variables['EOcounts_corrected']
    
    E_correct_0 = np.zeros((len_time_correct,3),'f')
    E_correct_90 = np.zeros((len_time_correct,3),'f')
    E_correct_0[:,:] = E_correct_0_orig[:,1:4]
    E_correct_90[:,:] = E_correct_90_orig[:,1:4]
    
    BLC_flux = nc_fid_BLC.variables['BLC_Flux']
    #E_correct_0_new = np.zeros((timeindex,3),'f')
    #E_correct_90_new = np.zeros((timeindex,3),'f')
    E_correct_Gmean = np.zeros((timeindex,3),'f')
    #BLCflux_new = np.zeros((timeindex,5),'f')
    
    print("filename_Pettit_0")
    print(filename_Pettit_0)
    print("filename_Pettit_90")
    print(filename_Pettit_90)
    
    #print("Josh_time")
    #for i in range(len_time_correct):
    # if i >= POES_indexstart and i < POES_indexend:
    # print("%.f"%time_correct[i])

    print("Josh_E_correct_90")
    for i in range(len_time_correct):
        for j in range(3):
            if E_correct_0[i,j] <= 0:
                E_correct_0[i,j] = 0.01
            if E_correct_90[i,j] <= 0:
                E_correct_90[i,j] = 0.01
    print("finished with non-zeros")

    len_2s = len_time_correct
    #m_time_2sec = np.array([time_correct[0] + 2*1000*i for i in range(len_2s)])
    m_time_2sec = time_correct
    #pettit_time = [dt.datetime.fromtimestamp(int(m_time_2sec[i]/1000-86400), tz=timezone.utc) for i in range(len_2s)]  # pettit time is a day ahead  - IDL?
    pettit_time = [dt.datetime.fromtimestamp(int(m_time_2sec[i]/1000), tz=timezone.utc) for i in range(len_2s)]  # corrected 20220908 - no longer a day ahead!
    pstart = pettit_time[0]
    print("start")
    print(pstart)
    pend = pettit_time[len_2s-1]
    print("end")
    print(pend)
    flag1 = 0
    flag2 = 0
    poes_datetime = POEStime_orig
    for ntime in range(timeindex):
        poes_datetime[ntime] = poes_datetime[ntime].replace(tzinfo=timezone.utc)
    for ntime in range(timeindex):
        if flag1 == 0 and poes_datetime[ntime] > pstart:
            index_pstart = ntime
            flag1 = 1
        if flag2 == 0 and poes_datetime[ntime] > pend:
            index_pend = ntime
            flag2 = 1
    print("index_pstart")
    print(index_pstart)
    print("index_pend")
    print(index_pend)
    timespan = index_pend-index_pstart

    #pettit_interp0 = interp1d(time_correct,E_correct_0,kind='linear',axis=0,bounds_error=False,fill_value=0.)
    #pettit_interp90 = interp1d(time_correct,E_correct_90,kind='linear',axis=0,bounds_error=False,fill_value=0.)
    #E_correct_0_new[index_pstart:index_pend+1,:] = pettit_interp0(m_time_2sec)
    #E_correct_90_new[index_pstart:index_pend+1,:] = pettit_interp90(m_time_2sec)
    #E_correct_Gmean[index_pstart:index_pend+1,:] = np.sqrt(E_correct_0_new[index_pstart:index_pend+1,:] * E_correct_90_new[index_pstart:index_pend+1,:])

    print("calculating Gmean")
    E_correct_Gmean[index_pstart:index_pend+1,:] = np.sqrt(E_correct_0[index_pstart:index_pend+1,:] * E_correct_90[index_pstart:index_pend+1,:])
    print("finished calculating Gmean")

    BLC_flux_ln = np.zeros((len_time_correct,27),'f')
    BLCflux_int_ln = np.zeros((len_time_correct,5), 'f')
    BLCflux_int = np.zeros((len_time_correct,5), 'f')

    for ntime in range(len_time_correct):
        for nen in range(27):
            if BLC_flux[ntime,nen] > 0.:
                BLC_flux_ln[ntime,nen] = np.log(BLC_flux[ntime,nen])
            else:
                BLC_flux_ln[ntime,nen] = np.log(1.e-5)

    print("interpolating BLC flux to FIREBIRD energies")

    interpBLC_energy = interp1d(energy_27,BLC_flux_ln,kind='linear',axis=1,bounds_error=False,fill_value=-11.)
    BLCflux_int_ln[:,:] = interpBLC_energy(FB_en)
    for ntime in range(len_time_correct):
        BLCflux_int[ntime,:] = np.exp(BLCflux_int_ln[ntime,:])

    print(BLCflux_int[POES_indexstart:POES_indexend,3])
    #pettit_interpBLC = interp1d(time_correct,BLCflux_int,kind='linear',axis=0,bounds_error=False,fill_value=0.)
    #BLCflux_new[index_pstart:index_pend+1,:] = pettit_interpBLC(m_time_2sec)

    nc_fid_0.close()
    nc_fid_90.close()
    nc_fid_BLC.close()
    #return E_correct_0, E_correct_90, E_correct_Gmean, BLCflux_new
    return E_correct_0, E_correct_90, E_correct_Gmean, BLCflux_int

# calculate counts in each Yando energy bin
# for ntime in timestamp_list:

def calculate_counts(nlines, energy_FB, J0, E0, table4, table2):
    fluxall = np.zeros((nlines,5), dtype = float)
    FB_POES_calc_counts = np.zeros((nlines,4),dtype=float)
    for ntime in range(nlines):
        flux_en_lev = np.zeros((23), dtype = float)
        flux_en_bin = np.zeros((22), dtype = float)
        counts_en_bin = np.zeros((22,4), dtype = float)
        sum_E1 = 0.
        sum_E2 = 0.
        sum_E3 = 0.
        sum_P6 = 0.
        for nenergy in range(5):
            #fluxall[ntime, nenergy] =  J0[ntime]*np.exp(-energy_FB[nenergy]/1000./E0[ntime])
            fluxall[ntime, nenergy] =  J0[ntime]*np.exp(-energy_FB[nenergy]/E0[ntime])
        #timeprint = 121
        #if ntime == timeprint:
            #print("J0 = ",J0[timeprint])
            #print("E0 = ",E0[timeprint])
            #print("energy_levels = ", energy_FB)
            #print("fluxall = ", fluxall[timeprint,:])
        for en in range(22):
            #flux_en_lev[en] = J0[ntime]*np.exp(-table4[en,0]/1000./E0[ntime])
            flux_en_lev[en] = J0[ntime]*np.exp(-table4[en,0]/E0[ntime])
            #if ntime == timeprint and en == 21:
                #print("flux_en_lev = ", flux_en_lev)
        flux_en_lev[22] = J0[ntime]*np.exp(-5./E0[ntime])
        for en in range(22):
            #flux_en_bin[en] = ((flux_en_lev[en]+flux_en_lev[en+1])/2.)*(table4[en+1,0]-table4[en,0])
            flux_en_bin[en] =  integral_efolding(table4[en,0],table4[en+1,0], J0[ntime], E0[ntime])
            counts_en_bin[en,0] = flux_en_bin[en]*table4[en,2]
            counts_en_bin[en,1] = flux_en_bin[en]*table4[en,3]
            counts_en_bin[en,2] = flux_en_bin[en]*table4[en,4]
            counts_en_bin[en,3] = flux_en_bin[en]*table2[en,2]
            sum_E1 = sum_E1+counts_en_bin[en,0]
            sum_E2 = sum_E2+counts_en_bin[en,1]
            sum_E3 = sum_E3+counts_en_bin[en,2]
            #if ntime == timeprint:
                #print("energy = ", table4[en,0], " sum_E3 == ", sum_E3)
            sum_P6 = sum_P6+counts_en_bin[en,3]
        FB_POES_calc_counts[ntime,0] = sum_E1
        FB_POES_calc_counts[ntime,1] = sum_E2
        FB_POES_calc_counts[ntime,2] = sum_E3
        FB_POES_calc_counts[ntime,3] = sum_P6
    return FB_POES_calc_counts, fluxall

# Print values from FIREBIRD equivalent counts to output file
def save_POES_equiv_counts(filename, nlines, conj_index, Lshell, FB_POES_counts):
    output = open(filename, 'a')
    print("\nTitle: ", title_plot, file = output)
    print("\nFIREBIRD ", FBfig_txt_conj, file = output)
    print("POES ", POESfig_txt_conj, file = output)
    print("\nFIREBIRD time plot range", FBfig_text, file = output)
    print("POES time plot range", POESfig_text, file = output)
    print("FIREBIRD lat, lon range", FBfig_text_latlon, file = output)
    print("POES lat, lon range", POESfig_text_latlon, file = output)
    print("\nFIREBIRD energies", energy_legend, file = output)
    print("\nFIREBIRD EQUIVALENT COUNTS", file = output)
    flag3 = 0
    flag4 = 0
    flag5_minus = 0
    flag5 = 0
    flag5_plus = 0
    flag6 = 0
    flag7 = 0
    flag8 = 0
    flag9 = 0
    flag10 = 0
    L_conj = Lshell[conj_index]
    if Lshell[0] < Lshell[nlines-1]:
        for ntime in range(0,nlines):
            if flag3 == 0 and Lshell[ntime] >= 3:
                index_3L = ntime
                flag3 = 1
                print("   L = 3 values E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_3L,0], FB_POES_counts[index_3L,1], FB_POES_counts[index_3L,2], FB_POES_counts[index_3L,3]), file = output)
            if flag4 == 0 and Lshell[ntime] >= 4:
                index_4L = ntime
                flag4 = 1
                print("   L = 4 values E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_4L,0], FB_POES_counts[index_4L,1], FB_POES_counts[index_3L,2], FB_POES_counts[index_4L,3]), file = output)
            if flag5_minus == 0 and Lshell[ntime] >= 5.-0.25:
                index_5L_low = ntime
                flag5_minus = 1
            if flag5 == 0 and Lshell[ntime] >= 5.:
                index_5L = ntime
                flag5 = 1
                print("   L = 5 values E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_5L,0], FB_POES_counts[index_5L,1], FB_POES_counts[index_5L,2], FB_POES_counts[index_5L,3]), file = output)
            if flag5_plus == 0 and Lshell[ntime] >= 5.+0.25:
                index_5L_high = ntime
                flag5_plus = 1
                print("   ",file = output)
                print("    Average FB eq 0.5 L bin conjunction values at L = 5   E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (np.nanmean(FB_POES_counts[index_5L_low:index_5L_high,0]),np.nanmean(FB_POES_counts[index_5L_low:index_5L_high,1]), np.nanmean(FB_POES_counts[index_5L_low:index_5L_high,2]), np.nanmean(FB_POES_counts[index_5L_low:index_5L_high,3])), file = output)
            if flag6 == 0 and Lshell[ntime] >= 6:
                index_6L = ntime
                flag6 = 1
                print("   L = 6 values E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_6L,0], FB_POES_counts[index_6L,1], FB_POES_counts[index_6L,2], FB_POES_counts[index_6L,3]), file = output)
            if flag7 == 0 and Lshell[ntime] >= 6.98:
                index_7L = ntime
                flag7 = 1
                print("   L = 7 values E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_7L,0], FB_POES_counts[index_7L,1], FB_POES_counts[index_7L,2], FB_POES_counts[index_7L,3]), file = output)
            if flag8 == 0 and Lshell[ntime] <= L_conj-0.25:
                index_conjL_low = ntime
                flag8 = 1
            if flag9 == 0 and Lshell[ntime] <= L_conj:
                print("   Conjunction time values at L = ", L_conj ," E1, E2, E3, P6", file = output)
                index_conjL = ntime
                flag9 = 1
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_conjL,0],FB_POES_counts[index_conjL,1], FB_POES_counts[index_conjL,2], FB_POES_counts[index_conjL,3]), file = output)
            if flag10 == 0 and Lshell[ntime] <= L_conj+0.25:
                index_conjL_high = ntime
                flag8 = 1
    else:
        for ntime in range(0,nlines):
            if flag7 == 0 and Lshell[ntime] <= 7:
                index_7L = ntime
                flag7 = 1
                print("   L = 7 values E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_7L,0], FB_POES_counts[index_7L,1], FB_POES_counts[index_7L,2], FB_POES_counts[index_7L,3]), file = output)
            if flag6 == 0 and Lshell[ntime] <= 6:
                index_6L = ntime
                flag6 = 1
                print("   L = 6 values E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_6L,0], FB_POES_counts[index_6L,1], FB_POES_counts[index_6L,2], FB_POES_counts[index_6L,3]), file = output)
            if flag5_plus == 0 and Lshell[ntime] <= 5.+0.25:
                index_5L_high = ntime
                flag5_plus = 1
            if flag5 == 0 and Lshell[ntime] <= 5.:
                index_5L = ntime
                flag5 = 1
                print("   L = 5 values E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_5L,0], FB_POES_counts[index_5L,1], FB_POES_counts[index_5L,2], FB_POES_counts[index_5L,3]), file = output)
            if flag5_minus == 0 and Lshell[ntime] <= 5.-0.25:
                index_5L_low = ntime
                flag5_minus = 1
                print("   ",file = output)
                print("    Average FB eq 0.5 L bin conjunction values at L = 5   E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (np.nanmean(FB_POES_counts[index_5L_high:index_5L_low,0]),np.nanmean(FB_POES_counts[index_5L_high:index_5L_low,1]), np.nanmean(FB_POES_counts[index_5L_high:index_5L_low,2]), np.nanmean(FB_POES_counts[index_5L_high:index_5L_low,3])), file = output)
            if flag4 == 0 and Lshell[ntime] <= 4:
                index_4L = ntime
                flag4 = 1
                print("   L = 4 values E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_4L,0], FB_POES_counts[index_4L,1], FB_POES_counts[index_4L,2], FB_POES_counts[index_4L,3]), file = output)
            if flag3 == 0 and Lshell[ntime] <= 3.2:
                index_3L = ntime
                flag3 = 1
                print("   L = 7 values E1, E2, E3, P6", file = output)
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_3L,0], FB_POES_counts[index_3L,1], FB_POES_counts[index_3L,2], FB_POES_counts[index_3L,3]), file = output)
            if flag8 == 0 and Lshell[ntime] <= L_conj+0.25:
                index_conjL_low = ntime
                flag8 = 1
            if flag9 == 0 and Lshell[ntime] <= L_conj:
                print("   Conjunction time values at L = ", L_conj ," E1, E2, E3, P6", file = output)
                index_conjL = ntime
                flag9 = 1
                print("   %5.2f  %5.2f  %5.2f  %5.2f" % (FB_POES_counts[index_conjL,0],FB_POES_counts[index_conjL,1], FB_POES_counts[index_conjL,2], FB_POES_counts[index_conjL,3]), file = output)
            if flag10 == 0 and Lshell[ntime] <= L_conj-0.25:
                index_conjL_high = ntime
                flag8 = 1

    print("   ",file = output)
    print("    Average FB eq 0.5 L bin conjunction values at L = ", L_conj ," E1, E2, E3, P6", file = output)
    print("   %5.2f  %5.2f  %5.2f  %5.2f" % (np.nanmean(FB_POES_counts[index_conjL_low:index_conjL_high,0]),np.nanmean(FB_POES_counts[index_conjL_low:index_conjL_high,1]), np.nanmean(FB_POES_counts[index_conjL_low:index_conjL_high,2]), np.nanmean(FB_POES_counts[index_conjL_low:index_conjL_high,3])), file = output)
    print("   ",file = output)
    print("     max FB E1 = ", np.nanmax(FB_POES_counts[:,0]), file = output)
    print("     max FB E2 = ", np.nanmax(FB_POES_counts[:,1]), file = output)
    print("     max FB E3 = ", np.nanmax(FB_POES_counts[:,2]), file = output)
    print("     max FB P6 = ", np.nanmax(FB_POES_counts[:,3]), file = output)
    return


# identify time index for FB conjunction
def FB_conj_index(ntimestamps, timestamp, J0, E0, Lshell, MLT, conj, filename):
    output = open(filename, 'a')
    flag1 = 0
    ntime = 0
    index_start = 0
    index_end = 0
    FB_conjun_Lshell = 0.
    FB_conjun_MLT = 0.
    for ntime in range(ntimestamps):
        if flag1 == 0 and timestamp[ntime] >= conj:
            FB_index_conj = ntime
            flag1 = 1
    print("Conjunction time index for FB = ", FB_index_conj)
    print("At conjunction time =  ", timestamp[FB_index_conj], "J0 = %5.2f" %(J0[FB_index_conj]), "E0 = %5.2f" %(E0[FB_index_conj]), "L = %3.2f" %(Lshell[FB_index_conj]), file = output)
    #print("\nAt conjunction time = ", timestamp[FB_index_conj], "J0 = %5.2f" %J0[FB_index_conj], "E0 = %5.2f" %E0[FB_index_conj], file = output)
    print("   At Prior timestep = ", timestamp[FB_index_conj-1], "J0 = %5.2f" %J0[FB_index_conj-1], "E0 = %5.2f" %E0[FB_index_conj-1], "L = %3.2f" %(Lshell[FB_index_conj-1]), file = output)
    print("   At Next timestep = ", timestamp[FB_index_conj+1], "J0 = %5.2f" %J0[FB_index_conj+1], "E0 = %5.2f" %E0[FB_index_conj+1], "L = %3.2f" %(Lshell[FB_index_conj+1]), file = output)
    FB_conjun_Lshell = Lshell[FB_index_conj]
    FB_conjun_MLT = MLT[FB_index_conj]
    return FB_index_conj,FB_conjun_Lshell,FB_conjun_MLT


# identify time interval for POES plots
def time_interval(POESntime, start, end, conj, POEStime, L_POES):
    flag1 = 0
    flag2 = 0
    flag3 = 0
    ntime = 0
    index_start = 0
    index_end = 0
    poes_datetime = POEStime
    for ntime in range(POESntime):
        if flag1 == 0 and POEStime[ntime] > start:
            index_start = ntime
            flag1 = 1
        if flag2 == 0 and POEStime[ntime] > end:
            index_end = ntime
            flag2 = 1
        if flag3 == 0 and POEStime[ntime] >= conj:
            index_conj = ntime
            flag3 = 1
    print("POES indexstart = ", index_start)
    print("POES indexend = ", index_end)
    print("POES L of indexstart = ",L_POES[index_start])
    print("POES L of indexend = ",L_POES[index_end])
    print("Conjuntion time index for POES = ", index_conj)
    return index_start, index_end, index_conj

'''
# Identify time interval for corrected  counts
#def time_intervalcor(ntime_correct, start, end, flux, time_correct, Lshell_cor):
def time_intervalcor(ntime_correct, start, end, flux, time_correct):
    flag1 = 0
    flag2 = 0
    ntime = 0
    index_start = 0
    index_end = 0
    for ntime in range(ntime_correct):
        for channel in range(4):
            if flux[ntime,channel] < 1.:
                flux[ntime,channel] = 'nan'
            if flag1 == 0 and time_correct[ntime] > start:
                index_start = ntime
                flag1 = 1
            if flag2 == 0 and time_correct[ntime] > end:
                index_end = ntime
                flag2 = 2
    print("indexstart correct = ", index_start)
    print("indexend correct = ", index_end)
    #print(Lshell_cor[index_start])
    #print(Lshell_cor[index_end])
    return start, end, flux
'''

def save_POES_obs(filename, start, end, FB_index_conj, counts, correct_counts, L_POES, POES_tel): # Print POES values into output file
    output = open(filename, 'a')
    print("\n"+POES_tel, file = output)
    L_POES = abs(L_POES) #get rid of negative values of L
    flag3 = 0
    flag4 = 0
    flag5_minus = 0
    flag5 = 0
    flag5_plus = 0
    flag6 = 0
    flag7 = 0
    flag8 = 0
    flag9 = 0
    flag10 = 0
    L_conj = L_POES[FB_index_conj]
    if L_POES[start] < L_POES[end]:
        for index_L in range(start, end):
            for channel in range(4):
                if counts[index_L,channel] < 1.:
                    counts[index_L,channel] = 'nan'
            if flag3 == 0 and L_POES[index_L] >= 3:
                index_3L = index_L
                flag3 = 1
                print("Raw L = 3 values E1, E2, E3, P6", file = output)
                print(counts[index_3L,0],counts[index_3L,1], counts[index_3L,2], counts[index_3L,3], file = output)
                print("Correct L = 3 values  E2, E3, P6", file = output)
                print(correct_counts[index_3L,0],correct_counts[index_3L,1], correct_counts[index_3L,2],  file = output)
            if flag4 == 0 and L_POES[index_L] >= 4:
                index_4L = index_L
                flag4 = 1
                print("Raw L = 4 values E1, E2, E3, P6", file = output)
                print(counts[index_4L,0],counts[index_4L,1], counts[index_4L,2], counts[index_4L,3], file = output)
                print("Correct L = 4 values  E2, E3, P6", file = output)
                print(correct_counts[index_4L,0],correct_counts[index_4L,1], correct_counts[index_4L,2],file = output)
            if flag5_minus == 0 and L_POES[index_L] >= 4.75:
                index_5L_low = index_L
                flag5_minus = 1
                #print("L = 4.75 values E1, E2, E3, P6", file = output)
                #print(counts[index_5L_low,0],counts[index_5L_low,1], counts[index_5L_low,2], counts[index_5L_low,3], file = output)
            if flag5 == 0 and L_POES[index_L] >= 5.:
                index_5L = index_L
                flag5 = 1
                print("Raw L = 5 values E1, E2, E3, P6", file = output)
                print(counts[index_5L,0],counts[index_5L,1], counts[index_5L,2], counts[index_5L,3], file = output)
                print("Correct L = 5 values  E2, E3, P6", file = output)
                print(correct_counts[index_5L,0],correct_counts[index_5L,1], correct_counts[index_5L,2],  file = output)
            if flag5_plus == 0 and L_POES[index_L] >= 5.25 :
                index_5L_high = index_L
                #print("L = 5.25 values E1, E2, E3, P6", file = output)
                #print("index_5L_low", index_5L_low)
                #print("index_5L_high", index_5L_high)
                #print(counts[index_5L_high,0],counts[index_5L_high,1], counts[index_5L_high,2], counts[index_5L_high,3], file = output)
                print("   ",file = output)
                print("Raw Average POES 0.5 L bin conjunction values at L = 5   E1, E2, E3, P6 ",file = output)
                print(np.nanmean(counts[index_5L_low:index_5L_high,0]),np.nanmean(counts[index_5L_low:index_5L_high,1]), np.nanmean(counts[index_5L_low:index_5L_high,2]), np.nanmean(counts[index_5L_low:index_5L_high,3]), file = output)
                print("   ",file = output)
                print("   ",file = output)
                print("Correct Average POES 0.5 L bin conjunction values at L = 5    E2, E3, P6 ",file = output)
                print(np.nanmean(correct_counts[index_5L_low:index_5L_high,0]),np.nanmean(correct_counts[index_5L_low:index_5L_high,1]), np.nanmean(correct_counts[index_5L_low:index_5L_high,2]),  file = output)
                print("   ",file = output)
                flag5_plus = 1
            if flag6 == 0 and L_POES[index_L] >= 6:
                index_6L = index_L
                flag6 = 1
                print("Raw L = 6 values E1, E2, E3, P6", file = output)
                print(counts[index_6L,0],counts[index_6L,1], counts[index_6L,2], counts[index_6L,3], file = output)
                print("Correct L = 6 values  E2, E3, P6", file = output)
                print(correct_counts[index_6L,0],correct_counts[index_6L,1], correct_counts[index_6L,2],  file = output)
            if flag7 == 0 and L_POES[index_L] >= 7:
                index_7L = index_L
                flag7 = 1
                print("Raw L = 7 values E1, E2, E3, P6", file = output)
                print(counts[index_7L,0],counts[index_7L,1], counts[index_7L,2], counts[index_7L,3], file = output)
                print("Correct L = 7 values  E2, E3, P6", file = output)
                print(correct_counts[index_7L,0],correct_counts[index_7L,1], correct_counts[index_7L,2],  file = output)
            if flag8 == 0 and L_POES[index_L] >= L_conj-0.25:
                index_conjL_low = index_L
                flag8 = 1
            if flag9 == 0 and L_POES[index_L] >= L_conj:
                index_conjL = index_L
                flag9 = 1
                print("Raw Conjunction time values at L = ", L_conj ," E1, E2, E3, P6", file = output)
                print(counts[index_conjL,0],counts[index_conjL,1], counts[index_conjL,2], counts[index_conjL,3], file = output)
                print("Correct Conjunction time values at L = ", L_conj ,"  E2, E3, P6", file = output)
                print(correct_counts[index_conjL,0],correct_counts[index_conjL,1], correct_counts[index_conjL,2], file = output)
            if flag10 == 0 and L_POES[index_L] >= L_conj+0.25:
                index_conjL_high = index_L
                flag10 = 1
                print("   ",file = output)
                print("Raw Average POES 0.5 L bin conjunction values at L = ", L_conj ," E1, E2, E3, P6", file = output)
                print(np.nanmean(counts[index_conjL_low:index_conjL_high,0]),np.nanmean(counts[index_conjL_low:index_conjL_high,1]), np.nanmean(counts[index_conjL_low:index_conjL_high,2]), np.nanmean(counts[index_conjL_low:index_conjL_high,3]), file = output)
                print("   ",file = output)
                print("Correct Average POES 0.5 L bin conjunction values at L = ", L_conj ," E2, E3, P6", file = output)
                print("   ",file = output)
                print(np.nanmean(correct_counts[index_conjL_low:index_conjL_high,0]),np.nanmean(correct_counts[index_conjL_low:index_conjL_high,1]), np.nanmean(correct_counts[index_conjL_low:index_conjL_high,2]),  file = output)
                print("   ",file = output)
        #print("index_3L = ", index_3L, file = output)
        #print("index_7L = ", index_7L, file = output)

        print("raw max POES E1 = ", np.nanmax(counts[index_3L:index_7L,0]), file = output)
        print("raw max POES E2 = ", np.nanmax(counts[index_3L:index_7L,1]), file = output)
        print("raw max POES E3 = ", np.nanmax(counts[index_3L:index_7L,2]), file = output)
        print("raw max POES P6 = ", np.nanmax(counts[index_3L:index_7L,3]), file = output)
        print("   ",file = output)
        print("correct max POES E2 = ", np.nanmax(correct_counts[index_3L:index_7L,0]), file = output)
        print("correct max POES E3 = ", np.nanmax(correct_counts[index_3L:index_7L,1]), file = output)
        print("correct max POES P6 = ", np.nanmax(correct_counts[index_3L:index_7L,2]), file = output)
    else:
        for index_L in range(end, start, -1):
            for channel in range(4):
                if counts[index_L,channel] < 1.:
                    counts[index_L,channel] = 'nan'
                if flag3 == 0 and L_POES[index_L] >= 3:
                    index_3L = index_L
                    flag3 = 1
                    print("Raw L = 3 values E1, E2, E3, P6", file = output)
                    print(counts[index_3L,0],counts[index_3L,1], counts[index_3L,2], counts[index_3L,3], file = output)
                    print("Correct L = 3 values E2, E3, P6", file = output)
                    print(correct_counts[index_3L,0],correct_counts[index_3L,1], correct_counts[index_3L,2], file = output)
                if flag4 == 0 and L_POES[index_L] >= 4:
                    index_4L = index_L
                    flag4 = 1
                    print("Raw L = 4 values E1, E2, E3, P6", file = output)
                    print(counts[index_4L,0],counts[index_4L,1], counts[index_4L,2], counts[index_4L,3], file = output)
                    print("Correct L = 4 values  E2, E3, P6", file = output)
                    print(correct_counts[index_4L,0],correct_counts[index_4L,1], correct_counts[index_4L,2], file = output)
                if flag5_minus == 0 and L_POES[index_L] >= 4.75:
                    index_5L_low = index_L
                    flag5_minus = 1
                    #print("L = 4.75 values E1, E2, E3, P6", file = output)
                    #print(counts[index_5L_low,0],counts[index_5L_low,1], counts[index_5L_low,2], counts[index_5L_low,3], file = output)
                if flag5 == 0 and L_POES[index_L] >= 5.:
                    index_5L = index_L
                    flag5 = 1
                    print("Raw L = 5 values E1, E2, E3, P6", file = output)
                    print(counts[index_5L,0],counts[index_5L,1], counts[index_5L,2], counts[index_5L,3], file = output)
                    print("Correct L = 5 values E2, E3, P6", file = output)
                    print(correct_counts[index_5L,0],correct_counts[index_5L,1], correct_counts[index_5L,2], file = output)
                if flag5_plus == 0 and L_POES[index_L] >= 5.25:
                    index_5L_high = index_L
                    #print("L = 5.25 values E1, E2, E3, P6", file = output)
                    #print(counts[index_5L_low,0],counts[index_5L_low,1], counts[index_5L_low,2], counts[index_5L_low,3], file = output)
                    print("   ",file = output)
                    print("Raw Average POES 0.5 L bin conjunction values at L = 5   E1, E2, E3, P6 ",file = output)
                    print(np.nanmean(counts[index_5L_high:index_5L_low,0]),np.nanmean(counts[index_5L_high:index_5L_low,1]), np.nanmean(counts[index_5L_high:index_5L_low,2]), np.nanmean(counts[index_5L_high:index_5L_low,3]), file = output)
                    print("   ",file = output)
                    print("   ",file = output)
                    print("Correct Average POES 0.5 L bin conjunction values at L = 5   E2, E3, P6 ",file = output)
                    print(np.nanmean(correct_counts[index_5L_high:index_5L_low,0]),np.nanmean(correct_counts[index_5L_high:index_5L_low,1]), np.nanmean(correct_counts[index_5L_high:index_5L_low,2]),  file = output)
                    print("   ",file = output)
                    flag5_plus = 1
                if flag6 == 0 and L_POES[index_L] >= 6:
                    index_6L = index_L
                    flag6 = 1
                    print("Raw L = 6 values E1, E2, E3, P6", file = output)
                    print(counts[index_6L,0],counts[index_6L,1], counts[index_6L,2], counts[index_6L,3], file = output)
                    print("Correct L = 6 values E2, E3, P6", file = output)
                    print(correct_counts[index_6L,0],correct_counts[index_6L,1], correct_counts[index_6L,2], file = output)
                if flag7 == 0 and L_POES[index_L] >= 7:
                    index_7L = index_L
                    flag7 = 1
                    print("Raw L = 7 values E1, E2, E3, P6", file = output)
                    print(counts[index_7L,0],counts[index_7L,1], counts[index_7L,2], counts[index_7L,3], file = output)
                    print("Correct L = 7 values E2, E3, P6", file = output)
                    print(correct_counts[index_7L,0],correct_counts[index_7L,1], correct_counts[index_7L,2],  file = output)
                if flag8 == 0 and L_POES[index_L] >= L_conj-0.25:
                    index_conjL_low = index_L
                    flag8 = 1
                if flag9 == 0 and L_POES[index_L] >= L_conj:
                    index_conjL = index_L
                    flag9 = 1
                    print("Raw Conjunction time values at L = ", L_conj ," E1, E2, E3, P6", file = output)
                    print(counts[index_conjL,0],counts[index_conjL,1], counts[index_conjL,2], counts[index_conjL,3], file = output)
                    print("Correct Conjunction time values at L = ", L_conj ,"  E2, E3, P6", file = output)
                    print(correct_counts[index_conjL,0],correct_counts[index_conjL,1], correct_counts[index_conjL,2],  file = output)
                if flag10 == 0 and L_POES[index_L] >= L_conj+0.25:
                    index_conjL_high = index_L
                    flag10 = 1
                    print("   ",file = output)
                    print("Raw Average POES 0.5 L bin conjunction values at L = ", L_conj ," E1, E2, E3, P6", file = output)
                    print(np.nanmean(counts[index_conjL_high:index_conjL_low,0]),np.nanmean(counts[index_conjL_high:index_conjL_low,1]), np.nanmean(counts[index_conjL_high:index_conjL_low,2]), np.nanmean(counts[index_conjL_high:index_conjL_low,3]), file = output)
                    print("   ",file = output)
                    print("Correct Average POES 0.5 L bin conjunction values at L = ", L_conj ,"  E2, E3, P6", file = output)
                    print("   ",file = output)
                    print(np.nanmean(correct_counts[index_conjL_high:index_conjL_low,0]),np.nanmean(correct_counts[index_conjL_high:index_conjL_low,1]), np.nanmean(correct_counts[index_conjL_high:index_conjL_low,2]),  file = output)
                    print("   ",file = output)
        #print("index_3L = ", index_3L)
        #print("index_7L = ", index_7L)
        print("raw max POES E1 = ", np.nanmax(counts[index_7L:index_3L,0]), file = output)
        print("raw max POES E2 = ", np.nanmax(counts[index_7L:index_3L,1]), file = output)
        print("raw max POES E3 = ", np.nanmax(counts[index_7L:index_3L,2]), file = output)
        print("raw max POES P6 = ", np.nanmax(counts[index_7L:index_3L,3]), file = output)
        print("   ",file = output)
        print("correct max POES E2 = ", np.nanmax(correct_counts[index_7L:index_3L,0]), file = output)
        print("correct max POES E3 = ", np.nanmax(correct_counts[index_7L:index_3L,1]), file = output)
        print("correct max POES P6 = ", np.nanmax(correct_counts[index_7L:index_3L,2]), file = output)
    return

def save_flux(filename, start, end, FB_index_conj, POES_BLCflux, L_POES, title, energy_legend, POEStime):
    output = open(filename, 'a')
    lenPOEStime = len(POEStime)
    print("\n"+title, file = output)
    L_POES = abs(L_POES) #get rid of negative values of L
    flag3 = 0
    flag4 = 0
    flag5_minus = 0
    flag5 = 0
    flag5_plus = 0
    flag6 = 0
    flag7 = 0
    flag8 = 0
    flag9 = 0
    flag10 = 0
    L_conj = L_POES[FB_index_conj]
    if L_POES[start] < L_POES[end]:
        for index_L in range(start, end):
            for channel in range(5):
                if POES_BLCflux[index_L,channel] < 0.001:
                    POES_BLCflux[index_L,channel] = 'nan'
            if flag3 == 0 and L_POES[index_L] >= 3:
                index_3L = index_L
                flag3 = 1
                print("\n BLC flux L = 3 for energies",energy_legend, file = output)
                print(POES_BLCflux[index_3L,0],POES_BLCflux[index_3L,1], POES_BLCflux[index_3L,2], POES_BLCflux[index_3L,3], POES_BLCflux[index_3L,4],file = output)
            if flag4 == 0 and L_POES[index_L] >= 4:
                index_4L = index_L
                flag4 = 1
                print("\n BLC flux L = 4 for energies ",energy_legend, file = output)
                print(POES_BLCflux[index_4L,0],POES_BLCflux[index_4L,1], POES_BLCflux[index_4L,2], POES_BLCflux[index_4L,3], POES_BLCflux[index_4L,4],file = output)
            if flag5_minus == 0 and L_POES[index_L] >= 4.75:
                index_5L_low = index_L
                flag5_minus = 1
                #print("L = 4.75 values E1, E2, E3, P6", file = output)
                #print(counts[index_5L_low,0],counts[index_5L_low,1], counts[index_5L_low,2], counts[index_5L_low,3], file = output)
            if flag5 == 0 and L_POES[index_L] >= 5.:
                index_5L = index_L
                flag5 = 1
                print("\n BLC flux L = 5 for energies ",energy_legend, file = output)
                print(POES_BLCflux[index_5L,0],POES_BLCflux[index_5L,1], POES_BLCflux[index_5L,2], POES_BLCflux[index_5L,3], POES_BLCflux[index_5L,4],file = output)
            if flag5_plus == 0 and L_POES[index_L] >= 5.25 :
                index_5L_high = index_L
                #print("L = 5.25 values E1, E2, E3, P6", file = output)
                #print("index_5L_low", index_5L_low)
                #print("index_5L_high", index_5L_high)
                #print(counts[index_5L_high,0],counts[index_5L_high,1], counts[index_5L_high,2], counts[index_5L_high,3], file = output)
                print("\n BLC flux Average 0.5 L bin around L = 5 for energies  ",energy_legend, file = output)
                print(np.nanmean(POES_BLCflux[index_5L_low:index_5L_high,0]),np.nanmean(POES_BLCflux[index_5L_low:index_5L_high,1]), np.nanmean(POES_BLCflux[index_5L_low:index_5L_high,2]), np.nanmean(POES_BLCflux[index_5L_low:index_5L_high,3]),np.nanmean(POES_BLCflux[index_5L_low:index_5L_high,4]), file = output)
                print("   ",file = output)
                flag5_plus = 1
            if flag6 == 0 and L_POES[index_L] >= 6:
                index_6L = index_L
                flag6 = 1
                print("\n BLC flux L = 6 for energies  ",energy_legend,  file = output)
                print(POES_BLCflux[index_6L,0],POES_BLCflux[index_6L,1], POES_BLCflux[index_6L,2], POES_BLCflux[index_6L,3], POES_BLCflux[index_6L,4],file = output)
            if flag7 == 0 and L_POES[index_L] >= 7:
                index_7L = index_L
                flag7 = 1
                print("\n BLC flux L = 7 for energies  ",energy_legend,  file = output)
                print(POES_BLCflux[index_7L,0],POES_BLCflux[index_7L,1], POES_BLCflux[index_7L,2], POES_BLCflux[index_7L,3], POES_BLCflux[index_7L,4],file = output)
            if flag8 == 0 and L_POES[index_L] >= L_conj-0.25:
                index_conjL_low = index_L
                flag8 = 1
            if flag9 == 0 and L_POES[index_L] >= L_conj:
                index_conjL = index_L
                flag9 = 1
                print("\n BLC flux Conjunction time values at L = ", L_conj ," for energies  ",energy_legend,  file = output)
                print(POES_BLCflux[index_conjL,0],POES_BLCflux[index_conjL,1], POES_BLCflux[index_conjL,2], POES_BLCflux[index_conjL,3], POES_BLCflux[index_conjL,4],file = output)
                print("   ",file = output)
            if flag10 == 0 and L_POES[index_L] >= L_conj+0.25:
                index_conjL_high = index_L
                flag10 = 1
                print("\n BLC flux Average 0.5 L bin around conjunction values at L = ", L_conj ," for energies  ",energy_legend,  file = output)
                print(np.nanmean(POES_BLCflux[index_conjL_low:index_conjL_high,0]),np.nanmean(POES_BLCflux[index_conjL_low:index_conjL_high,1]), np.nanmean(POES_BLCflux[index_conjL_low:index_conjL_high,2]), np.nanmean(POES_BLCflux[index_conjL_low:index_conjL_high,3]), np.nanmean(POES_BLCflux[index_conjL_low:index_conjL_high,4]),file = output)
                print("   ",file = output)
        #print("index_3L = ", index_3L, file = output)
        #print("index_7L = ", index_7L, file = output)
        print("\n BLC flux max POES for energies  ",energy_legend, "\n", np.nanmax(POES_BLCflux[start:end,0]),np.nanmax(POES_BLCflux[start:end,1]), np.nanmax(POES_BLCflux[start:end,2]), np.nanmax(POES_BLCflux[start:end,3]), np.nanmax(POES_BLCflux[start:end,4]), file = output)
        print("   ",file = output)
    else:
        for index_L in range(end, start, -1):
            for channel in range(5):
                if POES_BLCflux[index_L,channel] < 0.001:
                    POES_BLCflux[index_L,channel] = 'nan'
            if flag3 == 0 and L_POES[index_L] >= 3:
                index_3L = index_L
                flag3 = 1
                print("\n BLC flux L = 3 values for energies ",energy_legend,  file = output)
                print(POES_BLCflux[index_3L,0],POES_BLCflux[index_3L,1], POES_BLCflux[index_3L,2], POES_BLCflux[index_3L,3], POES_BLCflux[index_3L,4], file = output)
            if flag4 == 0 and L_POES[index_L] >= 4:
                index_4L = index_L
                flag4 = 1
                print("\n BLC flux L = 4 values for energies ",energy_legend,  file = output)
                print(POES_BLCflux[index_4L,0],POES_BLCflux[index_4L,1], POES_BLCflux[index_4L,2], POES_BLCflux[index_4L,3], POES_BLCflux[index_4L,4],file = output)
            if flag5_minus == 0 and L_POES[index_L] >= 4.75:
                index_5L_low = index_L
                flag5_minus = 1
                    #print("L = 4.75 values E1, E2, E3, P6", file = output)
                    #print(POES_BLCflux[index_5L_low,0],POES_BLCflux[index_5L_low,1], POES_BLCflux[index_5L_low,2], POES_BLCflux[index_5L_low,3], POES_BLCflux[index_5L_low,4],file = output)
            if flag5 == 0 and L_POES[index_L] >= 5.:
                index_5L = index_L
                flag5 = 1
                print("\n BLC flux L = 5 values for energies ",energy_legend,  file = output)
                print(POES_BLCflux[index_5L,0],POES_BLCflux[index_5L,1], POES_BLCflux[index_5L,2], POES_BLCflux[index_5L,3], POES_BLCflux[index_5L,4],file = output)
            if flag5_plus == 0 and L_POES[index_L] >= 5.25:
                index_5L_high = index_L
                #print("L = 5.25 values E1, E2, E3, P6", file = output)
                #print(POES_BLCflux[index_5L_low,0],POES_BLCflux[index_5L_low,1], POES_BLCflux[index_5L_low,2], POES_BLCflux[index_5L_low,3], fPOES_BLCflux[index_5L_low,4], file = output
                print("\n BLC flux Average 0.5 L bin around L = 5 for energies ",energy_legend, file = output)
                print(np.nanmean(POES_BLCflux[index_5L_high:index_5L_low,0]),np.nanmean(POES_BLCflux[index_5L_high:index_5L_low,1]), np.nanmean(POES_BLCflux[index_5L_high:index_5L_low,2]), np.nanmean(POES_BLCflux[index_5L_high:index_5L_low,3]), np.nanmean(POES_BLCflux[index_5L_high:index_5L_low,4]),file = output)
                print("   ",file = output)
                flag5_plus = 1
            if flag6 == 0 and L_POES[index_L] >= 6:
                index_6L = index_L
                flag6 = 1
                print("\n BLC flux L = 6 values for energies ",energy_legend, file = output)
                print(POES_BLCflux[index_6L,0],POES_BLCflux[index_6L,1], POES_BLCflux[index_6L,2], POES_BLCflux[index_6L,3], POES_BLCflux[index_6L,4], file = output)
            if flag7 == 0 and L_POES[index_L] >= 7:
                index_7L = index_L
                flag7 = 1
                print("\n BLC flux L = 7 values for energies ",energy_legend,  file = output)
                print(POES_BLCflux[index_7L,0],POES_BLCflux[index_7L,1], POES_BLCflux[index_7L,2], POES_BLCflux[index_7L,3], POES_BLCflux[index_7L,4],file = output)
            if flag8 == 0 and L_POES[index_L] >= L_conj-0.25:
                index_conjL_low = index_L
                flag8 = 1
            if flag9 == 0 and L_POES[index_L] >= L_conj:
                index_conjL = index_L
                flag9 = 1
                print("\n BLC flux Conjunction time values at L = ", L_conj ," for energies ",energy_legend,  file = output)
                print(POES_BLCflux[index_conjL,0],POES_BLCflux[index_conjL,1], POES_BLCflux[index_conjL,2], POES_BLCflux[index_conjL,3], POES_BLCflux[index_conjL,4], file = output)
                print("   ",file = output)
            if flag10 == 0 and L_POES[index_L] >= L_conj+0.25:
                index_conjL_high = index_L
                flag10 = 1
                print("\n BLC flux Average 0.5 L bin around conjunction values at L = ", L_conj ," for energies ",energy_legend, file = output)
                print(np.nanmean(POES_BLCflux[index_conjL_high:index_conjL_low,0]),np.nanmean(POES_BLCflux[index_conjL_high:index_conjL_low,1]), np.nanmean(POES_BLCflux[index_conjL_high:index_conjL_low,2]), np.nanmean(POES_BLCflux[index_conjL_high:index_conjL_low,3]), np.nanmean(POES_BLCflux[index_conjL_high:index_conjL_low,4]),file = output)
                print("   ",file = output)
                    #print("index_3L = ", index_3L)
                    #print("index_7L = ", index_7L)
        print("\n BLC flux max POES for energies ",energy_legend, "\n", np.nanmax(POES_BLCflux[start:end,0]),  np.nanmax(POES_BLCflux[start:end,1]), np.nanmax(POES_BLCflux[start:end,2]), np.nanmax(POES_BLCflux[start:end,3]), np.nanmax(POES_BLCflux[start:end,4]),file = output)
        print("   ",file = output)
    return

def read_POES_Lshell(filename):
    dataList = []
    try:
        file = open(filename, 'r')
    except:
        print('The file may not exist or the program may have not have been able to open it.')
        return([])
    else:
        file = open(filename, 'r')
        for line in file:
            line = line.strip().split(',')
            dataList.append(line)
    file.close()
    datetime_list = []
    Lshell_list = []
    MLT_list = []
    for row in dataList:
        time = dt.datetime.strptime(row[0], '%Y-%m-%d %H:%M:%S.%f')
        Lshell = np.abs(float(row[1]))
        MLT = float(row[2])
        datetime_list.append(time)
        Lshell_list.append(Lshell)
        MLT_list.append(MLT)
    Lshell_array = np.array(Lshell_list)
    MLT_array = np.array(MLT_list)
    return Lshell_array, MLT_array

#create_plot_Lshell(Lshell, flux_all, energy_legend,FBfig_text, FBfig_text_latlon, FB_POES_counts, indexstart, indexend, L_T89_POES, POES_obs_counts,POESsc,POESfig_text, POESfig_text_latlon, POES_90_obs_counts, filename_plot_L, FB_conj_L, FB_conj_MLT)

def create_plot_Lshell(Lshell, FB_flux_all, energy_legend,FBfig_text, FBfig_text_latlon, FB_POES_counts, indexstart, indexend, L_T89_POES, POES_obs_counts,POESsc,POESfig_text, POESfig_text_latlon, POES_90_obs_counts, POES_Gmean_obs_counts, filename_plot_L, Lshell_conjunction, MLT_conjunction, Correct_0_counts, Correct_90_counts, Correct_mean_counts, BLCflux, timeindex):
    
    POES_0_all = np.zeros((timeindex,7),'f')
    POES_90_all = np.zeros((timeindex,7),'f')
    POES_mean_all = np.zeros((timeindex,7),'f')
    POES_0_all[:,0:4] = POES_obs_counts[:,:]
    POES_0_all[:,4:8] = Correct_0_counts[:,:]
    POES_90_all[:,0:4] = POES_90_obs_counts[:,:]
    POES_90_all[:,4:8] = Correct_90_counts[:,:]
    POES_mean_all[:,0:4] = POES_Gmean_obs_counts[:,:]
    POES_mean_all[:,4:8] = Correct_mean_counts[:,:]
    
    Ltitle = str(Lshell_conjunction)
    MLTtitle = str(MLT_conjunction)
    # Make panel plots
    #fig, axs = plt.subplots(5)
    fig, axs = plt.subplots(7)
    #fig.suptitle(title_plot+' L = '+Ltitle+' FB-MLT = '+MLTtitle)
    fig.suptitle(FBfig_text+'\n '+FBfig_text_latlon+'\n L = '+Ltitle+', MLT = '+MLTtitle)
    fig.set_size_inches(horz_size, vert_size_L)
    fig.frameon=False
    colorsFB = ['blue','red','green']
    colorsPOES = ['#D3D3D3','lightblue','pink','lightgreen','blue', 'red', 'green']
    colorsflux = ['orange','red','purple','green','cyan']
    axs[0].set_prop_cycle(color=colorsflux)
    axs[1].set_prop_cycle(color=colorsFB)
    axs[2].set_prop_cycle(color=colorsPOES)
    axs[3].set_prop_cycle(color=colorsPOES)
    axs[4].set_prop_cycle(color=colorsPOES)
    axs[5].set_prop_cycle(color=colorsflux)
    axs[6].set_prop_cycle(color=colorsflux)

    # Plot FIREBIRD flux
    axs[0].plot(Lshell,FB_flux_all, marker='.', markersize='3', linestyle='none')
    axs[0].set_yscale('log')
    axs[0].set_ylim([0.01,1.e2])
    axs[0].set_yticks([1.e-2,1.e-1, 1.e0, 1.e1, 1.e2])
    #locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=5)
    #axs[0].yaxis.set_minor_locator(locmin)
    #axs[0].set_xlim([3,7])
    axs[0].set_xlim([3,8])
    axs[0].set_xlabel('L shell')
    axs[0].set_ylabel('FIREBIRD-II \n '+FBsc+' Flux')
    axs[0].legend((energy_legend), loc='upper right', fontsize = 'small')
    #axs[0].text(3.1, 30, FBfig_text, {'color': 'k', 'fontsize': 8})
    #axs[0].text(3.1, 10, FBfig_text_latlon, {'color': 'k', 'fontsize': 8})
    #axs[0].text(3.1, 3, FBfig_txt_conj,{'color': 'k', 'fontsize': 8})
    axs[0].axvline(x=Lshell_conjunction, color='k')
    axs[0].grid(True)
    
    # Plot FIREBIRD equivalent counts
    axs[1].plot(Lshell,FB_POES_counts[:,1:4], marker='.', markersize='3', linestyle='none')
    axs[1].set_yscale('log')
    axs[1].set_ylim([1.,1.e4])
    axs[1].set_yticks([1., 1.e1, 1.e2, 1.e3, 1.e4])
    #axs[1].set_xlim([3,7])
    axs[1].set_xlim([3,8])
    axs[1].set_xlabel('L shell')
    axs[1].set_ylabel('Predicted \n Counts \n  ')
    axs[1].legend(('> 100 keV', '> 300 keV', '> 600 keV'), loc='upper right', fontsize = 'small')
    axs[1].axvline(x=Lshell_conjunction, color='k')
    axs[1].grid(True)
    
    # Plot POES counts 0tel
    #axs[2].plot(L_T89_POES[indexstart:indexend], POES_obs_counts[indexstart:indexend,:])
    axs[2].plot(L_T89_POES[indexstart:indexend], POES_0_all[indexstart:indexend,:])
    
    #axs[2].set_xlim([3,7])
    axs[2].set_xlim([3,8])
    axs[2].set_yscale('log')
    axs[2].set_ylim([1.,1.e4])
    axs[2].set_yticks([1., 1.e1, 1.e2, 1.e3, 1.e4])
    axs[2].set_xlabel('L shell')
    axs[2].set_ylabel(POESsc+'\n tel0 \n Counts \n ')
    axs[2].legend(('> 40 keV', '> 100 keV', '> 300 keV', '> 600 keV'),loc='upper right', fontsize = 'small')
    #axs[2].text(3.1, 2000, POESfig_text, {'color': 'k', 'fontsize': 8})
    #axs[2].text(3.1, 400, POESfig_text_latlon, {'color': 'k', 'fontsize': 8})
    #axs[2].text(3.1, 100, POESfig_txt_conj,{'color': 'k', 'fontsize': 8})
    axs[2].axvline(x=Lshell_conjunction, color='k')
    axs[2].grid(True)
    
    # Plot POES counts Geometric mean
    axs[3].plot(L_T89_POES[indexstart:indexend], POES_mean_all[indexstart:indexend,:])
    #axs[3].set_xlim([3,7])
    axs[3].set_xlim([3,8])
    axs[3].set_yscale('log')
    axs[3].set_ylim([1.,1.e4])
    axs[3].set_yticks([1., 1.e1, 1.e2, 1.e3, 1.e4])
    axs[3].set_xlabel('L shell')
    axs[3].set_ylabel(POESsc+'\n Geometric \n Mean')
    axs[3].legend(('> 40 keV', '> 100 keV', '> 300 keV', '> 600 keV'), loc='upper right', fontsize = 'small')
    axs[3].axvline(x=Lshell_conjunction, color='k')
    axs[3].grid(True)
    
    # Plot POES counts 90tel
    axs[4].plot(L_T89_POES[indexstart:indexend], POES_90_all[indexstart:indexend,:])
    #axs[4].set_xlim([3,7])
    axs[4].set_xlim([3,8])
    axs[4].set_yscale('log')
    axs[4].set_ylim([1.,1.e4])
    axs[4].set_yticks([1., 1.e1, 1.e2, 1.e3, 1.e4])
    axs[4].set_xlabel('L shell')
    axs[4].set_ylabel(POESsc+'\n tel90 \n Counts \n  ')
    axs[4].legend(('> 40 keV', '> 100 keV', '> 300 keV', '> 600 keV'), loc='upper right', fontsize = 'small')
    axs[4].axvline(x=Lshell_conjunction, color='k')
    axs[4].grid(True)
    
    '''
        # Plot Corrected POES counts
        #axs[4].plot(POES_time[indexstartcor:indexendcor], E_correct[indexstartcor:indexendcor,:], marker=".")
        axs[4].plot(L_T89_POES[indexstartcor:indexendcor], E_correct[indexstartcor:indexendcor,:], marker=".")
        #axs[4].set_xlim([dt.datetime(2018,9,22,19,43,00), dt.datetime(2018,9,22,19,55,30)])axs[4].set_xlim([3,7])
        axs[4].set_yscale('log')
        axs[4].set_ylim([1.,1.e4])
        axs[4].set_yticks([1., 1.e1, 1.e2, 1.e3])
        #axs[4].set_xlabel('time')
        axs[4].set_xlabel('L shell')
        axs[4].set_ylabel('corrected')
        axs[4].legend(('> 40 keV', '> 100 keV', '> 300 keV', '> 600 keV'),
        loc='upper right', fontsize = 'small')
        '''
    
    
    axs[5].plot(L_T89_POES[indexstart:indexend], BLCflux[indexstart:indexend,:])
    #axs[5].set_xlim([3,7])
    axs[5].set_xlim([3,8])
    axs[5].set_yscale('log')
    axs[5].set_ylim([0.01,1.e3])
    axs[5].set_yticks([1.e-2,1.e-1, 1.e0, 1.e1, 1.e2, 1.e3])
    axs[5].set_xlabel('L shell')
    axs[5].set_ylabel('POES BLC_flux')
    axs[5].legend((energy_legend), loc='upper right', fontsize = 'small')
    axs[5].axvline(x=Lshell_conjunction, color='k')
    axs[5].grid(True)
    
    #FB_BLCflux = FB_flux_all
    #FB_BLCflux[:,:] = FB_flux_all[:,:]*3.758
    FB_BLCflux = FB_flux_all*3.758
    
    # Plot FIREBIRD BLC flux
    axs[6].plot(Lshell,FB_BLCflux, marker='.', markersize='3', linestyle='none')
    #axs[6].set_xlim([3,7])
    axs[6].set_xlim([3,8])
    axs[6].set_yscale('log')
    axs[6].set_ylim([0.01,1.e3])
    axs[6].set_yticks([1.e-2,1.e-1, 1.e0, 1.e1, 1.e2, 1.e3])
    axs[6].set_xlabel('L shell')
    axs[6].set_ylabel('FIREBIRD-II \n '+FBsc+'BLC Flux')
    axs[6].legend((energy_legend), loc='upper right', fontsize = 'small')
    axs[6].axvline(x=Lshell_conjunction, color='k')
    axs[6].grid(True)
    
    fig.autofmt_xdate()
    #plt.savefig(filename_plot_L, orientation='landscape', format='ps')   #MLT
    #plt.savefig(filename_plot_L, orientation='portrait', format = 'ps')   #MLT
    plt.savefig(filename_plot_L, format = 'ps', bbox_inches='tight')  
    #plt.show()

def create_plot_Lshell_all(Lshell, flux_all, energy_legend,FBfig_text, FBfig_text_latlon, FB_POES_counts, indexstart, indexend, L_T89_POES, POES_obs_counts,POESsc,POESfig_text, POESfig_text_latlon, POES_90_obs_counts, POES_Gmean_obs_counts, filename_plot_L, Lshell_conjunction, MLT_conjunction, Correct_0_counts, Correct_90_counts, Correct_mean_counts, BLCflux, timeindex):
    
    Ltitle = str(Lshell_conjunction)
    MLTtitle = str(MLT_conjunction)
    # Make panel plots
    fig, axs = plt.subplots(2)
    #fig.suptitle(title_plot+' L = '+Ltitle+' FB-MLT = '+MLTtitle)
    fig.suptitle(FBfig_text+'\n '+FBfig_text_latlon+'\n L = '+Ltitle+', MLT = '+MLTtitle)
    fig.set_size_inches(horz_size, vert_size_L)
    fig.frameon=False
    colorsFB = ['blue','red','green']
    colorsPOES = ['red','darkblue', 'blue','cornflowerblue','darkgray','gray','lightgray']
    #colorsPOES = ['#D3D3D3','lightblue','pink','lightgreen','blue', 'red', 'green']
    colorsflux = ['orange','red','purple','green','cyan']
    axs[0].set_prop_cycle(color=colorsflux)
    axs[1].set_prop_cycle(color=colorsPOES)
    
    # Plot FIREBIRD flux
    axs[0].plot(Lshell,flux_all, marker='.', markersize='3', linestyle='none')
    axs[0].set_yscale('log')
    axs[0].set_ylim([0.01,1.e2])
    axs[0].set_yticks([1.e-2,1.e-1, 1.e0, 1.e1, 1.e2])
    #locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=5)
    #axs[0].yaxis.set_minor_locator(locmin)
    axs[0].set_xlim([3,7])
    axs[0].set_xlabel('L shell')
    axs[0].set_ylabel('FIREBIRD-II \n '+FBsc+' Flux')
    axs[0].legend((energy_legend), loc='upper right', fontsize = 'small')
    #axs[0].text(3.1, 30, FBfig_text, {'color': 'k', 'fontsize': 8})
    #axs[0].text(3.1, 10, FBfig_text_latlon, {'color': 'k', 'fontsize': 8})
    #axs[0].text(3.1, 3, FBfig_txt_conj,{'color': 'k', 'fontsize': 8})
    axs[0].axvline(x=Lshell_conjunction, color='k')
    axs[0].grid(True)
    
    # Plot FIREBIRD equivalent counts
    #axs[1].plot(Lshell,FB_POES_counts[:,2], L_T89_POES[indexstart:indexend], POES_obs_counts[indexstart:indexend,2],L_T89_POES[indexstart:indexend], POES_Gmean_obs_counts[indexstart:indexend,2],L_T89_POES[indexstart:indexend], POES_90_obs_counts[indexstart:indexend,2])
    axs[1].plot(Lshell,FB_POES_counts[:,2], L_T89_POES[indexstart:indexend], Correct_0_counts[indexstart:indexend,1], L_T89_POES[indexstart:indexend], Correct_mean_counts[indexstart:indexend,1], L_T89_POES[indexstart:indexend], Correct_90_counts[indexstart:indexend,1], L_T89_POES[indexstart:indexend], POES_obs_counts[indexstart:indexend,2], L_T89_POES[indexstart:indexend], POES_Gmean_obs_counts[indexstart:indexend,2], L_T89_POES[indexstart:indexend], POES_90_obs_counts[indexstart:indexend,2])    #pettit May 2023
    axs[1].set_yscale('log')
    axs[1].set_ylim([1.,1.e4])
    axs[1].set_yticks([1., 1.e1, 1.e2, 1.e3, 1.e4])
    axs[1].set_xlim([3,7])
    axs[1].set_xlabel('L shell')
    axs[1].set_ylabel('Predicted \n Counts (>300 keV) \n  ')
    axs[1].legend(('FB', '0 tel', 'geo mean','90 tel'), loc='upper right', fontsize = 'small')
    axs[1].axvline(x=Lshell_conjunction, color='k')
    axs[1].grid(True)

    
    fig.autofmt_xdate()
    #plt.savefig(filename_plot_L, orientation='landscape', format='ps')   #MLT
    #plt.savefig(filename_plot_L, orientation='portrait', format = 'ps')   #MLT
    plt.savefig(filename_plot_L_all, format = 'ps', bbox_inches='tight')
#plt.show()



def create_plot_time(Lshell, timestamp, flux_all, energy_legend, FBfig_text, FBfig_text_latlon, FB_POES_counts, starttime_plot, endtime_plot, L_T89_POES, MLT_T89_POES, POES_time, POES_obs_counts,POESsc,POESfig_text, POESfig_text_latlon, POES_90_obs_counts, POES_Gmean_obs_counts, filename_plot_time, MLT_FB, Lshell_conjunction, MLT_conjunction, Correct_0_counts, Correct_90_counts, Correct_mean_counts, BLCflux, timeindex):
    
    POES_0_all = np.zeros((timeindex,7),'f')
    POES_90_all = np.zeros((timeindex,7),'f')
    POES_mean_all = np.zeros((timeindex,7),'f')
    POES_0_all[:,0:4] = POES_obs_counts[:,:]
    POES_0_all[:,4:8] = Correct_0_counts[:,:]
    POES_90_all[:,0:4] = POES_90_obs_counts[:,:]
    POES_90_all[:,4:8] = Correct_90_counts[:,:]
    POES_mean_all[:,0:4] = POES_Gmean_obs_counts[:,:]
    POES_mean_all[:,4:8] = Correct_mean_counts[:,:]
    
    Ltitle = str(Lshell_conjunction)
    MLTtitle = str(MLT_conjunction)
    # Make panel plots
    #fig, axs = plt.subplots(5)  #MLT
    fig, axs = plt.subplots(7)  #MLT
    #fig.suptitle(title_plot+' L = '+Ltitle+' FB-MLT = '+MLTtitle)
    fig.suptitle(FBfig_text+'\n '+FBfig_text_latlon+'\n L = '+Ltitle+', MLT = '+MLTtitle)
    fig.set_size_inches(horz_size, vert_size_time)
    fig.frameon=False
    colorsFB = ['blue','red','green']
    #colorsPOES = ['#D3D3D3','blue', 'red', 'green']
    colorsPOES = ['#D3D3D3','lightblue','pink','lightgreen','blue', 'red', 'green']
    colorsflux = ['orange','red','purple','green','cyan']
    colors_L = ['purple','orange']
    axs[0].set_prop_cycle(color=colorsflux)
    axs[1].set_prop_cycle(color=colorsFB)
    axs[2].set_prop_cycle(color=colorsPOES)
    axs[3].set_prop_cycle(color=colorsPOES)
    axs[4].set_prop_cycle(color=colors_L)
    axs[5].set_prop_cycle(color=colorsflux)
    axs[6].set_prop_cycle(color=colorsflux)
    #myFmt = mdates.DateFormatter('%H:%M')
    
    # Plot FIREBIRD flux
    axs[0].plot(timestamp,flux_all, marker='.', markersize='3', linestyle='none')
    axs[0].set_yscale('log')
    axs[0].set_ylim([0.01,1.e2])
    axs[0].set_yticks([1.e-2,1.e-1, 1.e0, 1.e1, 1.e2])
    #locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=5)
    #axs[0].yaxis.set_minor_locator(locmin)
    axs[0].set_xlim([starttime_plot, endtime_plot])
    axs[0].set_xlabel('time')
    axs[0].set_ylabel('FIREBIRD-II \n '+FBsc+' Flux')
    axs[0].legend((energy_legend), loc='upper right', fontsize = 'small')
    axs[0].grid(True)
    #axs[0].text(3.1, 20, FBfig_text, {'color': 'k', 'fontsize': 8})
    #axs[0].text(3.1, 5, FBfig_text_latlon, {'color': 'k', 'fontsize': 8})
    
    # Plot FIREBIRD equivalent counts
    axs[1].plot(timestamp,FB_POES_counts[:,1:4], marker='.', markersize='3', linestyle='none')
    axs[1].set_yscale('log')
    axs[1].set_ylim([1.,1.e4])
    axs[1].set_yticks([1., 1.e1, 1.e2, 1.e3, 1.e4])
    #locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=5)
    #axs[1].yaxis.set_minor_locator(locmin)
    axs[1].set_xlim([starttime_plot, endtime_plot])
    axs[1].set_xlabel('time')
    axs[1].set_ylabel('Predicted \n Counts \n  ')
    axs[1].legend(('> 100 keV', '> 300 keV', '> 600 keV'), loc='upper right', fontsize = 'small')
    axs[1].grid(True)
    
    # Plot POES counts tel0
    #axs[2].plot(POES_time, POES_obs_counts)
    axs[2].plot(POES_time, POES_0_all)
    axs[2].set_xlim([starttime_plot, endtime_plot])
    axs[2].set_yscale('log')
    axs[2].set_ylim([1.,1.e4])
    axs[2].set_yticks([1., 1.e1, 1.e2, 1.e3, 1.e4])
    axs[2].set_xlabel('time')
    axs[2].set_ylabel(POESsc+'\n tel0 \n counts \n  ')
    axs[2].legend(('> 40 keV', '> 100 keV', '> 300 keV', '> 600 keV'),loc='upper right', fontsize = 'small')
    #axs[2].text(3.1, 2000, POESfig_text, {'color': 'k', 'fontsize': 8})
    #axs[2].text(3.1, 400, POESfig_text_latlon, {'color': 'k', 'fontsize': 8})
    axs[2].grid(True)
 
    # Plot POES counts tel90
    #axs[3].plot(POES_time, POES_90_obs_counts)
    axs[3].plot(POES_time, POES_90_all)
    axs[3].set_xlim([starttime_plot, endtime_plot])
    axs[3].set_yscale('log')
    axs[3].set_ylim([1.,1.e4])
    axs[3].set_yticks([1., 1.e1, 1.e2, 1.e3, 1.e4])
    axs[3].set_xlabel('time')
    axs[3].set_ylabel(POESsc+'\n tel90 \n counts \n   ')
    axs[3].legend(('40 keV', '> 100 keV', '> 300 keV', '> 600 keV'), loc='upper right', fontsize = 'small')
    axs[3].grid(True)
    
    '''
        # Plot Corrected POES counts
        #axs[4].plot(POES_time[indexstartcor:indexendcor], E_correct[indexstartcor:indexendcor,:], marker=".")
        axs[4].plot(L_T89_POES[indexstartcor:indexendcor], E_correct[indexstartcor:indexendcor,:], marker=".")
        #axs[4].set_xlim([dt.datetime(2018,9,22,19,43,00), dt.datetime(2018,9,22,19,55,30)])axs[4].set_xlim([3,7])
        axs[4].set_yscale('log')
        axs[4].set_ylim([1.,1.e4])
        axs[4].set_yticks([1., 1.e1, 1.e2, 1.e3])
        #axs[4].set_xlabel('time')
        axs[4].set_xlabel('L shell')
        axs[4].set_ylabel('corrected')
        axs[4].legend(('> 40 keV', '> 100 keV', '> 300 keV', '> 600 keV'),
        loc='upper right', fontsize = 'small')
        '''
    
    #Plot L-shell
    # mask points where y < 2
    ymask = np.ma.masked_where((L_T89_POES < 2), L_T89_POES)
    axs[4].plot(POES_time, ymask, color = 'purple')
    axs[4].plot(timestamp,Lshell, color = 'orange')
    axs[4].set_xlim([starttime_plot, endtime_plot])
    axs[4].set_ylim([3.,10.])
    axs[4].set_yticks([4., 6., 8., 10.])
    axs[4].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:00"))
    axs[4].xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M:00"))
    axs[4].set_xlabel('time')
    axs[4].set_ylabel('L-shell')
    axs[4].legend(('POES L-shell', 'FB L-shell'), loc='upper right', fontsize = 'small')
    axs[4].grid(True)
    
    
    axs[5].plot(POES_time, BLCflux)
    axs[5].set_xlim([starttime_plot, endtime_plot])
    axs[5].set_yscale('log')
    axs[5].set_ylim([0.01,1.e3])
    axs[5].set_yticks([1.e-2,1.e-1, 1.e0, 1.e1, 1.e2, 1.e3])
    axs[5].set_xlabel('time')
    axs[5].set_ylabel('POES BLC_flux')
    axs[5].legend((energy_legend), loc='upper right', fontsize = 'small')
    axs[5].grid(True)
    
    FB_BLCflux = flux_all*3.758


    # Plot FIREBIRD BLC flux
    axs[6].plot(timestamp,FB_BLCflux, marker='.', markersize='3', linestyle='none')
    axs[6].set_yscale('log')
    axs[6].set_ylim([0.01,1.e3])
    axs[6].set_yticks([1.e-2,1.e-1, 1.e0, 1.e1, 1.e2, 1.e3])
    axs[6].set_xlim([starttime_plot, endtime_plot])
    axs[6].set_xlabel('time')
    axs[6].set_ylabel('FIREBIRD-II \n '+FBsc+'BLC Flux')
    axs[6].legend((energy_legend), loc='upper right', fontsize = 'small')
    axs[6].grid(True)



    '''
        #MLT start
        ymaskMLT = np.ma.masked_where((L_T89_POES < 2), MLT_T89_POES)
        axs[5].plot(POES_time, ymaskMLT, color='purple')
        axs[5].plot(timestamp,MLT_FB,color='orange')
        axs[5].set_xlim([starttime_plot, endtime_plot])
        axs[5].set_ylim([0.,24.])
        axs[5].set_yticks([0., 6., 12., 18., 24.])
        axs[5].set_xlabel('time')
        axs[5].xaxis.set_major_formatter(myFmt)
        axs[5].set_ylabel('MLT')
        axs[5].legend(('POES MLT', 'FB MLT'), loc='upper right', fontsize='small')
        '''
    
    fig.autofmt_xdate()
    #plt.savefig(filename_plot_time, orientation='landscape', format = 'ps')  #MLT
    plt.savefig(filename_plot_time, orientation='portrait', format='ps', bbox_inches='tight')
    #plt.show()

def integral_efolding(E1, E2, J0, E0):
    return J0 * E0 * (np.exp(-E1/E0) - np.exp(-E2/E0))

def main():
    print("Read POES Gfactors")
    table4, table2 = read_POES_Gfactors(filename_POES_GfactorsE, filename_POES_GfactorsP)
    print("Read FB Flux")
    timestamp, J0, E0, Lshell, MLT_FB, nlines = read_FB_flux(filename_FB_flux, filename_FB_Lshell)    #MLT
    print("Read POES observations counts")
    POES_time, POES_obs_counts, POES_90_obs_counts, POES_Gmean_obs_counts, ntimePOES = read_POES_counts(filename_POES_counts, nlines)
    #timecorrect, E_correct, ntime_correct = read_POES_corrected(filename_correct)
    print("Read_POES_Lshell")
    L_T89_POES, MLT_T89_POES = read_POES_Lshell(filename_POES_L)
    
    FB_POES_counts = np.zeros((nlines,4),dtype=float)
    flux_all = np.zeros((nlines,5), dtype = float)
    FB_BLC_flux = np.zeros((nlines,5), dtype = float)

    print("Predicting what POES should see ")
    FB_POES_counts, flux_all = calculate_counts(nlines, energy_FB, J0, E0, table4, table2)
    open(filename_values, 'w')
    print("Finding conjuction time index for FB and POES")
    FBindexconj, FB_conj_L, FB_conj_MLT = FB_conj_index(nlines, timestamp, J0, E0, Lshell, MLT_FB, conj_time, filename_values)
    indexstart, indexend, indexconj = time_interval(ntimePOES, startPOES, endPOES, conj_time, POES_time, L_T89_POES)
    
    # Read Josh Pettit's's POES file - July 2022
    print("Read Josh Pettit's POES file")#
    #POES0_correct_counts, POES90_correct_counts, POESGmean_correct_counts, POES_BLC_flux = read_Pettit_POES_counts(filename_Pettit_nc, POES_time, energy_FB)
    
    # Read Josh Pettit's's POES file - May 2023
    #       Need to read BLC Flux from different file
    #       2 second timing
    #       "time" instead of "m_time"
    #       double EOcounts_corrected(time, electron_telescopes_and_E4)
    
    print("Read Josh Pettit's POES file May 2023")
    POES0_correct_counts, POES90_correct_counts, POESGmean_correct_counts, POES_BLC_flux = read_Pettit_POES_counts(filename_Pettit_0_nc, filename_Pettit_90_nc, filename_Pettit_BLC_nc, POES_time, energy_FB, indexstart, indexend)

    print("POES observed counts shape and dtype ")
    print(POES_obs_counts.shape, POES_obs_counts.dtype)

    print("correct counts shape and dtype ")
    print(POES0_correct_counts.shape, POES0_correct_counts.dtype)
    
    print("save output files")
    save_POES_equiv_counts(filename_values, nlines, FBindexconj, Lshell, FB_POES_counts)
    POES_tel = 'POES OBSERVED COUNTS 0 TELESCOPE'
    save_POES_obs(filename_values, indexstart, indexend, indexconj, POES_obs_counts, POES0_correct_counts, L_T89_POES, POES_tel)
    POES_tel = 'POES OBSERVED COUNTS 90 TELESCOPE'
    save_POES_obs(filename_values, indexstart, indexend, indexconj, POES_90_obs_counts, POES90_correct_counts, L_T89_POES, POES_tel)
    POES_tel = 'POES OBSERVED COUNTS GEOMETRIC MEAN'
    save_POES_obs(filename_values, indexstart, indexend, indexconj, POES_Gmean_obs_counts, POESGmean_correct_counts, L_T89_POES, POES_tel)
    
    print("save FB and POES flux")
    subtitle = ' --------  POES BOUNCE LOSS CONE FLUX -------- \n'
    save_flux(filename_values, indexstart, indexend, indexconj, POES_BLC_flux, L_T89_POES, subtitle, energy_legend, POES_time)

    FB_BLCflux = flux_all*3.758
    
    subtitle = ' --------  FB BOUNCE LOSS CONE FLUX -------- \n'
    len_FBflux = len(flux_all)-1
    save_flux(filename_values, 0, len_FBflux, FBindexconj, FB_BLCflux, Lshell, subtitle, energy_legend, timestamp)

    equiv_count_file = open(filename_output, 'w')
    print("date time, L, E0 equiv, E1 equiv, E2 equiv, P6 equiv",  file = equiv_count_file)
    np.savetxt(equiv_count_file, np.column_stack((timestamp, Lshell, FB_POES_counts)), fmt = '%s', delimiter=', ')
    #np.savetxt(filename_output, np.column_stack((timestamp, Lshell, FB_POES_counts)), fmt = '%s', delimiter=', ')
    poes_count_file_0tel = open(filename_output_POES_0tel,'w')
    print("date time, L, MLT, E0, E1, E2, P6",  file = poes_count_file_0tel)
    np.savetxt(poes_count_file_0tel, np.column_stack((POES_time[indexstart:indexend], L_T89_POES[indexstart:indexend], MLT_T89_POES[indexstart:indexend], POES_obs_counts[indexstart:indexend,:])), fmt = '%s', delimiter=', ')
    poes_count_file_90tel = open(filename_output_POES_90tel,'w')
    print("date time, L, MLT, E0, E1, E2, P6",  file = poes_count_file_90tel)
    np.savetxt(poes_count_file_90tel, np.column_stack((POES_time[indexstart:indexend], L_T89_POES[indexstart:indexend], MLT_T89_POES[indexstart:indexend], POES_90_obs_counts[indexstart:indexend,:])), fmt = '%s', delimiter=', ')
    poes_count_file_Gmean = open('OUTPUT/POES_counts_geometric_mean.txt','w')
    print("date time, L, MLT, E0, E1, E2, P6",  file = poes_count_file_Gmean)
    np.savetxt(poes_count_file_Gmean, np.column_stack((POES_time[indexstart:indexend], L_T89_POES[indexstart:indexend], MLT_T89_POES[indexstart:indexend], POES_Gmean_obs_counts[indexstart:indexend,:])), fmt = '%s', delimiter=', ')

    print("create plots")
    create_plot_Lshell(Lshell, flux_all, energy_legend,FBfig_text, FBfig_text_latlon, FB_POES_counts, indexstart, indexend, L_T89_POES, POES_obs_counts,POESsc,POESfig_text, POESfig_text_latlon, POES_90_obs_counts, POES_Gmean_obs_counts, filename_plot_L, FB_conj_L, FB_conj_MLT, POES0_correct_counts,POES90_correct_counts, POESGmean_correct_counts,POES_BLC_flux, ntimePOES)
    create_plot_time(Lshell, timestamp, flux_all, energy_legend,FBfig_text, FBfig_text_latlon, FB_POES_counts, starttime_plot, endtime_plot, L_T89_POES, MLT_T89_POES, POES_time, POES_obs_counts,POESsc,POESfig_text, POESfig_text_latlon, POES_90_obs_counts, POES_Gmean_obs_counts, filename_plot_time, MLT_FB, FB_conj_L, FB_conj_MLT, POES0_correct_counts,POES90_correct_counts, POESGmean_correct_counts,POES_BLC_flux,ntimePOES)   #MLT
    create_plot_Lshell_all(Lshell, flux_all, energy_legend,FBfig_text, FBfig_text_latlon, FB_POES_counts, indexstart, indexend, L_T89_POES, POES_obs_counts,POESsc,POESfig_text, POESfig_text_latlon, POES_90_obs_counts, POES_Gmean_obs_counts, filename_plot_L, FB_conj_L, FB_conj_MLT, POES0_correct_counts,POES90_correct_counts, POESGmean_correct_counts,POES_BLC_flux, ntimePOES)


main()








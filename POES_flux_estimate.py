'''
Filename: 
    POES_flux_estimate.py
    
Author: 
    Isabella M. Householder
    Modified by K. Duderstadt for POES
    
Date:
    05 September 2019
    Modified 27 Nov 2019 by K. Duderstadt to run through shell (sh FB-POES_counts.sh)
    Modified 12 Sep 2020 by K. Duderstadt to estimate differential flux from POES integral flux
    
Description:
    This program calculates the minimum parameters from the best curve fit of
    POES flux data.
    
Input: 
    poes_m02_20181231_proc.nc

Output: 
    flux-params_m02_2018-12-31.txt
    
To Run:
    python POES_flux_estimate.py
    
Modify:
    define filename and outfile in lines xx-xx
    
-------------------------------------------------------------------------------
'''

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import math
from math import log10
import datetime as dt
from datetime import datetime, timedelta
import netCDF4 as nc
from netCDF4 import Dataset
import sys
import ast

filename_POES_flux = 'DATA/POES/poes_m02_20181231_proc.nc'
outfile = 'OUTPUT/flux-params_poes_m02_20181231_diff_flux.txt'

startPOES=dt.datetime(2018,12,31,7,21,00)       #POES times to plot 3-7
endPOES=dt.datetime(2018,12,31,7,28,00)         #POES times to plot 3-7
conj_time = dt.datetime(2018,12,31,7,22,30)   # if a range, choose something in the middle

channels = np.array([40., 130., 287., 612.], dtype=float) #POES integral channels
#channel_middles = np.array([251.5, 333.5, 452.0, 620.5, 852.8], dtype=float) #FU4 energies
channel_middles = np.array([265.4,353.7,481.2,662.7,913.0], dtype=float) #FU3 energies


def read_POES_flux(filename_flux_nc):
    nc_fid = nc.Dataset(filename_flux_nc, mode = 'r')
    e1_0_flux_POES = nc_fid.variables['mep_ele_tel0_flux_e1'][:]
    e2_0_flux_POES = nc_fid.variables['mep_ele_tel0_flux_e2'][:]
    e3_0_flux_POES = nc_fid.variables['mep_ele_tel0_flux_e3'][:]
    e4_0_flux_POES = nc_fid.variables['mep_ele_tel0_flux_e4'][:]
    e1_90_flux_POES = nc_fid.variables['mep_ele_tel90_flux_e1'][:]
    e2_90_flux_POES = nc_fid.variables['mep_ele_tel90_flux_e2'][:]
    e3_90_flux_POES = nc_fid.variables['mep_ele_tel90_flux_e3'][:]
    e4_90_flux_POES = nc_fid.variables['mep_ele_tel90_flux_e4'][:]
    timeindex = len(e1_0_flux_POES)
    print("timeindex = ", timeindex)
    print(e1_0_flux_POES[0:10])
    print("formatting time for POES...this takes a couple minutes")
    timeclock = [dt.datetime.strptime(f'{nc_fid["year"][i]}, {nc_fid["day"][i]}', '%Y, %j') +timedelta(milliseconds=int(nc_fid['msec'][i])) for i in range(len(nc_fid['time'][:]))]
    e1_0_flux = np.zeros((timeindex),'f')
    e2_0_flux = np.zeros((timeindex),'f')
    e3_0_flux = np.zeros((timeindex),'f')
    e4_0_flux = np.zeros((timeindex),'f')
    e1_90_flux = np.zeros((timeindex),'f')
    e2_90_flux = np.zeros((timeindex),'f')
    e3_90_flux = np.zeros((timeindex),'f')
    e4_90_flux = np.zeros((timeindex),'f')

    e1_0_flux[:] = e1_0_flux_POES[:]
    e2_0_flux[:] = e2_0_flux_POES[:]
    e3_0_flux[:] = e3_0_flux_POES[:]
    e4_0_flux[:] = e4_0_flux_POES[:]
    e1_90_flux[:] = e1_90_flux_POES[:]
    e2_90_flux[:] = e2_90_flux_POES[:]
    e3_90_flux[:] = e3_90_flux_POES[:]
    e4_90_flux[:] = e4_90_flux_POES[:]
    POES_obs = np.zeros((timeindex,4),'f')
    POES_obs_90 = np.zeros((timeindex,4),'f')
    POES_obs[:,0] = e1_0_flux
    POES_obs[:,1] = e2_0_flux
    for ntime in range(timeindex):
        if e1_0_flux[ntime] <= 0.:
            POES_obs[ntime, 0] = 'nan'
        if e2_0_flux[ntime] <= 0.:
            POES_obs[ntime, 1] = 'nan'
        if e3_0_flux[ntime] > 250.:
            POES_obs[ntime,2] = e3_0_flux[ntime]
            if e4_0_flux[ntime] > 250.:
                POES_obs[ntime,3] = e4_0_flux[ntime]
            else:
                POES_obs[ntime,3] = 'nan'
        else:
            POES_obs[ntime,2] = 'nan'
            POES_obs[ntime,3] = 'nan'
    POES_obs_90[:,0] = e1_90_flux
    POES_obs_90[:,1] = e2_90_flux
    for ntime in range(timeindex):
        if e1_90_flux[ntime] <= 0.:
            POES_obs_90[ntime, 0] = 'nan'
        if e2_90_flux[ntime] <= 0.:
            POES_obs_90[ntime, 1] = 'nan'
        if e3_90_flux[ntime] > 250.:
            POES_obs_90[ntime,2] = e3_90_flux[ntime]
            if e4_90_flux[ntime] > 250.:
                POES_obs_90[ntime,3] = e4_90_flux[ntime]
            else:
                POES_obs_90[ntime,3] = 'nan'
        else:
            POES_obs_90[ntime,2] = 'nan'
            POES_obs_90[ntime,3] = 'nan'
    nc_fid.close()
    return timeclock, POES_obs, POES_obs_90, timeindex


def func(E, J0, E0):
    return J0 * np.exp(-E/E0)

def integral_efolding(E1, J0, E0):
    return J0 * E0 * (np.exp(-E1/(1000.*E0)) - np.exp(-1000./(1000.*E0)))


# identify time interval for POES plots
def time_interval(POESntime, start, end, conj, POEStime):
    flag1 = 0
    flag2 = 0
    flag3 = 0
    ntime = 0
    index_start = 0
    index_end = 0
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
    print("Conjuntion time index for POES = ", index_conj)
    return index_start, index_end, index_conj

       
def main():
    
    print("Read POES flux counts")
    POES_time, POES_obs_flux, POES_90_obs_flux, ntimePOES = read_POES_flux(filename_POES_flux)
    print("Finding conjuction time index for FB and POES")
    indexstart, indexend, indexconj = time_interval(ntimePOES, startPOES, endPOES, conj_time, POES_time)

    data = np.zeros((indexend-indexstart, 2), dtype=float)
    J_flux= np.zeros((indexend-indexstart, 5), dtype=float)
    #for timestep in range(n_timesteps):
    #for timestep in range(indexstart, indexend):
    
    print(POES_90_obs_flux[:,:])
    for timestep in range(indexstart, indexend):
        print("timestep = ",timestep)
        print("channels = ",channels[0:4])
        print("POES_90_obs_flux = ",POES_90_obs_flux[timestep,0:4])
        
        #---- curve fit ----#
        ## curve fit -- use 1000 keV (1 MeV) as upper bound
        
        if not math.isnan(POES_90_obs_flux[timestep,0]):
            if not math.isnan(POES_90_obs_flux[timestep,1]):
                if not math.isnan(POES_90_obs_flux[timestep,2]):
                    if not math.isnan(POES_90_obs_flux[timestep,3]):
                        print("all data")
                        popt, pcov = curve_fit(integral_efolding, channels[0:4], POES_90_obs_flux[timestep,0:4])
                    else:
                        print("> 600 keV is nan")
                        popt, pcov = curve_fit(integral_efolding, channels[0:3], POES_90_obs_flux[timestep,0:3])
                else:
                    popt, pcov = curve_fit(integral_efolding, channels[0:2], POES_90_obs_flux[timestep,0:2])
            else:
                popt[0] == 'nan'
                popt[1] == 'nan'
        else:
            popt[0] == 'nan'
            popt[1] == 'nan'
        coefficients = (popt[0], popt[1]*1000.)
        print('time J0 E0 :', timestep, coefficients)
            
            
    print(POES_obs_flux[:,:])
    for timestep in range(indexstart, indexend):
        print("timestep = ",timestep)
        print("channels = ",channels[0:4])
        print("POES_0_obs_flux = ",POES_obs_flux[timestep,0:4])
        
        #---- curve fit ----#
        ## curve fit -- use 1000 keV (1 MeV) as upper bound
        
        if not math.isnan(POES_obs_flux[timestep,0]):
            if not math.isnan(POES_obs_flux[timestep,1]):
                if not math.isnan(POES_obs_flux[timestep,2]):
                    if not math.isnan(POES_obs_flux[timestep,3]):
                        print("all data")
                        popt, pcov = curve_fit(integral_efolding, channels[0:4], POES_obs_flux[timestep,0:4],bounds=([0,0], [1.e5, 1.e3]))
                    else:
                        print("> 600 keV is nan")
                        popt, pcov = curve_fit(integral_efolding, channels[0:3], POES_obs_flux[timestep,0:3],bounds=([0,0], [1.e5, 1.e3]))
                else:
                    popt, pcov = curve_fit(integral_efolding, channels[0:2], POES_obs_flux[timestep,0:2],bounds=([0,0], [1.e5, 1.e3]))
            else:
                popt[0] == 'nan'
                popt[1] == 'nan'
        else:
            popt[0] == 'nan'
            popt[1] == 'nan'
        print("popt ", popt)
        popt[1] = popt[1]*1000.
        coefficients = (popt[0], popt[1])
        #print('time J0 E0 :', POES_time[timestep+indexstart], coefficients)
            
        data[timestep-indexstart,0] = popt[0]
        data[timestep-indexstart,1] = popt[1]
        
        J_flux[timestep-indexstart,:] = func(channel_middles, *popt )

    full_data = (POES_time[indexstart:indexend], data, J_flux)
    np.savetxt(outfile, np.column_stack(full_data), fmt='%s', delimiter=', ')
    print('Complete.')

main()


        






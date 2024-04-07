'''
Filename: 
    FB_flux_estimate_shell.py
    
Author: 
    Isabella M. Householder
    
Date:
    05 September 2019
    Modified 27 Nov 2019 by K. Duderstadt to run through shell (sh FB-POES_counts.sh)
    
Description:
    This program calculates the minimum parameters from the best curve fit of
    FIREBIRD high-resolution data.
    
Input: 
    FU3_Gfactors.txt
    sec-avg_FU3_2018-12-31_0720-0724.txt
    FU3-col.txt

Output: 
    flux-params_FU3_2018-12-31_0720-0724.txt
    
To Run:
    python FB_datetime_average_shell.py
    
Modify:
    define filename and outfile in lines 55-62
    
-------------------------------------------------------------------------------
'''

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from math import log10
import sys
import ast

filename1=sys.argv[1]   #$Gfactors
filename2=sys.argv[2]   #$FBavg
filename3=sys.argv[3]   #$GEANT
outfile=sys.argv[4]     #$fluxpar
outfile=sys.argv[4]           #$fluxpar
energy_FB_width=sys.argv[5]   #$energy_FB_width
energy_FB=sys.argv[6]         #$energy_FB

channel_width_list = [float(s) for s in energy_FB_width.strip('[]').split(',')]
channel_width = np.array(channel_width_list, dtype=float)
channel_middles_list = [float(s) for s in energy_FB.strip('[]').split(',')]
channel_middles = np.array(channel_middles_list, dtype=float)

'''
filename1 = 'G-factors/FU4_Gfactors.txt'
filename2 = 'OUTPUT/sec-avg_FU4_2018-09-27_0135-0142.txt'
filename3 = 'G-factors/FU4-col.txt'
outfile = 'OUTPUT/flux-params_FU4_2018-09-27_0135-0142.txt'
#channel_width = np.array([68.7, 107.9, 147.2, 215.9, 284.5], dtype=float)     #FU3
#channel_middles = np.array([265.4, 353.7, 481.2, 662.7, 913.0], dtype=float)  #FU3
channel_width = np.array([63.7, 100.2, 136.7, 200.4, 264.25], dtype=float)   #FU4
channel_middles = np.array([251.5, 333.5, 452.0, 620.5, 852.8], dtype=float) #FU4

#channel_middles = np.divide(channel_middles, 1000.)
'''

def read_file1(filename1):
    GfactorList = []
    try:
        file = open(filename1, 'r')
    except:
        print('The file may not exist or the program may have not have been able to open it.')
        return([])
    else:
        file = open(filename1, 'r')
        for line in file:
            line = line.strip().split(',')
            GfactorList.append(line)
    file.close()
    return GfactorList

def read_file2(filename2):
    countList = []
    try:
        file = open(filename2, 'r')
    except:
        print('The file may not exist or the program may have not have been able to open it.')
        return([])
    else:
        file = open(filename2, 'r')
        for line in file:
            line = line.strip().split(',')
            countList.append(line)
    file.close()
    return countList

def read_file3(filename3):
    colList = []
    try:
        file = open(filename3, 'r')
    except:
        print('The file may not exist or the program may have not have been able to open it.')
        return([])
    else:
        file = open(filename3, 'r')
        for line in file:
            line = line.strip().split(',')
            colList.append(line[0])
    file.close()
    colArray = np.array(colList, dtype=float)
    return colArray

def calculate_flux(GfactorList, countList, channel_width):
    cadence = 0.05
    G = 6.0
    timestepList = []
    counts0 = []
    counts1 = []
    counts2 = []
    counts3 = []
    counts4 = []
    
    for row in countList:
        timestepList.append(row[0])
        counts0.append(row[1])
        counts1.append(row[2])
        counts2.append(row[3])
        counts3.append(row[4])
        counts4.append(row[5])
    
    counts0 = np.array(counts0, dtype=float)
    counts1 = np.array(counts1, dtype=float)
    counts2 = np.array(counts2, dtype=float)
    counts3 = np.array(counts3, dtype=float)
    counts4 = np.array(counts4, dtype=float)
    
    ## --- measured counts array ---
    
    dim_counts = len(counts0)
    measuredcountsArray = np.zeros((dim_counts, 5), dtype=float)
    measuredcountsArray[:,0] = counts0[:]
    measuredcountsArray[:,1] = counts1[:]  
    measuredcountsArray[:,2] = counts2[:]
    measuredcountsArray[:,3] = counts3[:]
    measuredcountsArray[:,4] = counts4[:]      
    
    ## -------
    
    flux0 = counts0 / (cadence * G * channel_width[0])
    flux1 = counts1 / (cadence * G * channel_width[1])
    flux2 = counts2 / (cadence * G * channel_width[2])
    flux3 = counts3 / (cadence * G * channel_width[3])
    flux4 = counts4 / (cadence * G * channel_width[4])
    
    ntimesteps = np.size(flux0)
    
    fluxArray = np.zeros((ntimesteps, 5), dtype=float)
    fluxArray[:,0] = flux0
    fluxArray[:,1] = flux1
    fluxArray[:,2] = flux2
    fluxArray[:,3] = flux3
    fluxArray[:,4] = flux4

    return flux0, flux1, flux2, flux3, flux4, fluxArray, measuredcountsArray, timestepList, ntimesteps
    #flux=counts/(cadence*G*channel_width_keV) --> flux calculation formula

def func(E, J0, E0):
    return J0 * np.exp(-E/E0)

def integral_efolding(E1, E2, J0, E0):
    return J0 * E0 * (np.exp(-E1/E0) - np.exp(-E2/E0))

def midpoints(colArray):
    midpointArray = np.zeros((98), dtype=float)
    colwidthArray = np.zeros((98), dtype=float)
    
    n = 0
    while n < 98:
        midpointArray[n] = (np.divide((colArray[n] + colArray[n+1]), 2.))
        colwidthArray[n] = (colArray[n+1] - colArray[n])
        n += 1
    return midpointArray, colwidthArray

def calculate_counts (GfactorList, colArray, *popt):
    
    #Jflux = func(midpointArray, *popt)
    En_old = colArray
    En = np.zeros((99), dtype=float)
    En[:98] = En_old[:98]
    En[98] = 2000.

    newcountArray = np.zeros((99,5), dtype=float)
    estimatedcountsArray = np.zeros((5), dtype=float)
    cadence = 0.05
    GfactorArray = np.array(GfactorList, dtype=float)
    
    n = 0
    while n < 98:
        #newcountArray[n,0] = Jflux[n] * cadence * GfactorArray[n, 0] * (colwidthArray[n] * 1000.)
        #newcountArray[n,1] = Jflux[n] * cadence * GfactorArray[n, 1] * (colwidthArray[n] * 1000.)
        #newcountArray[n,2] = Jflux[n] * cadence * GfactorArray[n, 2] * (colwidthArray[n] * 1000.)
        #newcountArray[n,3] = Jflux[n] * cadence * GfactorArray[n, 3] * (colwidthArray[n] * 1000.)
        #newcountArray[n,4] = Jflux[n] * cadence * GfactorArray[n, 4] * (colwidthArray[n] * 1000.)
        
        newcountArray[n,0] = integral_efolding(En[n], En[n+1], *popt) * cadence * GfactorArray[n, 0]
        newcountArray[n,1] = integral_efolding(En[n], En[n+1], *popt) * cadence * GfactorArray[n, 1]
        newcountArray[n,2] = integral_efolding(En[n], En[n+1], *popt) * cadence * GfactorArray[n, 2]
        newcountArray[n,3] = integral_efolding(En[n], En[n+1], *popt) * cadence * GfactorArray[n, 3]
        newcountArray[n,4] = integral_efolding(En[n], En[n+1], *popt) * cadence * GfactorArray[n, 4]
        
        estimatedcountsArray[0] = estimatedcountsArray[0] + newcountArray[n,0]
        estimatedcountsArray[1] = estimatedcountsArray[1] + newcountArray[n,1]
        estimatedcountsArray[2] = estimatedcountsArray[2] + newcountArray[n,2]
        estimatedcountsArray[3] = estimatedcountsArray[3] + newcountArray[n,3]
        estimatedcountsArray[4] = estimatedcountsArray[4] + newcountArray[n,4]
        n += 1     
    return estimatedcountsArray
    #### counts = flux * cadence * Gfactor * colwidth

def wrapper(coeff, measuredcountsArray, GfactorList, colArray, timestep):
    coeff = list(coeff)
    sigfigs = 2
    order = [int(log10(i)) for i in coeff]
    for s in range(sigfigs):
        for i in [0, 1]:
            coeff[i] = round(coeff[i], sigfigs - 1 - order[i])
            if sigfigs - 1 - order[i] <= 0:
                coeff[i] = int(coeff[i])
        iter_range = [10 ** (x - s) for x in order]
        iter_step = [x / 10 for x in iter_range]
    
        e0_range = iter_range[0]
        j0_range = iter_range[1]
        e0_step = iter_step[0]
        j0_step = iter_step[1]
    
        min_params, chi_sqr = iterate_spectrum(coeff, measuredcountsArray, GfactorList, colArray, e0_range, j0_range, e0_step, j0_step, timestep)
        
        coeff = min_params
        
    return coeff, chi_sqr

def calculate_steps(coeff, e0_range, j0_range, e0_step, j0_step):
    steps = [[], []]
    e0i = coeff[0]
    j0i = coeff[1]
    if e0i > e0_range:
        steps[0] = np.arange(e0i - e0_range, e0i + e0_range, e0_step)
    else:
        steps[0] = np.arange(e0_step, e0i + e0_range, e0_step)
    if j0i > j0_range:
         steps[1] = np.arange(j0i - j0_range, j0i + j0_range, j0_step)
    else:
        steps[1] = np.arange(j0_step, j0i + j0_range, j0_step)

    return steps
       

def iterate_spectrum(coeff, measuredcountsArray, GfactorList, colArray, e0_range, j0_range, e0_step, j0_step, timestep):
    max_iterations = 50
    done = False
    loops = 0
    while not done:  ##while done == True
        if loops >= max_iterations:
            #print('Maximum iterations reached')
            #flag = 1
            break
        steps = calculate_steps(coeff, e0_range, j0_range, e0_step, j0_step)
        chi_sqr = find_chi_sqr(steps, measuredcountsArray, GfactorList, colArray, timestep)  # Calculate chi squared for each point
        idx = np.unravel_index(np.argmin(chi_sqr), chi_sqr.shape)
        min_params = np.array([steps[0][idx[0]], steps[1][idx[1]]]) 
        #print('Current minimum parameters, rounded:')
        #print(min_params)
        #print('Current iteration: E0: ',coeff[0],'J0: ',coeff[1])
        loops += 1
        if [round(x, 5) for x in min_params] == [round(x, 5) for x in coeff]:
            done = True
            #iteration_flag = 0
        else:
            # set coefficeints to the minimum point for the next loop
            coeff = min_params
            loops += 1
            
    return min_params, chi_sqr
            
    
def find_chi_sqr(steps, measuredcountsArray, GfactorList, colArray, timestep):
    dimension = (len(steps[0]), len(steps[1])) 
    chi_sqr = np.zeros(dimension) 
    for i, p1 in enumerate(steps[0]): 
        for j, p2 in enumerate(steps[1]): 
            estimatedcountsArray = calculate_counts(GfactorList, colArray, p1, p2)
            chi_sqr[i, j] = np.sum((measuredcountsArray[timestep,:] - estimatedcountsArray) ** 2) 
    return chi_sqr 
    
       
def main():
    Gfactor_var = read_file1(filename1)
    counts_var = read_file2(filename2)
    col_var = read_file3(filename3)*1000.

    
    flux_0, flux_1, flux_2, flux_3, flux_4, flux_array, measuredcounts_array, timestep_list, n_timesteps = calculate_flux(Gfactor_var, counts_var, channel_width)
    
    midpoint_array, colwidth_array = midpoints(col_var)
    
    '''
        t = 10        #arbitrary timestep to plots
        plt.figure()
        plt.suptitle('Energy (MeV) vs. Energetic Electron Flux', fontsize=11)
        plt.xlabel('Energy (MeV)')
        plt.ylabel('Energetic Electron Flux particles cm^-2 s^-1 st^-1 keV^-1')
        plt.xticks((channel_middles))
        '''
    
    data = np.zeros((n_timesteps, 2), dtype=float)
    for timestep in range(n_timesteps):
    #for timestep in range(120,123):
        
        '''
            if timestep == t:
            plt.plot(channel_middles[0], flux_0[timestep], 'rx')
            plt.plot(channel_middles[1], flux_1[timestep], 'rx')
            plt.plot(channel_middles[2], flux_2[timestep], 'rx')
            plt.plot(channel_middles[3], flux_3[timestep], 'rx')
            plt.plot(channel_middles[4], flux_4[timestep], 'rx')
            '''
        
        #---- curve fit ----#
        ## curve fit
        counts_np = np.array(counts_var)
        counts_time = counts_np[timestep,1:6]
        counts_energy = counts_time.astype(np.float64)
        print("counts")
        print(counts_energy)
        print("first guess flux")
        print(flux_array[timestep,:])
        if all(counts_energy[0:4]) >= 1.:
           print ("Four highest energy counts > 1")
            
           popt, pcov = curve_fit(func, channel_middles/1000., flux_array[timestep,:], bounds=(-np.inf, np.inf))
        
           coefficients = (popt[0], popt[1]*1000.)
           print('Initial parameters:')
           print(coefficients)
           if coefficients[0] > 0. and coefficients[1] > 0.:
                min_params, chi_sqr = wrapper(coefficients, measuredcounts_array, Gfactor_var, col_var, timestep)
                
           #J_flux = func(midpoint_array, *popt)
           
                print('Minimum parameters:')
                print(timestep, min_params)
                data[timestep,0] = min_params[0]
                data[timestep,1] = min_params[1]
           else:
                data[timestep,0] = 'nan'
                data[timestep,1] = 'nan'
        else:
           data[timestep,0] = 'nan'
           data[timestep,1] = 'nan'

    full_data = (timestep_list, data)
    np.savetxt(outfile, np.column_stack(full_data), fmt='%s', delimiter=', ')
    print('Complete.')

main()


        






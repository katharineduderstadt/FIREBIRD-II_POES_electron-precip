'''
Filename:
    FB_select_times_shell.py

Authors:
    Isabella M. Householder
    Katharine Duderstadt
    
Date:
    24 November 2019
    
Description:
    This program reads in .txt file downloaded from http://solar.physics.montana.edu/FIREBIRD_II/Data/
    Save column counts and L shell (this step was previously done using Autoplot)
    
Input Files:
    TimeSlice_FU3_2018-12-31_0720-0724.txt
    
Output files:
    TimeSlice_FU3_2018-12-31_0720-0724.txt
    
To Run:
    python FB_select_times.py
    
Modify:
    define filenames and dates

-------------------------------------------------------------------------------
'''
import datetime as dt
import numpy as np
import sys

filename=sys.argv[1]   #$FBhires
outfile=sys.argv[2]    #$FBslice
startarg=sys.argv[3]   #$startFB
endarg=sys.argv[4]     #$endFB

sd=[int(s) for s in startarg.split(',')]
ed=[int(s) for s in endarg.split(',')]
starttime = dt.datetime(sd[0],sd[1],sd[2],sd[3],sd[4],sd[5])
endtime = dt.datetime(ed[0],ed[1],ed[2],ed[3],ed[4],ed[5])

'''
filename = 'DATA/FIREBIRD/FU4_Hires_2018-09-27_L2.txt'
outfile = 'OUTPUT/TimeSlice_FU4_2018-09-27_0135-0142.txt'
starttime = dt.datetime(2018,9,27,1,35,00)
endtime = dt.datetime(2018,9,27,1,39,00)
'''

print('start time ', starttime)
print('end time ', endtime)

#------------------- reads in file and creates list of data -------------------
def read_data(filename):
    dataList = []
    try:
        file = open(filename, 'r')
    except:
        print('The file may not exist or the program may have not have been able to open it.')
        return([])
    else:
        file = open(filename, 'r')
        header_count = 0
        for line in file:
            if line.startswith("#"):
               header_count = header_count+1
            else:
                line = line.strip().split(' ')
                dataList.append(line)
    file.close()
    return dataList

def select_time(dataList):
    time_list = []
    col1_counts_list = []
    col2_counts_list = []
    col3_counts_list = []
    col4_counts_list = []
    col5_counts_list = []
    col6_counts_list = []
    MLT_list = []
    L_shell_list = []
    flag1 = 0
    flag2 = 0
    indexstart = 0
    indexend = 0
    timestep = 0
    for row in dataList:
        if len(row[0]) == 19:
            row[0] = row[0]+'.000000'
        time = dt.datetime.strptime(row[0], '%Y-%m-%dT%H:%M:%S.%f')
        time_list.append(time)
        col1_counts_list.append(float(row[13]))
        col2_counts_list.append(float(row[14]))
        col3_counts_list.append(float(row[15]))
        col4_counts_list.append(float(row[16]))
        col5_counts_list.append(float(row[17]))
        col6_counts_list.append(float(row[18]))
        MLT_list.append(float(row[30]))
        L_shell_list.append(float(row[26]))
        if flag1 == 0 and time >= starttime:
            indexstart = timestep
            flag1 = 1
            print('indexstart', indexstart)
        if flag2 == 0 and time >= endtime:
            indexend = timestep
            flag2 = 2
            print('indexend', indexend)
        timestep = timestep+1
    range = indexend-indexstart
    col_counts = np.zeros((range,6), 'f')
    L_shell = np.zeros((range), 'f')
    MLT = np.zeros((range), 'f')
    time_new = []
    time_new[:] = time_list[indexstart:indexend]
    col_counts[:,0] = col1_counts_list[indexstart:indexend]
    col_counts[:,1] = col2_counts_list[indexstart:indexend]
    col_counts[:,2] = col3_counts_list[indexstart:indexend]
    col_counts[:,3] = col4_counts_list[indexstart:indexend]
    col_counts[:,4] = col5_counts_list[indexstart:indexend]
    col_counts[:,5] = col6_counts_list[indexstart:indexend]
    L_shell[:] = L_shell_list[indexstart:indexend]
    MLT[:] = MLT_list[indexstart:indexend]
    return time_new, col_counts, L_shell, MLT, range

#---------- prints timestamps and average calculations to .txt file -----------
def main():
    data = read_data(filename)
    timestamp_list, col_counts, L_shell, MLT, timestep = select_time(data)
    timestamp = []
    for ntime in range(len(timestamp_list)):
        #print(timestamp_list[ntime])
        
        timestamp.append(timestamp_list[ntime].strftime("%Y-%m-%dT%H:%M:%S.%fZ"))
    
    np.savetxt(outfile, np.column_stack((timestamp, col_counts[:,0], col_counts[:,1], col_counts[:,2], col_counts[:,3], col_counts[:,4], col_counts[:,5], L_shell, MLT)), fmt = '%s', delimiter=', ')

    print('Files saved as ', outfile)
    print('      Time, Col_counts[0:5], L-shell, MLT')

    
main()



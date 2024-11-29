# -*- coding: utf-8 -*-
"""
Analyze range form different stares
"""

import os
cd=os.path.dirname(__file__)
import sys
sys.path.append('C:/Users/SLETIZIA/OneDrive - NREL/Desktop/PostDoc/utils')
import utils as utl
import xarray as xr
import numpy as np
import glob
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
plt.close('all')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 14

#%% Inputs
#dataset
source=os.path.join(cd,'data/sms.lidar.z01.a0/*0000*nc')
sdate='2024-08-28'#start date
edate='2024-11-30'#end date
power_cuts=['2024-09-04T23:31:13',
            '2024-10-24T16:54:15',
            '2024-11-05T18:06:36',
            '2024-11-11T21:23:35']

#stats
min_duration=3500#[s] minimum stare duration
min_ele=89#[deg] minimum elevation
min_range=500 #[m] minimum range for selection
max_range=1000#[m] maximum range for selection
max_beta=10**-5#max beta
min_range_beta=1000#[m] minimum range for beta
max_range_beta=2000#[m] maximum range for beta
window=5#time window for MA

#user
r=np.arange(0,3000)#[m] common range


#%% Initialization
files=glob.glob(source)

#zeroing
ppr=[]
dr=[]
snr=[]
rws=[]
beta=[]
rws_std=[]
time=np.array([],dtype='datetime64')
duration=[]

#%% Main
for f in files:
    Data=xr.open_dataset(f)   
    
    if Data.time.values[0]>np.datetime64(sdate+'T00:00:00') and\
       Data.time.values[0]<np.datetime64(edate+'T00:00:00') and\
        np.float(np.ptp(Data.time.values))/10**9>min_duration and\
            np.min(Data['Elevation (degrees)'])>min_ele:
            
        time=np.append(time,Data.time.values[0])
        duration=np.append(duration,np.float(np.ptp(Data.time.values))/10**9)
        ppr=np.append(ppr,Data.attrs['Pulses or ray'])
        dr=np.append(dr,Data.attrs['Range gate length (m)'])
        
        r1=Data.distance.values
        snr=utl.vstack(snr,np.interp(r,r1,Data['SNR'].median(dim='time').values))
        rws=utl.vstack(rws,np.interp(r,r1,Data['Radial_wind_speed'].median(dim='time').values))
        rws_std=utl.vstack(rws_std,np.interp(r,r1,Data['Radial_wind_speed'].std(dim='time').values))
        beta=utl.vstack(beta,np.interp(r,r1,Data.Beta.median(dim='time').values))


Output=xr.Dataset()
Output['snr']=xr.DataArray(data=snr,coords={'time':time,'range':r})
Output['rws']=xr.DataArray(data=rws,coords={'time':time,'range':r})
Output['rws_std']=xr.DataArray(data=rws_std,coords={'time':time,'range':r})
Output['beta']=xr.DataArray(data=beta,coords={'time':time,'range':r})
Output['ppr']=xr.DataArray(data=ppr,coords={'time':time})
Output['dr']=xr.DataArray(data=dr,coords={'time':time})

#exclude precipitation
Output['max_beta']=Output.beta.where(Output.range>=min_range_beta).where(Output.range<max_range_beta).max(dim='range')
Output=Output.where(Output['max_beta']<max_beta)

#resample
Output=Output.resample(time='24h').nearest(tolerance='12h')

#extract selection
snr_sel=Output.snr.where(Output.range>=min_range).where(Output.range<max_range).median(dim='range')
snr_sel_ma=snr_sel.rolling(time=window, center=True).construct("window").median(dim='window')

rws_std_sel=Output.rws_std.where(Output.range>=min_range).where(Output.range<max_range).median(dim='range')
rws_std_sel_ma=rws_std_sel.rolling(time=window, center=True).construct("window").median(dim='window')

#%% Plots
plt.close('all')
fig=plt.figure(figsize=(18,10))
gs = gridspec.GridSpec(3, 2, height_ratios=[5,2,2],width_ratios=[1,0.02])

ax=fig.add_subplot(gs[0,0])
pc=plt.pcolor(Output.time,Output.range,Output.snr.T,cmap='hot',vmin=-25,vmax=0)

plt.ylabel('Range [m]')
plt.grid()
plt.xlim([Output.time[0],Output.time[-1]])

for p in power_cuts:
    plt.plot([np.datetime64(p),np.datetime64(p)],Output.range.values[[0,-1]],'--m',linewidth=3)
    
ax=fig.add_subplot(gs[0,1])
fig.colorbar(pc,cax=ax,label='SNR [dB]')

ax=fig.add_subplot(gs[1,0])
plt.plot(Output.time,Output.ppr,'.k')
plt.ylabel('Pulses/ray')
plt.xlim([Output.time[0],Output.time[-1]])
plt.grid()

ax=fig.add_subplot(gs[2,0])
plt.plot(Output.time,Output.dr,'.k')
plt.ylabel('Gate length [m]')
plt.xlabel('Time (UTC)')
plt.xlim([Output.time[0],Output.time[-1]])
plt.grid()

fig=plt.figure(figsize=(18,10))
gs = gridspec.GridSpec(3, 2, height_ratios=[5,2,2],width_ratios=[1,0.02])

ax=fig.add_subplot(gs[0,0])
pc=plt.pcolor(Output.time,Output.range,np.log10(Output.beta.T),cmap='hot')

plt.ylabel('Range [m]')
plt.grid()
plt.xlim([Output.time[0],Output.time[-1]])

for p in power_cuts:
    plt.plot([np.datetime64(p),np.datetime64(p)],Output.range.values[[0,-1]],'--m',linewidth=3)
    
ax=fig.add_subplot(gs[0,1])
fig.colorbar(pc,cax=ax,label='log(Beta)')


ax=fig.add_subplot(gs[1,0])
plt.plot(Output.time,Output.ppr,'.k')
plt.ylabel('Pulses/ray')
plt.xlim([Output.time[0],Output.time[-1]])
plt.grid()

ax=fig.add_subplot(gs[2,0])
plt.plot(Output.time,Output.dr,'.k')
plt.ylabel('Gate length [m]')
plt.xlabel('Time (UTC)')
plt.xlim([Output.time[0],Output.time[-1]])
plt.grid()



fig=plt.figure(figsize=(18,10))
gs = gridspec.GridSpec(3, 2, height_ratios=[5,2,2],width_ratios=[1,0.02])

ax=fig.add_subplot(gs[0,0])
pc=plt.pcolor(Output.time,Output.range,Output.rws_std.T,cmap='hot')

plt.ylabel('Range [m]')
plt.grid()
plt.xlim([Output.time[0],Output.time[-1]])

for p in power_cuts:
    plt.plot([np.datetime64(p),np.datetime64(p)],Output.range.values[[0,-1]],'--m',linewidth=3)
    
ax=fig.add_subplot(gs[0,1])
fig.colorbar(pc,cax=ax,label=r'StDev of RWS [m s$^{-1}$]')

ax=fig.add_subplot(gs[1,0])
plt.plot(Output.time,Output.ppr,'.k')
plt.ylabel('Pulses/ray')
plt.xlim([Output.time[0],Output.time[-1]])
plt.grid()

ax=fig.add_subplot(gs[2,0])
plt.plot(Output.time,Output.dr,'.k')
plt.ylabel('Gate length [m]')
plt.xlabel('Time (UTC)')
plt.xlim([Output.time[0],Output.time[-1]])
plt.grid()


fig=plt.figure(figsize=(18,10))
gs = gridspec.GridSpec(4, 1, height_ratios=[5,5,2,2])

ax=fig.add_subplot(gs[0,0])
plt.plot(Output.time,snr_sel,'.r',markersize=10)
plt.plot(Output.time, snr_sel_ma,'-r',linewidth=2,label=f'{window}-day mean')
plt.ylabel(f'SNR between {min_range} m \n and {max_range} m')
plt.xlim([Output.time[0],Output.time[-1]])
plt.grid()

for p in power_cuts:
    plt.plot([np.datetime64(p),np.datetime64(p)],[-25,-10],'--m',linewidth=3)


ax=fig.add_subplot(gs[1,0])
plt.plot(Output.time,rws_std_sel,'.r',markersize=10,label='No average')
plt.plot(Output.time, rws_std_sel_ma,'-r',linewidth=2,label=f'{window}-day mean')
plt.ylabel(f'StDev of RWS \n between {min_range} m and {max_range} ' +r'[m s$^{-1}$')
plt.xlim([Output.time[0],Output.time[-1]])
plt.grid()
plt.legend(draggable=True)

for p in power_cuts:
    plt.plot([np.datetime64(p),np.datetime64(p)],[0,25],'--m',linewidth=3)

ax=fig.add_subplot(gs[2,0])
plt.plot(Output.time,Output.ppr,'.k')
plt.ylabel('Pulses/ray')
plt.xlim([Output.time[0],Output.time[-1]])
plt.grid()

ax=fig.add_subplot(gs[3,0])
plt.plot(Output.time,Output.dr,'.k')
plt.ylabel('Gate length [m]')
plt.xlabel('Time (UTC)')
plt.xlim([Output.time[0],Output.time[-1]])
plt.grid()




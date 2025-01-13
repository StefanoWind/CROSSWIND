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
import matplotlib
plt.close('all')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 18

#%% Inputs
source=os.path.join(cd,'data/sms.lidar.z01.a0/*20240904*nc')

#graphics
colors={12:'b',
        21:'g',
        30:'r'}

#%% Initialization
files=glob.glob(source)

#zeroing
ppr=[]
dr=[]
snr=[]
rws=[]
r=[]
rws_std=[]

#%% Main
for f in files:
    Data=xr.open_dataset(f)
    ppr=np.append(ppr,Data.attrs['Pulses or ray'])
    dr=np.append(dr,Data.attrs['Range gate length (m)'])
    snr=utl.vstack(snr,Data['SNR'].mean(dim='time'))
    rws=utl.vstack(rws,Data['Radial_wind_speed'].mean(dim='time'))
    rws_std=utl.vstack(rws_std,Data['Radial_wind_speed'].std(dim='time'))
    r=utl.vstack(r,Data.distance.values)

order=np.argsort(ppr+dr)

#%% Plots
plt.figure(figsize=(16,10))
plt.subplot(2,1,1)
for i in order:
    plt.semilogx(r[i,:],snr[i,:],color=colors[dr[i]],linewidth=ppr[i]/1000)
plt.ylabel(r'SNR [dB]')
plt.xlim([100,2000])
plt.ylim([-30,-10])
plt.grid()
plt.title(source[len(cd):])

plt.subplot(2,1,2)
for i in order:
    plt.loglog(r[i,:],rws_std[i,:],label=r'PPR$='+str(int(ppr[i]))+'$, $\Delta r='+str(int(dr[i]))+'$ m',color=colors[dr[i]],linewidth=ppr[i]/1000)
plt.legend(draggable=True)
plt.grid()
plt.xlabel('Range [m]')
plt.ylabel(r'RWS StDev [m s$^{-1}$]')
plt.xlim([100,2000])
plt.ylim([0.1,25])
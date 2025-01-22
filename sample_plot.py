# -*- coding: utf-8 -*-
"""
Plot sample of stare data
"""

import os
cd=os.path.dirname(__file__)
import xarray as xr
import numpy as np
import glob
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import warnings
plt.close('all')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 16
warnings.filterwarnings('ignore')

#%% Inputs
source=os.path.join(cd,'data/crosswind/nwtc.lidar.z01.a0/*nc')
fileID=0

#%% Initialization
files=glob.glob(source)
Data=xr.open_dataset(files[fileID])

#%% Plots
plt.figure(figsize=(18,10))
plt.pcolor(Data.time,Data.distance,Data.radial_wind_speed.T,vmin=-5,vmax=5,cmap='seismic')
plt.xlabel('Time (UTC)')
plt.ylabel('Range [m]') 
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.colorbar(label=r'Radial wind speed [m s$^{-1}$]')


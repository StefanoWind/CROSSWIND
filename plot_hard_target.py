# -*- coding: utf-8 -*-
"""
Plot hard target results
"""

import os
cd=os.path.dirname(__file__)
import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 13
matplotlib.rcParams['savefig.dpi']=300
plt.close("all")

#%% Inputs
source=os.path.join(cd,'data','crosswind','sgre.lidar.z01.a0','sgre.lidar.z01.a0.20251111.225629.user1.nc')
min_SNR=10#[dB] minimum SNR for hard taget
max_SNR_bloc=-22.5#[dB] maximum mean near range SNR for blocked beams
max_range_bloc=500#[m] near range averaging distance for blecked beams
min_r=60#[s] minimum range
max_r=300#[m] range so split far and near range
max_ele=45#deg] maximum elevation
scan_openings=[20]#[deg] scan half-opening angles

#%% Initialization
Data=xr.open_dataset(source)

#%% Main
Data['SNR']=Data.SNR.where(~np.isnan(Data.SNR),Data.SNR.min())#fill null SNR

#Cartesian coordinates
Data['x']=Data.distance*np.cos(np.radians(90-Data.azimuth))*np.cos(np.radians(Data.elevation))
Data['y']=Data.distance*np.sin(np.radians(90-Data.azimuth))*np.cos(np.radians(Data.elevation))
Data['z']=Data.distance*np.sin(np.radians(Data.elevation))

#drop exlcuded area
Data=Data.where((Data.distance>=min_r)*(Data.elevation<=max_ele),drop=True)

#fix coords
Data['azimuth']=Data.azimuth.isel(range_gate=0)
Data['elevation']=Data.elevation.isel(range_gate=0)
Data['distance']=Data.distance.isel(time=0)

#hard-target flag
Data['ht']=Data.SNR>min_SNR
blocked=Data.SNR.where(Data.distance<=max_range_bloc).mean(dim='range_gate')<max_SNR_bloc
Data['ht'].isel(range_gate=0)[:] = xr.where(blocked, True, Data.ht.isel(range_gate=0))
Data_ht_near=Data.where(Data.ht*(Data.distance<max_r),drop=True)
Data_ht_far=Data.where(Data.ht*(Data.distance>max_r),drop=True)

#%% Main

fig=plt.figure(figsize=(18,10))
ax = fig.add_subplot(projection='3d')
ax.scatter(Data_ht_near.x,Data_ht_near.y,Data_ht_near.z,s=3,c='k',alpha=0.5)
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$y$ [m]')
ax.set_zlabel(r'$z$ [m]')
ax.set_aspect('equal')

fig=plt.figure(figsize=(18,10))
ax = fig.add_subplot(projection='3d')
x=Data_ht_far.x.values
y=Data_ht_far.y.values
z=Data_ht_far.z.values
sc=ax.scatter(x[~np.isnan(x+y+z)],y[~np.isnan(x+y+z)],z[~np.isnan(x+y+z)],s=2,c=z[~np.isnan(x+y+z)],cmap='summer',vmin=-200,vmax=50)
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$y$ [m]')
ax.set_zlabel(r'$z$ [m]')
ax.set_aspect('equal')
plt.colorbar(sc,label='$z$ [m]')

plt.figure(figsize=(18,6))
plt.scatter(Data.azimuth,Data.elevation,
            s=10,c=Data.ht.isel(range_gate=0))
plt.xlabel(r'$\alpha$ [$^\circ$]')
plt.ylabel(r'$\beta$ [$^\circ$]')
plt.grid()
plt.xticks(np.arange(0,361,10))
plt.yticks(np.arange(-20,31,5))
plt.plot(90,0,'xr')
for so in scan_openings:
    plt.plot([-so+90,-so+90,so+90,so+90,-so+90],
             [-so,so,so,-so,-so],'r')
plt.gca().set_aspect('equal')
plt.xlim([0,360])
plt.ylim([-20,30])

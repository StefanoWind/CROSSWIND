# -*- coding: utf-8 -*-
"""
Hard target test
"""
import os
cd=os.path.dirname(__file__)
import lidargo as lg
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import xarray as xr
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.size'] = 12
mpl.rcParams['savefig.dpi']=300

plt.close('all')

#%% Inputs

source=os.path.join(cd,'data/crosswind/sgre.lidar.z01.00/User1_217_20251111_225629.hpl')
source_a0=os.path.join(cd,'data/crosswind/sgre.lidar.z01.a0/sgre.lidar.z01.a0.20251111.225629.user1.nc')
config_format={'model':'halo','site':'sgre','instrument_id':1,'data_level_out':'a0'}
#site
D=101#[m] rotor diameter
H=80#[m] hub height

#%% Initialization


#%% Main
lproc = lg.Format(source, config=config_format, verbose=True)
lproc.process_scan(replace=False, save_file=True,save_path=os.path.join(cd,'data/crosswind/sgre.lidar.z01.a0'))

Data=xr.open_dataset(source_a0)
Data['x']=Data.distance*np.cos(np.radians(90-Data.azimuth))
Data['y']=Data.distance*np.sin(np.radians(90-Data.azimuth))
Data['z']=Data.distance*np.sin(np.radians(Data.elevation))

#%% Plots
fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
ax.scatter(Data.x,Data.y,Data.z,s=1,c=Data.SNR)
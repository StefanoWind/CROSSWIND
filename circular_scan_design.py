# -*- coding: utf-8 -*-
"""
Design a ZX-like circular scan
"""
import os
cd=os.path.dirname(__file__)
from halo_suite.utilities import scan_file_compiler
import numpy as np
from matplotlib import pyplot as plt
import yaml
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.size'] = 12
mpl.rcParams['savefig.dpi']=300

plt.close('all')

#%% Inputs

#site
gamma=3#[deg] half-opening angle
dtheta=8#[deg]
r=300#[m]
path_config='C:/Users/sletizia/Software/FIEXTA/halo_suite/halo_suite/configs/config.217.yaml'

#%% Functions
def cos(x):
    return np.cos(np.radians(x))
def sin(x):
    return np.sin(np.radians(x))
def tan(x):
    return np.tan(np.radians(x))
def atan(x):
    return np.degrees(np.atan(x))

#%% Initialization
theta=np.arange(0,360,dtheta)
with open(path_config, 'r') as fid:
    config_lidar = yaml.safe_load(fid)  
    
#%% Main
alpha=atan(tan(gamma)*cos(theta))
ele=  atan(tan(theta)*sin(alpha))
ele[theta==90]=gamma
ele[theta==270]=-gamma

azi=(90-alpha)%360
x=r*cos(ele)*cos(90-azi)
y=r*cos(ele)*sin(90-azi)
z=r*sin(ele)

scan_file=scan_file_compiler('SSM',azi=azi,ele=ele,identifier=f'circular_{gamma}deg',
                             save_path=os.path.join(cd,'scans'),reset=True)

#%% Plots
fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
plt.plot(x,y,z,'.k')
plt.plot(0,0,0,'.r')
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$y$ [m]')
ax.set_zlabel(r'$z$ [m]')
ax.set_aspect('equal')

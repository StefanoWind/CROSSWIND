# -*- coding: utf-8 -*-
"""
Optimize scan
"""
import os
cd=os.path.dirname(__file__)
from lisboa import scan_optimizer as opt
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.size'] = 12
mpl.rcParams['savefig.dpi']=300

plt.close('all')

#%% Inputs

#user
full_scan_file=True

#site
D=108#[m] rotor diameter
H=80#[m] hub height

#lidar
azi_offset=0#[deg] difference between scan direction and x axis
ppr=1000#pulses per ray
dr=30#[m] gate length
rmin=100#minimum range
rmax=1000#maximum range

#Pareto inputs
scan_widths=np.array([15])# width of the scan
# ang_resolutions=np.array([1,5/4,5/3,5/2,5])#[deg] angular resolutions
num_ang=np.array([17])#number of beams
xmin=1#[D] min x
xmax=6#[D] max x
ymin=-1.5#[D] min y
ymax=1.5#[D] max y
zmax=1.5#[D] x max
Dn0=[1,0.25,0.25]#[D] fundamental half_wavelength
T=600#[s] scan duration
tau=4#[s] timescale in the wake

#%% Initialization
azi1=90-scan_widths
azi2=90+scan_widths
ele1=-scan_widths
ele2=+scan_widths
dazi=None
dele=None
num_azi=num_ang
num_ele=num_ang

coords='xyz'
path_config_lidar='C:/Users/sletizia/Software/FIEXTA/halo_suite/halo_suite/configs/config.217.yaml'
volumetric=True
mode='CSM'

#LiSBOA config
config={'sigma':0.25,
        'max_iter':5,
        'mins':[D*xmin,D*ymin,-H],
        'maxs':[D*xmax,D*ymax,D*zmax],
        'Dn0':np.array(Dn0)*D,
        'r_max':3,
        'dist_edge':1,
        'tol_dist':0.1,
        'grid_factor':0.25,
        'max_Dd':1}

#%% Main
scopt=opt.scan_optimizer(config,save_path=os.path.join(cd,'data','Pareto'),logfile=os.path.join(cd,'log','test.log'))

Pareto=scopt.pareto(coords,azi1,azi2,ele1,ele2,dazi,dele,num_azi,num_ele,
                    volumetric=volumetric,rmin=rmin,rmax=rmax, T=T,tau=tau,
                    mode=mode, ppr=ppr, dr=dr, path_config_lidar=path_config_lidar,
                    full_scan_file=full_scan_file)
# -*- coding: utf-8 -*-
"""
Extract performance indicators from stare data
"""

import os
cd=os.path.dirname(__file__)
import xarray as xr
import numpy as np
import glob
from matplotlib import pyplot as plt
import matplotlib
import warnings
plt.close('all')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 12
warnings.filterwarnings('ignore')

#%% Inputs
sources={'257':os.path.join(cd,'data/crosswind/nwtc.lidar.z01.a0/*nc'),
         '178':os.path.join(cd,'data/crosswind/nwtc.lidar.z02.a0/*nc')}

r=np.arange(0,3000,10)#[m] common range

#%% Functions
def vstack(a,b):
    '''
    Stack vertically vectors
    '''
    if len(a)>0:
        ab=np.vstack((a,b))
    else:
        ab=b
    return ab  

#%% Initialization

#%% Main
for s in sources:
    ppr=[]
    dr=[]
    snr=[]
    rws=[]
    beta=[]
    rws_std=[]
    time=np.array([],dtype='datetime64')
    dt=[]
    files=glob.glob(sources[s])
    for f in files:
        Data=xr.open_dataset(f)
        
        time0=Data.time.values[0]
        time=np.append(time,np.round(time0.astype('datetime64[m]').astype(float) / 60).astype('datetime64[h]'))
        dt=np.append(dt,np.float64(np.nanmedian(np.diff(Data.time.values)))/10**9)
        ppr=np.append(ppr,Data.attrs['Pulses per ray'])
        dr=np.append(dr,Data.attrs['Range gate length (m)'])
        
        r1=Data.distance.values
        snr=vstack(snr,np.interp(r,r1,Data['SNR'].median(dim='time').values))
        rws=vstack(rws,np.interp(r,r1,Data['Radial_wind_speed'].median(dim='time').values))
        rws_std=vstack(rws_std,np.interp(r,r1,Data['Radial_wind_speed'].std(dim='time').values))
        beta=vstack(beta,np.interp(r,r1,Data.Beta.median(dim='time').values))


    Output=xr.Dataset()
    Output['snr']=xr.DataArray(data=snr,coords={'time':time,'range':r})
    Output['rws']=xr.DataArray(data=rws,coords={'time':time,'range':r})
    Output['rws_std']=xr.DataArray(data=rws_std,coords={'time':time,'range':r})
    Output['beta']=xr.DataArray(data=beta,coords={'time':time,'range':r})
    Output['ppr']=xr.DataArray(data=ppr,coords={'time':time})
    Output['dr']=xr.DataArray(data=dr,coords={'time':time})
    Output['dt']=xr.DataArray(data=dt,coords={'time':time})
    
    Output.to_netcdf(os.path.join(cd,'data',os.path.dirname(f).split('/')[-1]+'.all.nc'))
    
    
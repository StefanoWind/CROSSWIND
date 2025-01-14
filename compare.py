# -*- coding: utf-8 -*-
"""
Compare performance of two lidars
"""

import os
cd=os.path.dirname(__file__)
import sys
sys.path.append('C:/Users/SLETIZIA/OneDrive - NREL/Desktop/Main/utils')
import utils as utl
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
plt.close('all')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] =26

#%% Inputs
sources={'257':os.path.join(cd,'data/nwtc.lidar.z01.a0.all.nc'),
         '178':os.path.join(cd,'data/nwtc.lidar.z02.a0.all.nc')}

colors={'257':'r','178':'k'}

#%% Initialization
Data={}
pprs=[]
drs=[]
for s in sources:
   Data[s]=xr.open_dataset(sources[s]).sortby('time')
   Data[s]=Data[s]
   pprs=np.append(pprs,Data[s].ppr)
   drs=np.append(drs,Data[s].dr)

pprs=[np.int(x) for x in np.unique(pprs)]
drs=[np.int(x) for x in np.unique(drs)]

# Select the common times for each dataset
common_times = xr.DataArray(sorted(set.intersection(*(set(Data[s].time.values) for s in sources))),dims="time")
for s in sources:
    Data[s] = Data[s].sel(time=common_times) 
    
#data availability
da=xr.DataArray(np.zeros((len(drs),len(pprs))),coords={'dr':drs,'ppr':pprs})
for ppr in pprs:
    for dr in drs:
        da.loc[dr,ppr]=np.sum(~np.isnan(Data[s]['dt'].where(Data[s].ppr==ppr).where(Data[s].dr==dr))).values
        

#%% Main
plt.close('all')

#data availability
fig=plt.figure(figsize=(12,10))
plt.pcolor(da.values,cmap='RdYlGn')
    
ippr=0
for ppr in pprs:
    idr=0
    for dr in drs: 
        plt.text(ippr+0.5-0.05,idr+0.5,f'{str(int(da.loc[dr,ppr].values))}',fontsize=26,fontweight='bold')
        idr+=1
    ippr+=1
plt.xticks(np.arange(len(pprs))+0.5,pprs)
plt.yticks(np.arange(len(drs))+0.5,drs)
plt.xlabel('Pulses/ray')
plt.ylabel('Gate length [m]')
plt.colorbar(label='Data availability [hrs]')

#SNR
fig=plt.figure(figsize=(30,20))
ippr=0
for ppr in pprs:
    idr=0
    for dr in drs:
        ax=plt.subplot(len(drs),len(pprs),(len(drs)-idr-1)*len(pprs)+ippr+1)
        for s in sources:
            Data_sel=Data[s].where(Data[s].ppr==ppr).where(Data[s].dr==dr)
            plt.plot(Data[s].range,Data_sel['snr'].mean(dim='time'),'-',
                     color=colors[s],label=s)
            ax.fill_between(Data[s].range,Data_sel['snr'].mean(dim='time')+Data_sel['snr'].std(dim='time'),
                            Data_sel['snr'].mean(dim='time')-Data_sel['snr'].std(dim='time'),color=colors[s],alpha=0.1)
            
        plt.text(2500*0.78,-5-1,f'PPR={str(int(ppr))}\nΔr={str(int(dr))} m', bbox={'edgecolor':'k','facecolor':'b','alpha':0.25})
        if idr==0:
            plt.xlabel('Range [m]')
        else:
            ax.set_xticklabels([])
        if ippr==0:
            plt.ylabel(r'SNR [dB]')
        else:
            ax.set_yticklabels([])
        plt.grid()
        plt.ylim([-30,0])
                       
        idr+=1
    ippr+=1
plt.legend(draggable=True,loc='upper left')
plt.tight_layout()


#RWS std
fig=plt.figure(figsize=(30,20))
ippr=0
for ppr in pprs:
    idr=0
    for dr in drs:
        ax=plt.subplot(len(drs),len(pprs),(len(drs)-idr-1)*len(pprs)+ippr+1)
        for s in sources:
            Data_sel=Data[s].where(Data[s].ppr==ppr).where(Data[s].dr==dr)
            plt.plot(Data[s].range,Data_sel['rws_std'].mean(dim='time'),'-',
                     color=colors[s],label=s)
            ax.fill_between(Data[s].range,Data_sel['rws_std'].mean(dim='time')+Data_sel['rws_std'].std(dim='time'),
                            Data_sel['rws_std'].mean(dim='time')-Data_sel['rws_std'].std(dim='time'),color=colors[s],alpha=0.1)
            
        plt.text(2500*0.78,5,f'PPR={str(int(ppr))},\nΔr={str(int(dr))} m', bbox={'edgecolor':'k','facecolor':'b','alpha':0.25})
        if idr==0:
            plt.xlabel('Range [m]')
        else:
            ax.set_xticklabels([])
        if ippr==0:
            plt.ylabel(r'StDev of radial wind speed [m s$^{-1}$]')
        else:
            ax.set_yticklabels([])
        plt.grid()
        plt.ylim([0,30])
                       
        idr+=1
    ippr+=1
plt.legend(draggable=True,loc='upper left')
plt.tight_layout()

#RWS mean
fig=plt.figure(figsize=(30,20))
ippr=0
for ppr in pprs:
    idr=0
    for dr in drs:
        ax=plt.subplot(len(drs),len(pprs),(len(drs)-idr-1)*len(pprs)+ippr+1)
        for s in sources:
            Data_sel=Data[s].where(Data[s].ppr==ppr).where(Data[s].dr==dr)
            plt.plot(Data[s].range,Data_sel['rws'].mean(dim='time'),'-',
                     color=colors[s],label=s)
            ax.fill_between(Data[s].range,Data_sel['rws'].mean(dim='time')+Data_sel['rws'].std(dim='time'),
                            Data_sel['rws'].mean(dim='time')-Data_sel['rws'].std(dim='time'),color=colors[s],alpha=0.1)
            
        plt.text(2500*0.78,-2,f'PPR={str(int(ppr))}\nΔr={str(int(dr))} m', bbox={'edgecolor':'k','facecolor':'b','alpha':0.25})
        if idr==0:
            plt.xlabel('Range [m]')
        else:
            ax.set_xticklabels([])
        if ippr==0:
            plt.ylabel(r'Mean of radial wind speed [m s$^{-1}$]')
        else:
            ax.set_yticklabels([])
        plt.grid()
        plt.ylim([-5,5])
                       
        idr+=1
    ippr+=1
plt.legend(draggable=True,loc='upper left')
plt.tight_layout()

#sampling time
fig=plt.figure(figsize=(20,17))
ippr=0
for ppr in pprs:
    idr=0
    for dr in drs:
        ax=plt.subplot(len(drs),len(pprs),(len(drs)-idr-1)*len(pprs)+ippr+1)
       
        plt.bar(list(sources.keys()),[Data[s].where(Data[s].ppr==ppr).where(Data[s].dr==dr).dt.mean() for s in sources],color='k')
        
        plt.text(0,1.15,f'PPR={str(int(ppr))}\nΔr={str(int(dr))} m', bbox={'edgecolor':'k','facecolor':'b','alpha':0.25})
        if idr==0:
            plt.xlabel('Lidar SN')
        else:
            ax.set_xticklabels([])
        if ippr==0:
            plt.ylabel(r'Sampling time [s]')
        else:
            ax.set_yticklabels([])
        plt.grid()
        plt.ylim([0,1.5])
                       
        idr+=1
    ippr+=1
plt.tight_layout()   
    

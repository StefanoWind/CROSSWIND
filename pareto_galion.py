# -*- coding: utf-8 -*-
"""
Optimize Galion scan with LiSBOA
"""
import numpy as np
from lisboa import statistics as stats
from matplotlib import pyplot as plt
from matplotlib import cm
import warnings
import matplotlib

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 12

warnings.filterwarnings('ignore')
plt.close('all')

#%% Inputs

#site
Nr=34#number of gates
dr=30#[m] gate length
dt=3.47#[s] sampling time assumed constant
D=108#[m] rotor diameter
H=80#[m] hub height
St=0.4#Strouhal number in the turbine wake (Letizia et al. 2021)
U_inf=8#[m/s] nominal wind speed (for timescale only)

#Pareto
T_tot=1800#[s] experiment duration

dangs_2d=[0.25,0.5,1,2]#[deg] tested angular resolutions (PPI)
max_azi_2d=[15,20,25,30] #[deg] tested widths (PPI)

dangs_3d=[1,2,3]#[deg] tested angular resolutions (volumetric)
max_angs_3d=[[10,10],[15,15],[20,20]] #[deg] tested widths (volumetric) [azi,ele]

#LiSBOA config (PPI)
config_2d={'sigma':0.25,
        'mins':[D,-1.5*D],
        'maxs':[8*D,1.5*D],
        'Dn0':[D,0.25*D],
        'r_max':3,
        'tol_dist':0.1,
        'grid_factor':0.25,
        'max_Dd':1}

#LiSBOA config (volumetric)
config_3d={'sigma':0.25,
        'mins':[D,-1.5*D,-H],
        'maxs':[8*D,1.5*D,1.5*D],
        'Dn0':[D,0.25*D,0.25*D],
        'r_max':3,
        'tol_dist':0.1,
        'grid_factor':0.25,
        'max_Dd':1}

#% graphics
markers=['o','v','s','p','*']
cmap = cm.get_cmap('viridis')

#%% Initialization
r=np.arange(Nr)*dr+dr/2#ranges

lproc=stats.statistics()#load lisboa

tau=St*D/U_inf#[s] integral timescale

#%%Main

#PPI optimization
epsilon1_2d=np.zeros((len(max_azi_2d),len(dangs_2d)))
epsilon2_2d=np.zeros((len(max_azi_2d),len(dangs_2d)))
fig= plt.figure(figsize=(18,8))
i_max_ang=0
for max_ang in max_azi_2d:
    
    i_dang=0
    for dang in dangs_2d:
        
        #geometry
        azi=np.unique(np.append(-np.arange(0,max_ang+dang/2,dang)[::-1],np.arange(0,max_ang+dang/2,dang)))+90
        ele=azi*0
        
        #sampling points
        x=np.outer(np.cos(np.radians(ele))*np.cos(np.radians(90-azi)),r)
        y=np.outer(np.cos(np.radians(ele))*np.sin(np.radians(90-azi)),r)
        z=np.outer(np.sin(np.radians(ele)),r)
    
        #LiSBOA
        x_exp=[x.ravel(),y.ravel()]
        X2,Dd,_,_,_=lproc.calculate_weights(config_2d, x_exp)
        
        #first cost function
        excl=Dd>config_2d['max_Dd']
        epsilon1_2d[i_max_ang,i_dang]=np.sum(excl)/np.size(Dd)*100
        
        #second cost function
        T=len(azi)*dt
        L=np.floor(T_tot/T)
        
        p=np.arange(1,L)
        epsilon2_2d[i_max_ang,i_dang]=(1/L+2/L**2*np.sum((L-p)*np.exp(-T/tau*p)))**0.5*100
    
        #plots
        ax=plt.subplot(len(max_azi_2d),len(dangs_2d),i_max_ang*len(dangs_2d)+i_dang+1)
        plt.pcolor(X2[0]/D,X2[1]/D,Dd,vmin=0,vmax=config_2d['max_Dd'],cmap='RdYlGn_r')
        plt.title(r'$\theta_{max} ='+str(max_azi_2d[i_max_ang])+r'^\circ$, $\Delta \theta='+str(dangs_2d[i_dang])+r'^\circ$, $\tau_s= '+str(int(T))+r'$ s',fontsize=12)
        plt.plot(x/D,y/D,'.k',markersize=1)
        plt.grid()
        plt.gca().set_aspect(1)
        plt.xlim([config_2d['mins'][0]/D,config_2d['maxs'][0]/D])
        plt.ylim([config_2d['mins'][1]/D,config_2d['maxs'][1]/D])
        if i_max_ang<len(max_azi_2d)-1:
            ax.set_xticklabels([])
        else:
            plt.xlabel(r'$x/D$')
        if i_dang>0:
            ax.set_yticklabels([])
        else:
            plt.ylabel(r'$y/D$')
        i_dang+=1
        print(i_dang)
        
    i_max_ang+=1
# raise BaseException()
    
#3D scans
epsilon1_3d=np.zeros((len(max_angs_3d),len(dangs_3d)))
epsilon2_3d=np.zeros((len(max_angs_3d),len(dangs_3d)))
fig= plt.figure(figsize=(18,8))
i_max_ang=0
for max_ang in max_angs_3d:
    
    i_dang=0
    for dang in dangs_3d:
        
        #geometry
        azi=np.unique(np.append(-np.arange(0,max_ang[0]+dang/2,dang)[::-1],np.arange(0,max_ang[0]+dang/2,dang)))+90
        ele=np.unique(np.append(-np.arange(0,max_ang[1]+dang/2,dang)[::-1],np.arange(0,max_ang[1]+dang/2,dang)))
        azi2,ele2=np.meshgrid(azi,ele)
        
        #sampling points
        x=np.outer(np.cos(np.radians(ele2.ravel()))*np.cos(np.radians(90-azi2.ravel())),r)
        y=np.outer(np.cos(np.radians(ele2.ravel()))*np.sin(np.radians(90-azi2.ravel())),r)
        z=np.outer(np.sin(np.radians(ele2.ravel())),r)
    
        #LiSBOA
        x_exp=[x.ravel(),y.ravel(),z.ravel()]
        X2,Dd,_,_,_=lproc.calculate_weights(config_3d, x_exp)
        
        #first cost function
        excl=Dd>config_3d['max_Dd']
        epsilon1_3d[i_max_ang,i_dang]=np.sum(excl)/np.size(Dd)*100
        
        #second cost function
        T=len(azi2.ravel())*dt
        L=np.floor(T_tot/T)
        
        p=np.arange(1,L)
        epsilon2_3d[i_max_ang,i_dang]=(1/L+2/L**2*np.sum((L-p)*np.exp(-T/tau*p)))**0.5*100
    
        #plots
        Dd[Dd>config_3d['max_Dd']]=np.nan
        ax=fig.add_subplot(len(max_angs_3d),len(dangs_3d),i_max_ang*len(dangs_3d)+i_dang+1,projection='3d',alpha=0.25)
        ax.scatter(X2[0]/D,X2[1]/D,X2[2]/D,s=1,c=Dd,vmax=config_2d['max_Dd'],cmap='RdYlGn_r')
        plt.title(r'$\theta_{max} ='+str(max_angs_3d[i_max_ang][0])+r'^\circ$, $\beta_{max} ='+str(max_angs_3d[i_max_ang][1])+r'^\circ$, $\Delta \theta = \Delta \beta='+str(dangs_3d[i_dang])+r'^\circ$, $\tau_s= '+str(int(T))+r'$ s',fontsize=12)
        plt.plot(x/D,y/D,z/D,'.k',markersize=1,alpha=0.5)
        plt.grid()
        plt.gca().set_aspect('equal')
        ax.set_xlim([config_3d['mins'][0]/D,config_3d['maxs'][0]/D])
        ax.set_ylim([config_3d['mins'][1]/D,config_3d['maxs'][1]/D])
        ax.set_zlim([config_3d['mins'][2]/D,config_3d['maxs'][2]/D])
        
        if i_max_ang<len(max_angs_3d)-1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel(r'$x/D$')
        if i_dang<len(dangs_3d)-1:
            ax.set_yticklabels([])
            ax.set_zticklabels([])
        else:
            ax.set_ylabel(r'$y/D$')
            ax.set_zlabel(r'$z/D$')
            
        i_dang+=1
        print(i_dang)
        
        # raise BaseException()
        
    i_max_ang+=1
    
    
#%% Plots

#pareto 2D
plt.figure()
for i_max_ang in range(len(max_azi_2d)):
    for i_dang in range(len(dangs_2d)):
        plt.plot(epsilon1_2d[i_max_ang,i_dang],epsilon2_2d[i_max_ang,i_dang],marker=markers[i_max_ang],
             color=cmap(i_dang/(len(dangs_2d)-1)))
        
#legend
for i_max_ang in range(len(max_azi_2d)):
    plt.plot(1000,1000,'.',markeredgecolor='k',marker=markers[i_max_ang],label=r'$\theta_{max}='+str(max_azi_2d[i_max_ang])+'^\circ$',color='k')

for i_dang in range(len(dangs_2d)):
    plt.plot(1000,1000,'.',marker='o',markeredgecolor='k',color=cmap(i_dang/(len(dangs_2d)-1)),label=r'$\Delta \theta = '+str(dangs_2d[i_dang])+'^\circ$')

plt.legend(draggable=True)
plt.xlim([0,100])
plt.ylim([0,100])
plt.grid()
plt.xlabel(r'$\epsilon_I$ [%]')
plt.ylabel(r'$\epsilon_{II}$ [%]')


#pareto 3D
plt.figure()
for i_max_ang in range(len(max_angs_3d)):
    for i_dang in range(len(dangs_3d)):
        plt.plot(epsilon1_3d[i_max_ang,i_dang],epsilon2_3d[i_max_ang,i_dang],marker=markers[i_max_ang],
             color=cmap(i_dang/(len(dangs_3d)-1)))
        
#legend
for i_max_ang in range(len(max_angs_3d)):
    plt.plot(1000,1000,'.',markeredgecolor='k',marker=markers[i_max_ang],
             label=r'$\theta_{max} ='+str(max_angs_3d[i_max_ang][0])+r'^\circ$, $\beta_{max} ='+str(max_angs_3d[i_max_ang][1])+'^\circ$',color='k')

for i_dang in range(len(dangs_3d)):
    plt.plot(1000,1000,'.',marker='o',markeredgecolor='k',color=cmap(i_dang/(len(dangs_3d)-1)),
             label=r'$\Delta \theta =\Delta \beta ='+str(dangs_3d[i_dang])+'^\circ$')

plt.legend(draggable=True)
plt.xlim([0,100])
plt.ylim([0,100])
plt.grid()
plt.xlabel(r'$\epsilon_I$ [%]')
plt.ylabel(r'$\epsilon_{II}$ [%]')
        



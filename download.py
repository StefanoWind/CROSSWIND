# -*- coding: utf-8 -*-
"""
Download and format raw lidar data
"""

import os
cd=os.path.dirname(__file__)
import sys
from matplotlib import pyplot as plt
from datetime import datetime
from datetime import timedelta
import matplotlib
import yaml
import warnings
warnings.filterwarnings('ignore')
plt.close('all')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 14

#%% Inputs

path_config=os.path.join(cd,'configs','config.yaml')

#dataset
channels=['crosswind/nwtc.lidar.z01.a0','crosswind/nwtc.lidar.z02.a0']
t_start='2024-12-17'#start date
t_end='2025-01-15'#end date
file_types=['stare']

#processing
time_increment=24 # [hours] chunk of data downloaded and or processed at every cycle
download=True

#%% Initialization

#configs
with open(path_config, 'r') as fid:
    config = yaml.safe_load(fid)
    
#imports
sys.path.append(config['path_dap_py'])
from doe_dap_dl import DAP

#DAP setup
if download:
    a2e = DAP('a2e.energy.gov',confirm_downloads=False)
    a2e.setup_two_factor_auth(username=config['username'], password=config['password'])

N_periods=(datetime.strptime(t_end, '%Y-%m-%d')-datetime.strptime(t_start, '%Y-%m-%d'))/timedelta(hours=time_increment)
time_bin=[datetime.strptime(t_start, '%Y-%m-%d') + timedelta(hours=time_increment*x) for x in range(int(N_periods)+1)]

#%% Main
for t1,t2 in zip(time_bin[:-1],time_bin[1:]):

    for channel in channels:
        save_path=os.path.join(cd,'data',channel.replace('b0','b1'))
        
        if download:
            _filter = {
                'Dataset': channel,
                'date_time': {
                    'between':  [datetime.strftime(t1, '%Y%m%d%H%M%S'),
                                 datetime.strftime(t2, '%Y%m%d%H%M%S')]
                },
                'file_type': 'nc',
                'ext1':file_types, 
            }
            
            #check if MFA was already in place, otherwise ask to authenticate
            a2e.download_with_order(_filter, path=os.path.join(cd,'data',channel), replace=False)
            
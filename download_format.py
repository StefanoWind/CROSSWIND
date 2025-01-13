# -*- coding: utf-8 -*-
"""
Download and format raw lidar data
"""

import os
cd=os.path.dirname(__file__)
import sys
# import utils as utl
import xarray as xr
import numpy as np
import glob
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
channels=['crosswind/nwtc.lidar.z01.00','crosswind/nwtc.lidar.z02.00']
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

#%% Functions
def read(filename):
    with open(filename, "r") as f:
        lines = []
        for line_num in range(11):
            lines.append(f.readline())

        # read metadata into strings
        metadata = {}
        for line in lines:
            metaline = line.split(":")
            if "Start time" in metaline:
                metadata["Start time"] = metaline[1:]
            else:
                metadata[metaline[0]] = metaline[1]  # type: ignore

        # convert some metadata
        num_gates = int(metadata["Number of gates"])  # type: ignore

        # Read some of the label lines
        for line_num in range(6):
            f.readline()

        # initialize arrays
        time = []
        azimuth = []
        elevation = []
        pitch = []
        roll = []
        doppler = []
        intensity = []
        beta = []

        while True:
            a = f.readline().split()
            if not len(a):  # is empty
                break

            time.append(float(a[0]))
            azimuth.append(float(a[1]))
            elevation.append(float(a[2]))
            pitch.append(float(a[3]))
            roll.append(float(a[4]))

            doppler.append([0] * num_gates)
            intensity.append([0] * num_gates)
            beta.append([0] * num_gates)

            for _ in range(num_gates):
                b = f.readline().split()
                range_gate = int(b[0])
                doppler[-1:][0][range_gate] = float(b[1])
                intensity[-1:][0][range_gate] = float(b[2])
                beta[-1:][0][range_gate] = float(b[3])

    # convert date to np.datetime64
    start_time_string = "{}-{}-{}T{}:{}:{}".format(
        metadata["Start time"][0][1:5],  # year
        metadata["Start time"][0][5:7],  # month
        metadata["Start time"][0][7:9],  # day
        "00",  # hour
        "00",  # minute
        "00.00",  # second
    )

    # find times where it wraps from 24 -> 0, add 24 to all indices after
    new_day_indices = np.where(np.diff(time) < -23)
    for new_day_index in new_day_indices[0]:
        time[new_day_index + 1 :] += 24.0

    start_time = np.datetime64(start_time_string)
    datetimes = [
        start_time + np.timedelta64(int(3600 * 1e6 * dtime), "us") for dtime in time
    ]

    dataset = xr.Dataset(
        {
            "Decimal time (hours)": (("time"), time),
            "Azimuth (degrees)": (("time"), azimuth),
            "Elevation (degrees)": (("time"), elevation),
            "Pitch (degrees)": (("time"), pitch),
            "Roll (degrees)": (("time"), roll),
            "Radial_wind_speed": (("time", "range_gate"), doppler),
            "Intensity": (("time", "range_gate"), intensity),
            "Beta": (("time", "range_gate"), beta),
        },
        coords={"time": np.array(datetimes), "range_gate": np.arange(num_gates)},
        attrs={"Range gate length (m)": float(metadata["Range gate length (m)"])},  # type: ignore
    )

    # Save some attributes
    dataset.attrs["Range gate length (m)"] = float(
        metadata["Range gate length (m)"]  # type: ignore
    )
    dataset.attrs["Number of gates"] = float(metadata["Number of gates"])  # type: ignore
    dataset.attrs["Scan type"] = str(metadata["Scan type"]).strip()
    dataset.attrs["Pulses per ray"] = float(metadata["Pulses/ray"])  # type: ignore
    dataset.attrs["System ID"] = int(metadata["System ID"])  # type: ignore
    dataset.attrs["Filename"] = str(metadata["Filename"])[1:-5]
    dataset["distance"] = (
        "range_gate",
        dataset.coords["range_gate"].data * dataset.attrs["Range gate length (m)"]
        + dataset.attrs["Range gate length (m)"] / 2,
    )
    dataset["distance_overlapped"] = (
        "range_gate",
        dataset.coords["range_gate"].data * 1.5
        + dataset.attrs["Range gate length (m)"] / 2,
    )
    intensity = dataset.Intensity.data.copy()
    intensity[intensity <= 1] = np.nan
    dataset['SNR']=xr.DataArray(data=10 * np.log10(intensity - 1),coords={"time": np.array(datetimes), "range_gate": np.arange(num_gates)})


    # Dynamically add scan type and z-id (z02, z03, etc) to dataset metadata
    # loc_id, instrument, z02/z03, data level, date, time, scan type, extension
    raw_basename = filename.replace("\\", "/").rsplit("/")[-1]
    if ".z" in raw_basename:
        _, _, z_id, _, _, _, scan_type, _ = raw_basename.lower().split(".")
    else:  # local NREL tsdat-ing
        z_id = str(dataset.attrs["System ID"])
        scan_type = ""
        if "user" in raw_basename.lower():
            scan_type = raw_basename.split("_")[0].lower()
        elif "stare" in raw_basename.lower():
            scan_type = "stare"
        elif "vad" in raw_basename.lower():
            scan_type = "vad"
        elif "wind_profile" in raw_basename.lower():
            scan_type = "wind_profile"
        elif "rhi" in raw_basename.lower():
            scan_type = "rhi"

    valid_types = ["user", "stare", "vad", "wind_profile", "rhi"]
    if not any(valid_type in scan_type for valid_type in valid_types):
        raise NameError(f"Scan type '{scan_type}' not supported.")

    dataset.attrs["scan_type"] = scan_type
    dataset.attrs["z_id"] = z_id

    return dataset

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
                'file_type': 'hpl',
                'ext1':file_types, 
            }
            
            #check if MFA was already in place, otherwise ask to authenticate
            a2e.download_with_order(_filter, path=os.path.join(cd,'data',channel), replace=False)

        #list files to process
        files=glob.glob(os.path.join(cd,'data',channel,'*hpl'))
    
        for f in files:
            try:
                dataset=read(f)
                save_path=f.replace('.00.','.a0.').replace('.00\\','.a0\\').replace('.00/','.a0/').replace('hpl','nc')
                os.makedirs(os.path.dirname(save_path),exist_ok=True)
                dataset.to_netcdf(save_path)
            except:
                print(f+' failed')
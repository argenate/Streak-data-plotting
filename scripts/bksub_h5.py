#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from pathlib import Path
from scipy import ndimage
from scipy.constants import c, h, e, nano
from skimage.measure import block_reduce
import h5py
import json


# Import json with experiment parameters. You have to create the experiment json with experiment.py before running this.
with open('experiment.json', 'r', encoding='utf-8') as f:
    data = json.load(f)
    search = data['sample']
    names = data['names']
    integration = data['integration']
print(names)
print(search)

# Find raw and background data files. They should be stored in the "raw" folder
p = Path('.')
file = list(p.glob(f"raw/*{search}*.txt"))
file = sorted(file)

lengthf = len(file)
files = []
for i in range(lengthf):
    print(f'{i:02d}' + ': ' + str(file[i]))  # Print the file names.
    files.append(file[i])

files2 = files[:]
files3 = files[:]
bks = files2[1::2]
datas = files2[0::2]
print(len(bks))
print(len(datas))

"""
The following code is used to do various corrections
to the data and substract the background.

Many of the commented out lines are used for finding bugs. You may need these.
"""

f = h5py.File(f"corr/{search}_cwl500_700ps_cor_spectra.hdf5", "w")
t_l, t_u = input('Input the time indexes to look for the IRF (e.g 50 150):').split()

for x in range(len(bks)):  # Can use just 1 for testing.
    data = np.loadtxt(str(datas[x]), skiprows=16)
    np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
    # print(data.shape)
    wavelength = data[1:, 0]
    intensity = data[1:, 1:]
    # Photoswitch wavelength calibration. Double check this.
    wavelength = 675 - 0.152 * wavelength - 1.525 * 10**-5 * wavelength**2
    times = data[0, 1:]
    # Photoswitch time calibration
    pix_per_ps = 1.8177 + 0.002 * times - 8 * 10**-7 * times**2
    times = times / pix_per_ps

    bk = np.loadtxt(str(bks[x]), skiprows=16)
    wavelength_bk = bk[1:, 0]
    times_bk = bk[0, 1:]
    intensity_bk = bk[1:, 1:]

    intensity = intensity - intensity_bk  # Subtract background.
    wavelength = wavelength[500:1700]
    times = times[0:1450]
    intensity = intensity[500:1700, 0:1450]
    # print(wavelength.shape, times.shape, intensity.shape)
    # print(np.max(intensity))

    intensity -= np.mean(intensity[10:110, 10:110])
    average1 = np.mean(intensity)
    energy = 1239.88286 / wavelength  # nm2ev

    # shear correction
    shear_m = np.array([[1, 0],
                        [-1E-2, 1]])
    c_out = 0.5 * np.array(intensity.shape)
    c_int = 0.5 * np.array(intensity.shape)
    offset = c_int - c_out.dot(shear_m)
    intensity = ndimage.affine_transform(intensity, shear_m.T, offset=offset)

    # bin width
    # print(times)
    dt = times[1:] - times[:-1]
    # print(dt)
    # assert all(dt>0)
    new_dt = np.zeros_like(times)
    new_dt[:-1] = dt
    new_dt[-1] = dt[-1]
    dwl = energy[1:] - energy[:-1]
    new_dwl = np.zeros_like(energy)
    new_dwl[:-1] = -dwl
    new_dwl[-1] = -dwl[-1]
    intensity /= new_dt[np.newaxis, :]
    intensity /= new_dwl[:, np.newaxis]

    # find t0
    ev1 = (np.abs(energy - 3.0)).argmin()
    ev2 = (np.abs(energy - 3.2)).argmin()
    print(ev1, ev2)
    print(int(t_l), int(t_u))
    irf = intensity[ev1:ev2, int(t_l):int(t_u)].sum(axis=0) * -1  # No idea why this is negative. Have to times by -1.
    # print(f"irf: {irf}")
    max_pix = np.argmax(irf) + int(t_l)
    print(f"max_pix: {max_pix}")
    t0 = times[max_pix]
    print(f"t0: {t0}")
    times -= t0
    # print(f"new t0: {times[0]}")
    average2 = np.mean(intensity)
    intensity *= average1 / average2

    fig2, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(12, 9))
    ax1.plot(times[int(t_l):int(t_u)], irf)
    ax1.set_xlabel('Time (ps)')
    fig2.savefig(f"raw/irf_{search}_{x}.png", dpi=200, facecolor='white', transparent=False)

    # Find shift that keeps peak value closest to t0
    n = 2
    time_mask = slice(0, 1800)
    energy_mask = slice(0, 2100)
    shifts = []
    for i in range(n):
        irf_s = np.sum(irf[i:(-n + i)].reshape((-1, n)), axis=0)
        loc = np.argmax(irf_s)
        t_s = np.mean(times[i:(-n + i)].reshape((-1, n)), axis=0)
        shifts.append(t_s[loc])
    shift = np.argmin(np.abs(shifts))
    # masks
    time_mask = slice(time_mask.start + shift,
                      time_mask.stop - n + shift,
                      time_mask.step)
    # apply shifts
    times = times[time_mask]
    energy = energy[energy_mask]
    intensity = intensity[energy_mask, time_mask]

    # Find the maximum area and normalize.
    ev1 = (np.abs(energy - 2.30)).argmin()  # You may have to change these. Should surround the peak intensity.
    ev2 = (np.abs(energy - 2.46)).argmin()
    time0 = 0
    time1 = (np.abs(times + 5)).argmin()
    time2 = (np.abs(times - 10)).argmin()
    # print(times)
    intensity -= np.average(intensity[ev1:ev2, 0:70])
    intensity /= integration[i]  # Fix for integration time.
    inten_max = np.max(intensity[ev1:ev2, time1:time2])
    intensity /= inten_max

    # Clip the intensity to a certain range to remove spikes. Doesn't include the 3.1eV IRF. Fix the maximum intensity.
    ev1 = (np.abs(energy - 2.0)).argmin()
    ev2 = (np.abs(energy - 3.0)).argmin()  # Will have to change this value for OPA pumped experiments.
    # print(np.max(intensity), np.min(intensity))
    intensity[ev1:ev2, :] = np.clip(intensity[ev1:ev2, :], -0.5, 2)
    # print(np.max(intensity), np.min(intensity))
    intensity *= inten_max

    # Slice the data to remove unnecessary areas.
    ev1 = (np.abs(energy - 2.0)).argmin()
    ev2 = (np.abs(energy - 3.4)).argmin()
    time1 = (np.abs(times + 50)).argmin()
    time2 = (np.abs(times - 300)).argmin()
    intensity = intensity[ev1:ev2, time1:time2]
    energy = energy[ev1:ev2]
    times = times[time1:time2]

    print(f'{x:02d}' + '\n')

    norm = colors.SymLogNorm(linthresh=inten_max / 10, vmin=0, vmax=inten_max, base=10)
    fig3, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(12, 9))
    im = ax1.pcolormesh(times, energy, intensity, cmap='inferno', shading='auto', norm=norm)
    ax1.set_ylabel('Energy (eV)')
    ax1.set_xlabel('Time (ps)')
    fig3.colorbar(im, ax=ax1, format='%g')
    fig3.tight_layout()
    fig3.savefig(f"raw/corr_trace_{search}_{x}.png", dpi=200, facecolor='white', transparent=False)

    # Repack the time and intensity with the trace for output.
    energy = energy.reshape(-1, 1)
    a = np.nan
    times = np.insert(times, 0, a)
    times = times.reshape(-1, 1)
    times = times.T
    # print(energy.shape, times.shape, intensity.shape)
    data = np.concatenate((energy, intensity), axis=1)
    data = np.concatenate((times, data), axis=0)
    print(x)

    dset = f.create_dataset(f'{x}', data=data)
    plt.close('all')

#!/usr/bin/env python
# coding: utf-8

import numpy as np
import json
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.colors as colors
import matplotlib
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy import ndimage, stats
from scipy import signal
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.interpolation import spline_filter
from scipy import fftpack
from scipy.constants import pi, c, h, e, nano, k, Boltzmann
from skimage.measure import block_reduce
from skimage.restoration import denoise_bilateral
import h5py
from scipy import sparse
from scipy.sparse.linalg import spsolve
import numdifftools as nd
import time as timing

from lmfit.models import ConstantModel
from lmfit.model import save_modelresult
from lmfit.models import GaussianModel
from lmfit.models import SkewedGaussianModel
from lmfit.models import VoigtModel
from lmfit.models import LorentzianModel
from lmfit.models import SkewedVoigtModel
from lmfit.models import ExponentialModel, StepModel


"""
This script plots the traces, spectra and decays for each dataset. It should be the first analysis done.

The energy and time values will likely have to be changed for each experiment.

Also exports txt files for easy plotting later.
"""


matplotlib.rcParams.update({'font.size': 25})
matplotlib.use('Agg')
# %matplotlib qt
# matplotlib.use('Qt5Agg')

time0 = timing.perf_counter()
print(time0)

with open('experiment.json', 'r', encoding='utf-8') as f:
    data = json.load(f)
    search = data['sample']
    names = data['names']
    fluences = data['fluences']
    temperatures = data['temperatures']
    integration = data['integration']
print(names)

f = h5py.File(f'{search}_cwl500_700ps_cor_spectra.hdf5', 'r')
print(f.keys())

traces = {}
energys = {}
times = {}

for i in range(len(names)):
    print(i)
    data = f.get(f'{i}')
    data = np.array(data)
    print(data.shape)
    energy = data[1:, 0]
    time = data[0, 1:]
    intensity = data[1:, 1:]

    traces[i] = intensity[:, :]
    energys[i] = energy[:]
    times[i] = time[:]


def batch_trace_plot(list_index):
    """For batch plotting streak traces, spectra, and decays."""
    intensity = traces[i]
    time = times[i]
    evs = energys[i]
    label = names[i]
    print(label)
    print(np.max(intensity), np.min(intensity))

    ev1 = (np.abs(evs - 2.0)).argmin()
    ev2 = (np.abs(evs - 3.3)).argmin()
    ev3 = (np.abs(evs - 2.8)).argmin()
    time1 = (np.abs(time + 50)).argmin()
    time2 = (np.abs(time - 300)).argmin()
    spectrum = np.sum(intensity[ev1:ev2, time1:time2], axis=1)
    decay = np.sum(intensity[ev1:ev3, time1:time2], axis=0)
    evs = evs[ev1:ev2]
    time = time[time1:time2]
    inten = intensity[ev1:ev2, time1:time2]
    # np.savetxt(f'data/decays/{label}_decay.txt', decay)
    # np.savetxt(f'data/decays/{label}_decayt.txt', time)

    # ev1 = (np.abs(evs - 2.8)).argmin()
    # ev2 = (np.abs(evs - 3.0)).argmin()
    # time1 = (np.abs(time - 50)).argmin()
    # time2 = (np.abs(time - 150)).argmin()
    # inten -= np.average(inten[ev1:ev2, time1:time2])

    # # inten = ndimage.median_filter(inten, 5)
    ev1 = (np.abs(evs - 2.31)).argmin()
    ev2 = (np.abs(evs - 2.46)).argmin()
    time1 = (np.abs(time + 10)).argmin()
    time2 = (np.abs(time - 10)).argmin()
    print(np.max(inten[ev1:ev2, time1:time2]))
    inten /= np.max(inten[ev1:ev2, time1:time2])
    # Median filter and added noise to make the plot look a bit cleaner.
    inten = ndimage.median_filter(inten, 5) + np.random.normal(0, 0.05, inten.shape)

    # Batch plot traces
    norm = colors.SymLogNorm(linthresh=0.01, vmin=0.1, vmax=1.0, base=10)
    fig3, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(12, 9))
    im = ax1.pcolormesh(time, evs, inten, cmap='turbo', shading='auto', norm=norm, rasterized=False)
    ax1.set_ylabel('Energy (eV)')
    ax1.set_xlabel('Time (ps)')
    ax1.set_xlim(-5, 50)
    ax1.set_ylim(2.2, 2.7)
    fig3.colorbar(im, ax=ax1, format='%g', ticks=[1, 0.5, 0.25, 0.1])
    fig3.savefig(f"traces/trace_{search}_{label}_2.png", dpi=200, facecolor='white', transparent=False)

    # Turn energy and time into first column and row with NaN as [0,0] placeholder.
    energy = evs.reshape(-1, 1)
    a = np.nan
    time2 = np.insert(time, 0, a)
    time2 = time2.reshape(-1, 1).T
    print(energy.shape, time2.shape, inten.shape)
    data = np.concatenate((energy, inten), axis=1)
    data = np.concatenate((time2, data), axis=0)
    # np.savetxt(f'data/traces/trace_{search}_{label}.txt', data)

    np.savetxt(f'processed/decays/decay_ps_{search}_{label}.txt', decay)
    np.savetxt(f'processed/decays/decay_ps_{search}_{label}t.txt', time)

    # Batch plot spectra and decays.
    fig4, axs = plt.subplots(nrows=1, ncols=1, figsize=(9, 9))
    axs.plot(evs, spectrum)
    axs.set_xlim(2.2, 2.7)
    # axs.set_ylim(-0.001, 1.1)
    axs.set_ylabel('Intensity')
    axs.set_xlabel('Energy (eV)')
    plt.tight_layout()
    fig4.savefig(f"spectra/{search}_{label}_1.png", dpi=200, facecolor='white', transparent=False)

    fig5, axs = plt.subplots(nrows=1, ncols=1, figsize=(9, 9))
    axs.plot(time, decay / np.max(decay))
    axs.set_yscale('log')
    axs.set_xlim(-10, 250)
    axs.set_ylim(0.001, 1.1)
    axs.set_ylabel('Intensity')
    axs.set_xlabel('Time (ps)')
    plt.tight_layout()
    fig5.savefig(f"decays/{search}_{label}_2.png", dpi=200, facecolor='white', transparent=False)
    plt.close('all')

    if (i == 4 or i == 6):  # These should be the indexes for the high and low fluences for comparison.
        comp = spectrum
    else:
        comp = 0

    return comp, evs


compare = {}

for i in range(len(names)):
    # i=6
    compare[i], evs = batch_trace_plot(i)

# Compare and subtract the high and low fluence.
low_flu = compare[4] / np.max(compare[4]) * 0.8  # Multiply to get the correct intensity ratio.
high_flu = compare[6] / np.max(compare[6])
subtract = high_flu - low_flu
subtract /= np.max(subtract)

fig5, axs = plt.subplots(nrows=1, ncols=1, figsize=(12, 9))
axs.plot(evs, high_flu, label='High fluence')
axs.plot(evs, low_flu, label='Low fluence')
axs.plot(evs, subtract, label='Subtracted')
# axs.set_yscale('log')
axs.set_ylim(-0.001, 1.2)
axs.set_xlim(2.2, 2.7)
axs.legend()
fig5.savefig("spectra/spectra_compare.png", dpi=200, facecolor='white',
             transparent=False)
plt.close('all')

time1 = timing.perf_counter()
print(time1 - time0)

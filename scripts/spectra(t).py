#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.colors as colors
import matplotlib
from scipy import ndimage
from scipy import signal
from scipy.constants import pi, c, h, e, nano, k, Boltzmann
from skimage.measure import block_reduce
import h5py
from scipy import sparse
from scipy.sparse.linalg import spsolve
import numdifftools as nd
import time as timing
import json

from lmfit.models import ConstantModel
from lmfit.model import save_modelresult
from lmfit.models import GaussianModel
from lmfit.models import SkewedGaussianModel
from lmfit.models import VoigtModel
from lmfit.models import LorentzianModel
from lmfit.models import SkewedVoigtModel
from lmfit.models import ExponentialModel, StepModel
import corr_spec

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


"""For correcting and putting streaked files into a dictionary."""
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

# ex = pd.read_csv(f'data/{search}_PL.csv', sep=',', names=('nm', 'int'), usecols=(0, 1), header=None, engine='python',
#                  error_bad_lines=False, warn_bad_lines=False)
# ex = ex.dropna(how='any')
# ex = ex[pd.to_numeric(ex['nm'], errors='coerce').notnull()]
# ex = ex.to_numpy(dtype='float32')
# emission = ex[:, 1]
# evs2 = ex[:, 0]
# E = 1239.84193
# evs2 = E / evs2  #  - 0.03
# emission = [z/c**2 for z, c in zip(map(float, emission), map(float, (evs2)))]  # pl nm2ev
# emission /= np.max(emission)


def batch_spectra_t_plot(i):
    """Plotting slices of streak data along time
    (i.e. spectra in time)."""
    print(i)
    intensity = traces[i]
    time = times[i]
    evs = energys[i]
    label = names[i]

    ev1 = (np.abs(evs - 2.8)).argmin()
    ev2 = (np.abs(evs - 3.0)).argmin()
    time1 = (np.abs(time - 50)).argmin()
    time2 = (np.abs(time - 150)).argmin()
    intensity -= np.average(intensity[ev1:ev2, time1:time2])

    ev1 = (np.abs(evs - 2.1)).argmin()
    ev2 = (np.abs(evs - 2.8)).argmin()
    time1 = (np.abs(time + 10)).argmin()
    time2 = (np.abs(time - 200)).argmin()
    evs = evs[ev1:ev2]
    time = time[time1:time2]
    inten = intensity[ev1:ev2, time1:time2]
    # inten = ndimage.median_filter(inten, size=10)
    # label = names[i]

    t_bars = [-10, 5, 10, 30, 60, 100, 250]
    number = len(t_bars)
    t_barsi = []
    for x in range(number):
        t_barsi.append((np.abs(time - t_bars[x])).argmin())
    col = np.flip(plt.cm.jet(np.linspace(0, 0.9, number)), axis=0)

    fig4, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(9, 9))
    ax1.set_xlim(2.2, 2.7)
    for x in range(number - 1):
        spectrum = np.sum(inten[:, t_barsi[x]:t_barsi[x + 1]], axis=1)
        # spectrum = ndimage.median_filter(spectrum, 3)
        # spectrum = signal.savgol_filter(spectrum, 5, 1)
        ev3 = (np.abs(evs - 2.31)).argmin()
        ev4 = (np.abs(evs - 2.46)).argmin()
        spectrum = corr_spec.corr(spectrum)
        spectrum /= np.max(spectrum[ev3:ev4])  # np.mean(decay[time1:time2])
        ax1.plot(evs, spectrum, label=f'{t_bars[x]} - {t_bars[x + 1]} ps', color=col[x])
    # ax1.plot(evs2, emission, label='static', color='black')
    ax1.set_ylim(0, 1.1)  # bottom=(np.mean(decay[time1:time2]) - 0.5), top=(max2 + 0.5))
    ax1.set_xlabel('Energy (eV)')
    ax1.set_ylabel('Intensity (AU)')
    ax1.legend(fontsize=16)
    plt.tight_layout()
    fig4.savefig(f"spectra(t)\\spectra_{label}_t_1.png", dpi=250, facecolor="w")
    plt.close('all')

    return


for i in range(len(names)):
    batch_spectra_t_plot(i)

time1 = timing.perf_counter()

print(time1 - time0)

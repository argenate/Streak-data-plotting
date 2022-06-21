#!/usr/bin/env python
# coding: utf-8

import numpy as np
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

from lmfit.models import *
from lmfit.model import save_modelresult
import lmfit

matplotlib.rcParams.update({'font.size': 25})
# matplotlib.use('Agg')
matplotlib.use('Qt5Agg')
timing0 = timing.perf_counter()

irf = 2.123
t0 = 0
amps = [0.5, 1]
taus = [80, 1000]
params = amps + taus + [t0, irf]


def convolve(arr, kernel):
    """
    Convolution of array with kernel.
    """
    npts = min(len(arr), len(kernel))
    pad = np.ones(npts)
    tmp = np.concatenate((pad*arr[0], arr, pad*arr[-1]))
    norm = np.sum(kernel)
    out = np.convolve(tmp, kernel, mode='valid')
    noff = int((len(out) - npts)/2)
    return out[noff:noff+npts] / norm


def gauss_kernel(x, t0, irf):
    """
    Gaussian convolution kernel.

    Parameters
    ----------
    x : array-like
        Independant variable
    t0 : array-like
        t0 offset
    irf : array-like
        Irf gaussian width (sigma)
    """
    midp = 0.5 * (np.max(x) + np.min(x))
    gauss = lmfit.lineshapes.gaussian(x, 1, midp + t0, irf)
    return gauss


def step(x):
    """Heaviside step function."""
    step = np.ones_like(x, dtype='float')
    step[x<0] = 0
    step[x==0] = 0.5
    return step


def biexp_decay(t, a0, a1, tau0, tau1):
    biexp = step(t)*a0*np.exp(-t/tau0) + step(t)*a1*np.exp(-t/tau1)
    return biexp


def biexp_conv(t, a0, a1, tau0, tau1, t0, irf):
    conv = convolve(biexp_decay(t, a0, a1, tau0, tau1), gauss_kernel(t, t0, irf))
    return conv


def fit_decays(times, decay, irf):
    """
    Function for fitting biexponential decays with an IRF.

    Returns two lifetimes and two amplitudes.

    IRF must be a gaussian sigma value (i.e. FWHM = 2.355*sigma).

    Three inputs are time and intensity and the IRF sigma value.
    """

    init_guess = [1, 1, 10, 100, 0, irf]  # a0, a1, tau0, tau1, t0, irf
    decay_model = lmfit.Model(biexp_conv)
    guess = decay_model.make_params()
    guess["a0"].set(value=init_guess[0], min=0, max=200000)
    guess["a1"].set(value=init_guess[1], min=0, max=200000)
    guess["tau0"].set(value=17, vary=True, min=0.1, max=50)
    guess["tau1"].set(value=100, vary=True, min=20, max=1000)
    guess["t0"].set(value=init_guess[4], min=-10, max=10)
    guess["irf"].set(value=init_guess[5], vary=False, min=1, max=3)
    result = decay_model.fit(decay, guess, t=times, method='dual_annealing',
                             nan_policy='propagate', weights=1/np.max(decay))

    print(result.fit_report())
    best_fit = result.best_fit
    amp0 = result.params['a0'].value
    amp1 = result.params['a1'].value
    tau0 = result.params['tau0'].value
    tau1 = result.params['tau1'].value

    # plt.plot(times, decay, 'bo')
    # plt.plot(times, result.init_fit, 'k--', label='initial fit')
    # plt.plot(times, result.best_fit, 'r-', label='best fit')
    # plt.legend(loc='best')
    # plt.show()

    return tau0, tau1, amp0, amp1, best_fit

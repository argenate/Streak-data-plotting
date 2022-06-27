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


"""
This script should be imported into others for use to do reconvolution fits of fluorescence decays.

Right now the biexponential fit is the most up to date.
The triexponential "fit_tridecays" function must be updated.

See the functions for info on imported values.
"""


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


def triexp_decay(t, a0, a1, a2, tau0, tau1, tau2):
    triexp = step(t)*a0*np.exp(-t/tau0) + step(t)*a1*np.exp(-t/tau1) + step(t)*a2*np.exp(-t/tau2)
    return triexp


def biexp_conv(t, a0, a1, tau0, tau1, t0, irf):
    conv = convolve(biexp_decay(t, a0, a1, tau0, tau1), gauss_kernel(t, t0, irf))
    return conv


def triexp_conv(t, a0, a1, a2, tau0, tau1, tau2, t0, irf):
    conv = convolve(triexp_decay(t, a0, a1, a2, tau0, tau1, tau2), gauss_kernel(t, t0, irf))
    return conv


def fit_bidecays(times, decay, irf, values):
    """
    Function for fitting biexponential decays with an IRF.

    Returns two lifetimes and two amplitudes.

    IRF must be a gaussian sigma value (i.e. FWHM = 2.355*sigma).

    Three inputs are time and intensity and the IRF sigma value.
    
    Example imput of fit parameters:
        values = [amp1, True, life1, True, amp2, False, life2, False]
    where the boolean values decide if the value before are varied or not (False == fixed).
    This example allows amp1 and life1 to vary but fixes amp2 and life2.
    """

    # values = [amp1, True, life1, True, amp2, False, life2, False]


    init_guess = [values[0], values[4], values[2], values[6], 0, irf]  # a0, a1, tau0, tau1, t0, irf
    decay_model = lmfit.Model(biexp_conv)
    guess = decay_model.make_params()
    guess["a0"].set(value=init_guess[0], vary=values[1], min=0.001, max=4)
    guess["a1"].set(value=init_guess[1], vary=values[5], min=0.001, max=4)
    guess["tau0"].set(value=init_guess[2], vary=values[3], min=0.01, max=50)
    guess["tau1"].set(value=init_guess[3], vary=values[7], min=10, max=1000)
    guess["t0"].set(value=init_guess[4], min=-30, max=30)
    guess["irf"].set(value=init_guess[5], vary=True, min=0.1, max=8)
    result = decay_model.fit(decay, guess, t=times, method='dual_annealing',
                             nan_policy='propagate')

    print(result.fit_report())
    best_fit = result.best_fit
    amp0 = result.params['a0'].value
    amp1 = result.params['a1'].value
    tau0 = result.params['tau0'].value
    tau1 = result.params['tau1'].value

    fig = plt.figure(figsize=(10, 10))
    plt.plot(times, decay, 'bo')
    # plt.plot(times, result.init_fit, 'k--', label='initial fit')
    plt.plot(times, result.best_fit, 'r-', label='best fit')
    plt.legend(loc='best')
    fig.savefig('biexp_fit.png')

    return tau0, tau1, amp0, amp1, best_fit, result


def fit_tridecays(times, decay, irf):
    init_guess = [0, 0, 1, 0.05, 0.2, 337, 0, irf]  # a0, a1, a2, tau0, tau1, tau2, t0, irf
    decay_model = lmfit.Model(triexp_conv)
    guess = decay_model.make_params()
    guess["a0"].set(value=init_guess[0], vary=True, min=0.001, max=4)
    guess["a1"].set(value=init_guess[1], vary=True, min=0.001, max=4)
    guess["a2"].set(value=init_guess[2], min=0.001, max=4)
    guess["tau0"].set(value=init_guess[3], vary=True, min=10, max=50)
    guess["tau1"].set(value=init_guess[4], vary=True, min=75, max=200)
    guess["tau2"].set(value=init_guess[5], vary=True, min=300, max=1500)
    guess["t0"].set(value=init_guess[6], min=-30, max=30)
    guess["irf"].set(value=init_guess[7], vary=True, min=0.1, max=8)
    result = decay_model.fit(decay, guess, t=times, method='dual_annealing',
                             nan_policy='propagate')

    print(result.fit_report())
    best_fit = result.best_fit
    amp0 = result.params['a0'].value
    amp1 = result.params['a1'].value
    amp2 = result.params['a2'].value
    tau0 = result.params['tau0'].value
    tau1 = result.params['tau1'].value
    tau2 = result.params['tau2'].value

    # fig = plt.figure(figsize=(10, 10))
    # plt.plot(times, decay, 'bo')
    # # plt.plot(times, result.init_fit, 'k--', label='initial fit')
    # plt.plot(times, result.best_fit, 'r-', label='best fit')
    # plt.legend(loc='best')
    # fig.savefig('biexp_fit.png')

    return tau0, tau1, tau2, amp0, amp1, amp2, best_fit, result.params

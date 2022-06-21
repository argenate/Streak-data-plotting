#!/usr/bin/env python
# coding: utf-8

import numpy as np


def corr(y):
    num_average = 10
    left = np.average(y[:num_average])
    right = np.average(y[(-1 * num_average):])
    left_idx = num_average / 2
    right_idx = y.shape[0] - left_idx

    slope = (right - left) / (right_idx - left_idx)
    intercept = right - slope * right_idx

    y_corr = y - (slope * np.arange(0, y.shape[0]) + intercept)

    return y_corr

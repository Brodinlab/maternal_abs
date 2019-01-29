#!/usr/bin/env python

import numpy as np
import scipy as sp
import scipy.special
import scipy.optimize

params = ['theta', 'lambda']
params_bounds = [(0, None), (0, 1)]

def pmf(x, theta, lambda_):
    return np.exp(logpmf(x, theta, lambda_))
    # return theta * (theta + x*lambda_)**(x - 1) * np.exp(-theta - x*lambda_) / sp.misc.factorial(x)

def logpmf(x, theta, lambda_):
    return np.log(theta) + (x - 1)*np.log(theta + x * lambda_) - (theta + x*lambda_) - sp.special.gammaln(x+1)

def logsf(x, theta, lambda_):
    for chunksize in range(4, 8):
        xs = np.arange(x + 1, x + 10**chunksize)
        logpmf_chunk = logpmf(xs, theta, lambda_)
        accum = np.logaddexp.accumulate(logpmf_chunk)
        if accum[-1] == accum[-2]:
            return accum[-1]
    return np.nan

def lambda_MLE(counts):
    n = len(counts)
    x = np.arange(max(counts) + 1)
    n_x = np.bincount(counts)
    x_bar = np.mean(counts)
    
    is_unique = sum(n_x[2:] * x[2:] * (x[2:] - 1)) - n*(x_bar**2) > 0
    if not is_unique:
        raise ValueError

    return lambda lambda_: sum(n_x * (x*(x-1)/(x_bar + (x-x_bar)*lambda_))) - n*x_bar

def fit(counts):
    try:
        H = lambda_MLE(counts)
    except ValueError:
        return np.nan, np.nan
    # find root by brentq since we know 0 < lambda < 1
    lt1 = 1. - np.finfo(np.float64).epsneg
    lambda_ = sp.optimize.brentq(H, 0., lt1)
    theta = np.mean(counts) * (1 - lambda_)
    return theta, lambda_

def dlogpmf_dtheta(x, theta, lambda_):
    return 1/theta + (x - 1)/(theta + x*lambda_) - 1

def dlogpmf_dlambda(x, theta, lambda_):
    return (x - 1)/(theta + x*lambda_) * x - x

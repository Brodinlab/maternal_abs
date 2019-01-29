#!/usr/bin/env python

from __future__ import division

import numpy as np
import scipy.optimize

import generalizedpoisson as gp

params = ['theta', 'lambda', 'pi']
params_bounds = [(0, None), (0, 1), (0, 1)]

def pmf(x, theta, lambda_, pi):
    return np.where(x == 0,
                    pi + (1 - pi)*gp.pmf(0, theta, lambda_),
                    (1 - pi)*gp.pmf(x, theta, lambda_))


def logpmf(x, theta, lambda_, pi):
    return np.where(x == 0,
                    np.log(pi + (1 - pi)*gp.pmf(0, theta, lambda_)),
                    np.log(1 - pi) + gp.logpmf(x, theta, lambda_))


def logsf(x, theta, lambda_, pi):
    for chunksize in range(4, 8):
        xs = np.arange(x + 1, x + 10**chunksize)
        logpmf_chunk = logpmf(xs, theta, lambda_, pi)
        accum = np.logaddexp.accumulate(logpmf_chunk)
        if accum[-1] == accum[-2]:
            return accum[-1]
    return np.nan


def fit(counts):
    theta0, lambda0 = gp.fit(counts)
    pi0 = np.sum(counts == 0) / len(counts)
    res = scipy.optimize.fmin_l_bfgs_b(loglike, 
                                       x0=[theta0, lambda0, pi0],
                                       fprime=loglike_grad,
                                       # approx_grad=True,
                                       args=(np.bincount(counts),),
                                       bounds=params_bounds,
                                       factr=10.0,
                                       disp=0)
    if res[2]['warnflag'] != 0:
        raise ValueError('Did not converge')
    return res[0]


def loglike(params, n):
    theta, lambda_, pi = params
    m = n.shape[0]
    return -(n[0]*np.log(pi + (1 - pi)*gp.pmf(0, theta, lambda_))
            + (n[1:].sum())*np.log(1 - pi) 
            + np.sum(n[1:] * gp.logpmf(np.arange(1, m), theta, lambda_)))


def loglike_grad(params, n):
    return np.array([dloglike_dtheta(params, n),
                    dloglike_dlambda(params, n),
                    dloglike_dpi(params, n)])


def dloglike_dtheta(params, n):
    theta, lambda_, pi = params
    m = n.shape[0]
    a = pi + (1 - pi)*gp.pmf(0, theta, lambda_)
    b = -(1 - pi)*np.exp(-theta)
    c = np.sum(n[1:] * gp.dlogpmf_dtheta(np.arange(1, m), theta, lambda_))
    return - (n[0]/a*b + c)


def dloglike_dlambda(params, n):
    theta, lambda_, pi = params
    m = n.shape[0]
    return -np.sum(n[1:] * gp.dlogpmf_dlambda(np.arange(1, m), theta, lambda_))


def dloglike_dpi(params, n):
    theta, lambda_, pi = params
    a = pi + (1 - pi)*gp.pmf(0, theta, lambda_)
    b = 1 - np.exp(-theta)
    c = (n[1:].sum()) * -(1/(1 - pi))
    return -(n[0]/a*b + c)

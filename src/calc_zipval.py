#!/usr/bin/env python
from __future__ import division

import sys
import gzip
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import zigp


def fit_models(in_, out, model, **kwargs):
    # only try to fit when there are enough datapoints
    filtered = out.groupby(in_).filter(lambda df: len(df.unique()) > 4)
    # fit parameters for each unique input value
    grouped = filtered.groupby(in_.ix[filtered.index])
    params_fits = grouped.apply(lambda k: pd.Series(fit_model(k, model, **kwargs), index=model.params)).unstack()
    return params_fits


def fit_model(counts, model, remove_outlier=True):
    if remove_outlier:
        try:
            counts = filter_counts(counts, model)
        except ValueError:
            return [np.nan for _ in model.params]

    if len(counts) == 0:
        return [np.nan for _ in model.params]
    
    try:
        params = model.fit(counts)
    except ValueError:
        return [np.nan for _ in model.params]
    
    return params


def filter_counts(counts, model):
    obs = counts.value_counts().sort_index()
    fit = model.pmf(obs.index.values, *model.fit(counts)) * sum(obs)
    
    p = len(model.params)
    threshold = 4 / len(obs)    
    min_outlier = np.inf
    for i in range(1, 4):
        outlier = obs.index.values[-i]
        mask = counts != outlier
        fit1 = model.pmf(obs.index.values, *model.fit(counts[mask])) * sum(mask)
        
        D = sum((fit - fit1)**2) / (p * (sum((obs - fit)**2) / len(obs)))
        if D < threshold:
            break
        min_outlier = outlier
    
    mask = counts < min_outlier
    return counts[mask]


def regress_theta(series):
    mask = np.isfinite(series)
    coeffs = np.polyfit(series.index[mask], series[mask], 1)
    return lambda x: np.polyval(coeffs, x)


def regress_lambda(series):
    mask = np.isfinite(series)
    coeffs = np.polyfit(series.index[mask], series[mask], 1)
    return lambda x: np.polyval(coeffs, x)
    # constant = series.mean()
    # return lambda x: constant * np.ones(*x.shape)


def regress_pi(series):
    def f(x, a, b, c):
        return a * np.exp(-b*x) + c
    mask = np.isfinite(series)
    if mask.sum() < 5:
        return lambda x: np.zeros_like(x)
    x = np.asarray(series.index[mask])
    y = np.asarray(series[mask])
    y0 = [0, 0, 0]
    params, cov = curve_fit(f, x, y, y0, maxfev=int(1e6))
    return lambda x: np.clip(f(x, *params), 0, 1 - np.finfo(float).eps)


param_models = {'theta': regress_theta,
                'lambda': regress_lambda,
                'pi': regress_pi}


def calc_pvals(in_, out, model, param_models):
    if np.any(in_.index != out.index):
        raise ValueError("Input and output indexes are not matched.")    
    params_fits = fit_models(in_, out, model)
    # regress parameters on input count
    param_funcs = {}
    for param in model.params:
        param_funcs[param] = param_models[param](params_fits[param])
    # calculate -log10pval for unique (input count, output count) pairs
    pairs = pd.concat([in_, out], axis=1).astype(int)
    uniq_pairs = pairs.drop_duplicates()
    params_ests = [param_funcs[p](uniq_pairs[in_.name]) for p in model.params]
    uniq_pairs['pval'] = [calc_pval(k, model, *p) for k, p in zip(uniq_pairs[out.name], zip(*params_ests))]
    # merge back to original
    pvals = pd.merge(pairs.reset_index(), uniq_pairs).set_index('id')['pval'].sort_index()
    pvals.name = out.name
    return pvals


def calc_pval(count, model, *params):
    return -1/np.log(10) * model.logsf(count, *params)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_counts")
    parser.add_argument("input_counts")
    parser.add_argument("fig_dir")
    args = parser.parse_args()

    # load data
    out = pd.read_csv(gzip.open(args.sample_counts), index_col=0, squeeze=True)
    in_ = pd.read_csv(gzip.open(args.input_counts), index_col=0, squeeze=True)
    
    # calculate -log10pvals
    pvals = calc_pvals(in_, out, zigp, param_models)
    pvals.sort_index().to_csv(sys.stdout, header=True, float_format='%.2f')

    #
    # diagnostic plots 
    #
    read_count = out.sum()
    # scatter: sample_count vs input_count
    fig, ax = plt.subplots()
    ax.loglog(in_, out, 'o')

    ax.set_xlabel("input count")
    ax.set_ylabel("sample count")
    ax.set_title(out.name)

    fig.tight_layout()
    fig.savefig(args.fig_dir + '/{}.scatter.png'.format(out.name))

    # read count distribution
    fig, ax = plt.subplots()
    (out.value_counts() / len(in_)).sort_index().plot(style='o', ax=ax)

    read_count = out.sum()
    ax.text(0.98, 0.98, '{:0.2e} reads'.format(read_count), transform=ax.transAxes, va='top', ha='right')
    ax.set_title(out.name)
    ax.set_xlim(0, 50)

    fig.tight_layout()
    fig.savefig(args.fig_dir + '/{}.distribution.png'.format(out.name))

    # parameter regression
    model = zigp
    n = len(model.params)
    fig, axes = plt.subplots(1, n, figsize=(4*n,4))

    params_fits = fit_models(in_, out, model)
    param_funcs = {}
    for param in model.params:
        param_funcs[param] = param_models[param](params_fits[param])
    
    for ax, param, bounds in zip(axes, model.params, model.params_bounds):
        params_fits[param].plot(style='o', ax=ax)
        x = np.array(params_fits.index)
        y = np.array(param_funcs[param](params_fits.index))
        ax.plot(x, y, 'r--')
#        ax.set_ylim(min(y.min(), bounds[0]), max(y.max(), bounds[1]))
        ax.set_title(r'$\{}$'.format(param))
    
    fig.tight_layout()
    fig.savefig(args.fig_dir + '/{}.parameter_regression.png'.format(out.name))

    # fits
    min_y, max_y = np.inf, -np.inf
    max_incount = (in_.value_counts().sort_index() < 100).argmax()
    fig, axes = plt.subplots(ncols=2, figsize=(8,4))
    for incount in range(0, max_incount, int(max_incount/2)):
        values = out.ix[in_ == incount]
        obs_pmf = (values.value_counts() / len(values)).sort_index()
        min_y = min(obs_pmf.min(), min_y)
        max_y = max(obs_pmf.max(), max_y)

        obs_pmf.plot(style='o', label='input = {}'.format(incount), ax=axes[0])
        obs_pmf.plot(style='o', label='input = {}'.format(incount), logx=True, logy=True, ax=axes[1])
        color = axes[0].get_lines()[-1].get_color()

        params = [param_funcs[param](np.array([incount])) for param in model.params]
        
        x = np.arange(max(values) + 1)
        fit_pmf = np.exp(model.logpmf(x, *params))
        axes[0].plot(x, fit_pmf, '--', color=color)
        axes[1].loglog(x, fit_pmf, '--', color=color)

    axes[0].set_xlim(0, 25)
    axes[1].set_xlim(0, 100)
    axes[1].set_ylim(min_y, max_y)
    axes[0].set_title(out.name)
    axes[1].set_title(out.name)
    axes[0].legend(loc='best')
    axes[0].legend(loc='best')

    fig.tight_layout()
    fig.savefig(args.fig_dir + '/{}.null_fits.png'.format(out.name))


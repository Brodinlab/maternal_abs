#!/usr/bin/env python
from __future__ import division

import gzip
import argparse

import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("rep1")
parser.add_argument("rep2")
parser.add_argument("fig_dir")
parser.add_argument("threshold", type=float)
args = parser.parse_args()

# load data
rep1 = pd.read_csv(gzip.open(args.rep1), index_col=0, squeeze=True)
rep2 = pd.read_csv(gzip.open(args.rep2), index_col=0, squeeze=True)

# call hits
hits = (rep1 > args.threshold) & (rep2 > args.threshold)
# save
hits.name = rep1.name
hits.to_csv(sys.stdout, header=True)

# diagnostic plots
fig, ax = plt.subplots(figsize=(5,5))
ax.loglog(rep1, rep2, 'o')
ax.loglog(rep1[hits], rep2[hits], 'o', label='hits')

ax.legend(loc='best')
plt.tight_layout()

fig.savefig(args.fig_dir + '/{}.hit.png'.format(rep1.name))

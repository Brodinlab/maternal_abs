#!/usr/bin/env python

from __future__ import print_function

import sys
import gzip
import argparse

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('filenames', nargs='+')
args = parser.parse_args()

tables = []

for filename in args.filenames:
    # print("Loading {}".format(filename), file=sys.stderr)
    # sys.stdout.flush()

    if filename.endswith('.gz'):
        f = gzip.open(filename)
        filename = filename[:-3]
    else:
        f = open(filename)
    
    if filename.endswith('.csv'):
        read_fn = pd.read_csv
    elif filename.endswith('.tsv'):
        read_fn = pd.read_table
    else:
        raise IOError("Unrecognized file extension: {}".format(filename))

    table = read_fn(f, index_col=0, squeeze=True)
    tables.append(table)

concatenated = pd.concat(tables, axis=1)
concatenated.to_csv(sys.stdout, header=True)

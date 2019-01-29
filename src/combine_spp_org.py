#!/usr/bin/env python

from __future__ import print_function

import sys
import gzip
import argparse

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('virscores_spp')
parser.add_argument('virscores_org')
parser.add_argument('metadata')
args = parser.parse_args()

virscores_spp = pd.read_csv(gzip.open(args.virscores_spp), index_col=0)
virscores_org = pd.read_csv(gzip.open(args.virscores_org), index_col=0)
metadata = pd.read_csv(gzip.open(args.metadata))

full_length_viruses = metadata.ix[metadata['Protein names'].str.startswith('Genome polyprotein'), 'Species'].unique()
for species in full_length_viruses:
    organisms = metadata.ix[metadata['Species'] == species, 'Organism'].drop_duplicates()
    virscores_spp.ix[species] = virscores_org.ix[organisms].max()

virscores_spp.to_csv(sys.stdout, header=True)

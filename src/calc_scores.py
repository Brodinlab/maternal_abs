#!/usr/bin/env python
from __future__ import division

import sys
import gzip
import difflib
import argparse

# import regex
import pandas as pd

def is_novel_peptide(peptide, assigned_peptides, epitope_len):
    for assigned_peptide in assigned_peptides:        
        matcher = difflib.SequenceMatcher(None, peptide, assigned_peptide)
        match = matcher.find_longest_match(0, len(peptide), 0, len(assigned_peptide))
        if match.size >= epitope_len:
            return False
    return True

# def is_novel_peptide(peptide, assigned_peptides, epitope_len, max_mm):
#     for i in range(len(peptide) - epitope_len):
#         epitope = peptide[i:i + epitope_len]
#         fmt = r'({epitope}){{s<={max_mm}}}'.format(epitope=epitope, max_mm=max_mm)
#         r = regex.compile(fmt)
#         for assigned_peptide in assigned_peptides:
#             if r.search(assigned_peptide) is not None:
#                 return False
#     return True


def calc_virus_scores(series, level, epitope_len):
    # order viruses by decreasing number of hits
    nhits_per_virus = series.groupby(level=level).sum()
    # initialize all scores to 0
    virus_scores = pd.Series(index=nhits_per_virus.index, name=series.name).fillna(0).astype(int)

    assigned_peptides = set()
    print('hhh')
    print(nhits_per_virus[nhits_per_virus > 0].sort_values(ascending=False).index)
    grouped = series.groupby(level=level)
#    for virus in nhits_per_virus[nhits_per_virus > 0].sort(ascending=False).index:
    for virus in nhits_per_virus[nhits_per_virus > 0].sort_values(ascending=False).index:
        virus_hits = grouped.get_group(virus)

        score = 0
        peptides = virus_hits.index.get_level_values('peptide')
        for peptide in peptides:
            if is_novel_peptide(peptide, assigned_peptides, epitope_len):
            # if is_novel_peptide(peptide, assigned_peptides, epitope_len, max_mm):
                score += 1
                assigned_peptides.add(peptide)
        virus_scores[virus] = score
#        print(score)
    return virus_scores


if __name__ == '__main__':
    print('mmmmmmmm')
    parser = argparse.ArgumentParser()
    parser.add_argument("hits")
    parser.add_argument("oligo_metadata")
    parser.add_argument("beads_nhits")
    parser.add_argument("samps_nhits")
    parser.add_argument("level")
    parser.add_argument("epitope_len", type=int)
    # parser.add_argument("max_mismatch", type=int)
    args = parser.parse_args()

    hits = pd.read_csv(gzip.open(args.hits), index_col=0, squeeze=True)
    lib = pd.read_csv(gzip.open(args.oligo_metadata))
    beads_nhits = pd.read_csv(gzip.open(args.beads_nhits), index_col=0, squeeze=True)
    samps_nhits = pd.read_csv(gzip.open(args.samps_nhits), index_col=0, squeeze=True)
    hits[(beads_nhits > 2) | (samps_nhits < 2)] = 0

    columns = ['id', 'Species', 'Organism', 'Entry', 'peptide']
    hits2 = pd.merge(lib[columns], hits.reset_index(), on='id').set_index(columns)

    # virus_scores = calc_virus_scores(hits2[hits.name], args.level, args.epitope_len, args.max_mismatch)
    print(hits.name)
    virus_scores = calc_virus_scores(hits2[hits.name], args.level, args.epitope_len)
    virus_scores.to_csv(sys.stdout, header=True)

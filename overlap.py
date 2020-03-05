#!/usr/bin/python

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import csv
parser = argparse.ArgumentParser()
parser.add_argument('--sample_one', type=str)
parser.add_argument('--sample_one_name', type=str)
parser.add_argument('--sample_two', type=str)
parser.add_argument('--sample_two_name', type=str)
parser.add_argument('--output', type=str)
parser.add_argument('--union', type=str)
parser.add_argument('--intersection', type=str)
parser.add_argument('--sample_one_unique', type=str)
parser.add_argument('--sample_two_unique', type=str)

args = parser.parse_args()
one = set()
two = set()
with open(args.sample_one, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    fieldnames = list(reader.fieldnames)
    assert('Peptide' in fieldnames)
    for x in reader:
        one.add(x['Peptide'])
with open(args.sample_two, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    fieldnames = list(reader.fieldnames)
    assert('Peptide' in fieldnames)
    for x in reader:
        two.add(x['Peptide'])

venn2([one, two], [args.sample_one_name, args.sample_two_name])
matplotlib.pyplot.savefig(args.output, format='png')

with open(args.union, 'w') as f:
    for x in list(one.union(two)):
        f.write(x + '\n')
with open(args.intersection, 'w') as f:
    for x in list(one.intersection(two)):
        f.write(x + '\n')
with open(args.sample_one_unique, 'w') as f:
    for x in list(one):
        if x not in two:
            f.write(x + '\n')
with open(args.sample_two_unique, 'w') as f:
    for x in list(two):
        if x not in one:
            f.write(x + '\n')

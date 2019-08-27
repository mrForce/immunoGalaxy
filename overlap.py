#!/usr/bin/python

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
parser = argparse.ArgumentParser()
parser.add_argument('--sample_one', type=str)
parser.add_argument('--sample_two', type=str)
parser.add_argument('--output', type=str)
parser.add_argument('--union', type=str)
parser.add_argument('--intersection', type=str)
parser.add_argument('--union', type=str)
parser.add_argument('--sample_one_unique', type=str)
parser.add_argument('--sample_two_unique', type=str)

args = parser.parse_args()
one = set()
two = set()
with open(args.sample_one, 'r') as f:
    for x in f:
        line = x.strip()
        if len(line) > 0:
            one.add(line)
with open(args.sample_two, 'r') as f:
    for x in f:
        line = x.strip()
        if len(line) > 0:
            two.add(line)

venn2([one, two], ['Sample 1', 'Sample 2'])
matplotlib.pyplot.savefig(args.output)

with open(args.union, 'w') as f:
    for x in list(one + two):
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

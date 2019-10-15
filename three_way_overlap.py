#!/usr/bin/python

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
parser = argparse.ArgumentParser()
parser.add_argument('--sample_one', type=str)
parser.add_argument('--sample_one_name', type=str)
parser.add_argument('--sample_two', type=str)
parser.add_argument('--sample_two_name', type=str)
parser.add_argument('--sample_three', type=str)
parser.add_argument('--sample_three_name', type=str)
parser.add_argument('--output', type=str)
parser.add_argument('--union', type=str)
parser.add_argument('--intersection', type=str)
parser.add_argument('--sample_one_unique', type=str)
parser.add_argument('--sample_two_unique', type=str)
parser.add_argument('--sample_three_unique', type=str)

args = parser.parse_args()
one = set()
two = set()
three = set()
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

with open(args.sample_three, 'r') as f:
    for x in f:
        line = x.strip()
        if len(line) > 0:
            three.add(line)

venn3([one, two, three], [args.sample_one_name, args.sample_two_name, args.sample_three_name])
matplotlib.pyplot.savefig(args.output, format='png')

with open(args.union, 'w') as f:
    for x in list(one.union(two).union(three)):
        f.write(x + '\n')
with open(args.intersection, 'w') as f:
    for x in list(one.intersection(two).intersection(three)):
        f.write(x + '\n')
with open(args.sample_one_unique, 'w') as f:
    for x in list(one):
        if x not in two and x not in three:
            f.write(x + '\n')
with open(args.sample_two_unique, 'w') as f:
    for x in list(two):
        if x not in one and x not in three:
            f.write(x + '\n')

with open(args.sample_three_unique, 'w') as f:
    for x in list(three):
        if x not in one and x not in two:
            f.write(x + '\n')

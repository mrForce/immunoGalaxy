#!/usr/bin/python

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
parser = argparse.ArgumentParser()
parser.add_argument('--gal_input', type=str)
parser.add_argument('--uploaded', type=str)
parser.add_argument('--output', type=str)



args = parser.parse_args()
gal_input = set()
uploaded = set()
with open(args.gal_input, 'r') as f:
    for x in f:
        line = x.strip()
        if len(line) > 0:
            gal_input.add(line)
with open(args.uploaded, 'r') as f:
    for x in f:
        line = x.strip()
        if len(line) > 0:
            uploaded.add(line)

venn2([gal_input, uploaded], ['Galaxy', 'Uploaded'])
matplotlib.pyplot.savefig(args.output)


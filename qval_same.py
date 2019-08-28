#!/usr/bin/python
import sys
import argparse
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument('--run', action='append')
parser.add_argument('--output', type=str)

args = parser.parse_args()
master_dict = {}
count_dict = defaultdict(int)
num_runs = len(args.run)
for x in args.run:
    dictionary = {}
    with open(x, 'r') as f:
        for line in f:
            parts = line.split()
            peptide = parts[0].strip()
            q_val = parts[1].strip()
            
            if peptide in master_dict:
                if q_val != master_dict[peptide]:
                    with open(args.output, 'w') as g:
                        g.write('Not equal\n')
                    sys.exit()
                count_dict[peptide] += 1
for k, v in count_dict.items():
    assert v <= num_runs
    if v < num_runs:
        with open(args.output, 'w') as g:
            g.write('Not equal\n')
        sys.exit()
with open(args.output, 'w') as g:
    g.write('The are equal\n')

#!/usr/bin/python

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--output', type=str)
parser.add_argument('--input', type=str)
parser.add_argument('--threshold', type=float)



args = parser.parse_args()

with open(args.input, 'r') as f:
    with open(args.output, 'w') as g:
        header = f.readline()
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                peptide = parts[0]
                q_value = float(parts[1])
                if q_value <= args.threshold:
                    g.write(peptide + '\n')

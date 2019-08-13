#!/usr/bin/python

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--output', type=str)
parser.add_argument('--input', type=str)
parser.add_argument('--threshold', type=float)



args = parser.parse_args()

with open(args.input, 'r') as f:
    with open(args.output, 'w') as g:
        line = f.readline()
        g.write(line)

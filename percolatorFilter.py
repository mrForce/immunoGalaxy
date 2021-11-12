#!/usr/bin/python

import argparse
import os
import zipfile
import fnmatch
import re
import csv
import subprocess
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=str)
parser.add_argument('--input', type=str)
#peptide level Q-value cutoff
parser.add_argument('--threshold', type=float)

args = parser.parse_args()

with open(args.input, 'r') as f:
    with open(args.output, 'w') as g:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = list(reader.fieldnames)
        assert('percolator q-value' in reader.fieldnames)
        assert('sequence' in reader.fieldnames)
        for row in reader:
            q_value = float(row['percolator q-value'])
            peptide = row['sequence'].strip()
            if q_value <= args.threshold:
                g.write(peptide + '\n')

#!/usr/bin/python
from scipy import stats
import argparse
import matplotlib
import math
from Bio.Alphabet import IUPAC
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import tempfile
import subprocess
import re
import os
def read_scores_as_dict(file_path):
    d = dict()
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if ',' in line and len(line):
                split = line.split(',')
                peptide = split[0]
                score = float(split[1])
                assert(peptide not in d)
                d[peptide] = score
            else:
                break
    return d


ptm_removal_regex = re.compile('\[[^\]]*\]')
def create_qq(scores_dict, peptides, length, output_location):
    background_scores = []
    foreground_scores = []
    for peptide in peptides:
        if len(peptide) == length:
            if peptide in scores_dict:
                foreground_scores.append(math.log(scores_dict[peptide]))
            else:
                print('peptide: %s of length %d not in scores_dict' % (peptide, length))
    for k, v in scores_dict.items():
        if len(k) == length:
            background_scores.append(math.log(v))
    quants = list(np.linspace(0, 1, len(foreground_scores) + 2))[1:-1]
    assert(len(quants) == len(foreground_scores))
    background_scores = stats.mstats.mquantiles(background_scores, quants)
    background_scores.sort()
    foreground_scores.sort()
    plt.plot(foreground_scores, background_scores)
    plt.xlim((foreground_scores[0], foreground_scores[-1]))
    plt.ylim((background_scores[0], background_scores[-1]))
    plt.xlabel('Discovered %d-mer Scores' % length)
    plt.ylabel('Background %d-mer Scores' % length)
    plt.title('n = %d' % len(foreground_scores))
    plt.savefig(output_location, format='png')
    plt.clf()
parser = argparse.ArgumentParser()
parser.add_argument('--allele', type=str)
parser.add_argument('--base_project', type=str)
parser.add_argument('--peptides', type=str)
parser.add_argument('--output_eight_mers', type=str)
parser.add_argument('--output_nine_mers', type=str)
parser.add_argument('--output_ten_mers', type=str)

args = parser.parse_args()

eight = tempfile.NamedTemporaryFile()
eight_name = eight.name
nine = tempfile.NamedTemporaryFile()
nine_name = nine.name
ten = tempfile.NamedTemporaryFile()
ten_name = ten.name
base_project = args.base_project

tools_location = '/galaxy-prod/galaxy/tools-dependencies/bin/MSEpitope/tidePipeline'
p = subprocess.Popen(['python3', 'ExportNetMHCScores.py', base_project, args.allele, '8', eight_name], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)
p = subprocess.Popen(['python3', 'ExportNetMHCScores.py', base_project, args.allele, '9', nine_name], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)
p = subprocess.Popen(['python3', 'ExportNetMHCScores.py', base_project, args.allele, '10', ten_name], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)
peptides = set()
with open(args.peptides, 'r') as f:
    for line in f:
        line = line.strip()
        if len(line) > 0:
            cleaned = ptm_removal_regex.sub('', line)
            if all([x in IUPAC.protein.letters for x in cleaned]):
                peptides.add(cleaned)
            else:
                print('Peptide: %s removed, because it contained non-protein character' % cleaned)
print('Collected: %d peptides' % len(peptides))
eight_dict = read_scores_as_dict(eight_name)
eight.close()
nine_dict = read_scores_as_dict(nine_name)
nine.close()
ten_dict = read_scores_as_dict(ten_name)
ten.close()
create_qq(eight_dict, peptides, 8, args.output_eight_mers)
create_qq(nine_dict, peptides, 9, args.output_nine_mers)
create_qq(ten_dict, peptides, 10, args.output_ten_mers)

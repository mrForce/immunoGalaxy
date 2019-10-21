#!/usr/bin/python
from scipy import stats
import argparse
import matplotlib
import math
from Bio.Alphabet import IUPAC
from Bio import SeqIO
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import tempfile
import subprocess
import re
import os
import random

netchop_regex = re.compile('^(\s+(?P<pos>[0-9]+))(\s+[^\s]+){2}(\s+(?P<score>0\.[0-9]+))(\s+(?P<ident>.+))$')

"""
Pass path to a FASTA file, containing proteome. Also pass a list of peptides.

Search the proteome for the peptides, and run NetChop against the proteins that contain a peptide.

Extract the scores for the peptides.

For each protein, select a random set of positions, and get the scores for them.

Returns a tuple (positive_scores, control_scores)
"""
def run_netchop(proteome, peptides, netchop_location, directory):
    sequences = []
    counter = 1
    for record in SeqIO.parse(proteome, 'fasta'):
        sequence = record.seq
        #consists of the locations to 
        positions_to_score = []
        for x in peptides:
            start = sequence.seq.find(x, 0)
            while start > -1:
                positions_to_score.append(start + len(x))
                start = sequence.seq.find(x, start + len(x))
        if positions_to_score:
            sequence.id = str(counter)
            counter += 1
            positions_to_score = set(positions_to_score)
            negative_positions = set(range(1, len(sequence.seq) + 1)) - positions_to_score
            control_positions = set(random.sample(list(negative_positions), len(positions_to_score)))
            sequences[sequence.id] = (sequence, positions_to_score, control_positions)

    output_path = os.path.join(directory, 'netchop_input.fasta')
    output_fp = open(output_path, 'w')
    SeqIO.write([x[1] for x in sequences], output_fp, 'fasta')
    netchop_output_path = os.path.join(directory, 'netchop_output.txt')
    netchop_output_fp = open(netchop_output_path, 'r')
    subprocess.run([netchop_location + ' ' + output_path + ' > ' + netchop_output_path], shell=True)
    positive_scores = []
    control_scores = []
    with open(netchop_output_path, 'r') as f:
        for line in f:
            match = netchop_regex.match(line)
            if match:
                position = int(match.group('pos'))
                score = float(match.group('score'))
                identity = str(match.group('ident'))
                assert(identity in sequences)
                if position in sequences[identity][1]:
                    positive_scores.append(score)
            
    
        
            
    
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

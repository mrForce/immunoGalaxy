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
    sequences = {}
    counter = 1
    for record in SeqIO.parse(proteome, 'fasta'):
        sequence = record.seq
        #consists of the locations to 
        positions_to_score = []
        for x in peptides:
            start = sequence.find(x, 0)
            while start > -1:
                positions_to_score.append(start + len(x))
                start = sequence.find(x, start + len(x))
        if positions_to_score:
            record.id = str(counter)
            record.description = ''
            counter += 1
            positions_to_score = set(positions_to_score)
            negative_positions = set(range(8, len(sequence) + 1)) - positions_to_score
            #control_positions = set(random.sample(list(negative_positions), len(positions_to_score)))
            sequences[record.id] = [record, positions_to_score, negative_positions]
    sequences_values = list(sequences.values())
    min_multiplier = len(sequences_values[0][2])/len(sequences_values[1][1])
    for x in sequences_values:
        ratio = int(len(x[2])/len(x[1]))
        if min_multiplier > ratio:
            min_multiplier = ratio

    print('min multiplier: %d' % min_multiplier)
    assert(min_multiplier > 0)
    for k in sequences.keys():
        random_positions = set(random.sample(list(sequences[k][2]), min_multiplier*len(sequences[k][1])))
        sequences[k][2] = random_positions
    print('sequences')
    print(sequences)
    assert(sequences)
    output_path = os.path.join(directory, 'netchop_input.fasta')
    output_fp = open(output_path, 'w')
    SeqIO.write([x[0] for x in sequences.values()], output_fp, 'fasta')
    output_fp.close()
    netchop_output_path = os.path.join(directory, 'netchop_output.txt')
    subprocess.run([netchop_location + ' ' + os.path.abspath(output_path) + ' > ' + os.path.abspath(netchop_output_path)], shell=True)
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
                if position in sequences[identity][2]:
                    control_scores.append(score)
    
    return (positive_scores, control_scores)
            


ptm_removal_regex = re.compile('\[[^\]]*\]')
def create_qq(positive_scores, control_scores, output_location, num_peptides):
    assert(all([x <= 1.0 and x >= 0 for x in positive_scores]))
    assert(all([x <= 1.0 and x >= 0 for x in control_scores]))
    if len(positive_scores) < len(control_scores):
        quants = list(np.linspace(0, 1, len(positive_scores) + 2))[1:-1]
        assert(len(quants) == len(positive_scores))
        control_scores = stats.mstats.mquantiles(control_scores, quants)
    control_scores.sort()
    positive_scores.sort()
    plt.plot(positive_scores, control_scores)
    plt.xlim((0, 1))
    plt.ylim((0, 1))
    plt.xlabel('Positive Scores')
    plt.ylabel('Control Scores')
    plt.title('n = %d' % num_peptides)
    plt.savefig(output_location, format='png')
    plt.clf()
    
parser = argparse.ArgumentParser()
parser.add_argument('--proteome', type=str)
parser.add_argument('--peptides', type=str)
parser.add_argument('--qq_output', type=str)
parser.add_argument('--temp_directory', type=str)

args = parser.parse_args()

netchop_location = '/galaxy-prod/galaxy/tools-dependencies/bin/MSEpitope/netchop-3.1/netchop'
#netchop_location = '/home/jforce/netchop-3.1/netchop'
peptides = []
with open(args.peptides, 'r') as f:
    for x in f:
        if len(x.strip()) > 0:
            cleaned = ptm_removal_regex.sub('', x.strip())
            peptides.append(cleaned)
assert(len(peptides) > 0)
temp_directory = os.getcwd()
if args.temp_directory:
    temp_directory = args.temp_directory
positive_scores, control_scores = run_netchop(args.proteome, peptides, netchop_location, temp_directory)
create_qq(positive_scores, control_scores, args.qq_output, len(peptides))

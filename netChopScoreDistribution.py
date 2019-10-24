#!/usr/bin/python
from scipy import stats
import argparse
import matplotlib
import collections
import math
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import tempfile
import subprocess
import re
import os
import random


"""
source_string is a string.
position_to_insert_map is a dictionary that maps a position to the string to insert at that position.

Returns a tuple of the form (interleaved_string, positions)

interleaved_string is the new string with the positions inserted.

positions is the new locations of the keys of position_to_insert_map
"""
def interleave(source_string, position_to_insert_map):
    components = []
    temp = ''
    for i in range(0, len(source_string)):
        if i in position_to_insert_map:
            components.append(temp)
            components.append(position_to_insert_map[i])
            temp = source_string[i]
        else:
            temp += source_string[i]
    components.append(temp)
    new_string = ''.join(components)
    new_positions = []
    offset = 0
    for k, v in position_to_insert_map.items():
        offset += len(v)
        new_positions.append(offset + k)                             
        
    return (new_string, new_positions)




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
                positions_to_score.append((start + len(x), x))
                start = sequence.find(x, start + len(x))
        if positions_to_score:
            record.id = str(counter)
            record.description = ''
            counter += 1
            negative_positions = set(range(8, len(sequence) + 1)) - set([x[0] for x in positions_to_score])
            sequences[record.id] = [record, positions_to_score, set([x[0] for x in positions_to_score]), negative_positions]
    sequences_values = list(sequences.values())
    min_multiplier = int(len(sequences_values[0][3])/len(sequences_values[0][2]))
    for x in sequences_values:
        ratio = int(len(x[3])/len(x[2]))
        if min_multiplier > ratio:
            min_multiplier = ratio

    print('min multiplier: %d' % min_multiplier)
    assert(min_multiplier > 0)
    control_paste_sequences = {}
    for k in sequences.keys():
        random_positions = set(random.sample(list(sequences[k][3]), min_multiplier*len(sequences[k][2])))
        paste_random_positions = random.sample(list(sequences[k][3]), len(sequences[k][1]))
        selected_paste_positions = {}
        for i in range(0, len(sequences[k][1])):
            selected_paste_positions[paste_random_positions[i]] = sequences[k][1][i][1]
        new_string, new_positions = interleave(sequences[k][0].seq, selected_paste_positions)
        assert(len(new_positions) == len(list(selected_paste_positions.keys())))
        control_paste_sequences[k] = [SeqRecord(Seq(new_string), id=k, description=''), new_positions]
        sequences[k][1] = random_positions
    assert(sequences)
    output_path = os.path.join(directory, 'netchop_input.fasta')
    output_fp = open(output_path, 'w')
    SeqIO.write([x[0] for x in sequences.values()], output_fp, 'fasta')
    output_fp.close()
    netchop_output_path = os.path.join(directory, 'netchop_output.txt')
    subprocess.run([netchop_location + ' ' + os.path.abspath(output_path) + ' > ' + os.path.abspath(netchop_output_path)], shell=True)

    output_paste_path = os.path.join(directory, 'netchop_paste_input.fasta')
    output_paste_fp = open(output_paste_path, 'w')
    SeqIO.write([x[0] for x in control_paste_sequences.values()], output_paste_fp, 'fasta')
    output_paste_fp.close()
    netchop_paste_output_path = os.path.join(directory, 'netchop_paste_output.txt')
    subprocess.run([netchop_location + ' ' + os.path.abspath(output_paste_path) + ' > ' + os.path.abspath(netchop_paste_output_path)], shell=True)
    positive_scores = []
    control_scores = []
    control_paste_peptide_scores = []
    with open(netchop_output_path, 'r') as f:
        for line in f:
            match = netchop_regex.match(line)
            if match:
                position = int(match.group('pos'))
                score = float(match.group('score'))
                identity = str(match.group('ident'))
                assert(identity in sequences)
                if position in sequences[identity][2]:
                    positive_scores.append(score)
                if position in sequences[identity][1]:
                    control_scores.append(score)
    with open(netchop_paste_output_path, 'r') as f:
        for line in f:
            match = netchop_regex.match(line)
            if match:
                position = int(match.group('pos'))
                score = float(match.group('score'))
                identity = str(match.group('ident'))
                assert(identity in control_paste_sequences)
                if position in control_paste_sequences[identity][1]:
                    control_paste_peptide_scores.append(score)
    assert(sum([len(x[1]) for x in control_paste_sequences.values()]) == len(positive_scores))
    assert(sum([len(x[1]) for x in control_paste_sequences.values()]) == len(control_paste_peptide_scores))
    return (positive_scores, control_scores, control_paste_peptide_scores)
            


ptm_removal_regex = re.compile('\[[^\]]*\]')
def create_qq(positive_scores, control_scores, control_paste_scores, output_location, num_peptides):
    assert(all([x <= 1.0 and x >= 0 for x in positive_scores]))
    assert(all([x <= 1.0 and x >= 0 for x in control_scores]))
    if len(positive_scores) < len(control_scores):
        quants = list(np.linspace(0, 1, len(positive_scores) + 2))[1:-1]
        assert(len(quants) == len(positive_scores))
        control_scores = stats.mstats.mquantiles(control_scores, quants)
    control_scores.sort()
    positive_scores.sort()
    control_paste_scores.sort()
    plt.plot(positive_scores, control_scores, label='Random Position')
    plt.plot(positive_scores, control_paste_scores, label='Random position with peptide pasting')
    plt.xlim((0, 1))
    plt.ylim((0, 1))
    plt.xlabel('Positive Scores')
    plt.ylabel('Control Scores')
    plt.legend()
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
positive_scores, control_scores, control_paste_scores = run_netchop(args.proteome, peptides, netchop_location, temp_directory)
print('control paste scores')
print(control_paste_scores)
print('control scores')
print(control_scores)
print('positive scores')
print(positive_scores)
create_qq(positive_scores, control_scores, control_paste_scores, args.qq_output, len(peptides))

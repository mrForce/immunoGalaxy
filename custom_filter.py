#!/usr/bin/python

import argparse
import re
import itertools
from enum import Enum
import os
import csv
import sys
class CutoffType(Enum):
    Q_VALUE = 1
    FDR = 2

class FileType(Enum):
    COMBINED = 1
    TARGET = 2
    DECOY = 3
def fdr_cutoff(entries, cutoff, score_direction, cutoff_type, peptide_unique = True):
    """
    entries should be a list of dictionaries, each of the form {'score':, 'label':, 'peptide':}

    label is -1 for decoy, 1 for target
    """

    #Actually, first sort so decoys are before targets. That way if they score the same, the decoy will be included in the FDR/Q-value calculation.
    #entries = sorted(entries, key=lambda x: int(x['label']))
    
    #first, sort by score in descending order if score_direction is +, ascending order if score_direction is -
    assert(score_direction in ['+', '-'])
    assert(cutoff_type is CutoffType.Q_VALUE or cutoff_type is CutoffType.FDR)
    for i in range(0, len(entries)):
        entries[i]['index'] = i
    sorted_entries = []
    if score_direction == '+':
        print('reverse order')
        sorted_entries = sorted(entries, key=lambda x: float(x['score']), reverse=True)
    elif score_direction == '-':
        print('forward order')
        sorted_entries = sorted(entries, key=lambda x: float(x['score']))
    assert(sorted_entries)
    
    #if a peptide has multiple entries, take the one with the best score
    peptides = []
    unique_peptide_entries = []
    #assert(len(sorted_entries) > len(set([x['peptide'] for x in sorted_entries])))
    for x in sorted_entries:        
        if (not peptide_unique) or (x['peptide'] not in peptides):
            peptides.append(x['peptide'])
            unique_peptide_entries.append(x)
    if peptide_unique:
        assert(len(unique_peptide_entries) == len(set([x['peptide'] for x in unique_peptide_entries])))
    num_targets = 0
    num_decoys = 1
    indices = []
    for score, group in itertools.groupby(unique_peptide_entries, key=lambda x: float(x['score'])):
        group_list = list(group)
        print('group list')
        print(group_list)
        for entry in group_list:
            if entry['label'] == -1:
                num_decoys += 1
            elif entry['label'] == 1:
                num_targets += 1
            else:
                print('label: %d is an invalid value' % entry['label'])
                assert(False)
        if num_targets == 0:
            fdr = 1.0
        else:
            fdr = 1.0*num_decoys/num_targets
        if cutoff_type is CutoffType.FDR and fdr >= cutoff and num_decoys > 1:
            print('breaking out')
            break
        if fdr < cutoff:
            for entry in group_list:
                if entry['label'] == 1:
                    indices.append(entry['index'])
    return indices
    
def parse_peptide(peptide, peptide_regex, ptm_removal_regex = None):
    match = peptide_regex.match(peptide)
    if match and match.group('peptide'):
        peptide = match.group('peptide')
    if ptm_removal_regex:
        return ptm_removal_regex.sub('', peptide)
    else:
        return peptide
peptide_regex = re.compile('^[A-Z\-]\.(?P<peptide>.*)\.[A-Z\-]$')
ptm_removal_regex = re.compile('\[[^\]]*\]')

parser = argparse.ArgumentParser(description='Given TSV file(s), a peptide column, score column, score direction, target/decoy column, FDR threshold, and output directory, apply the peptide level FDR threshold with slight variations. The second variant is whether to use the first (FDR) or last threshold crossing (Q-value). Output files should be target rows from input TSV file(s) that pass the FDR/Q-value threshold')
parser.add_argument('--peptide_column', help='Which column contains the peptide')
parser.add_argument('--score_column', help='Which column contains the score')
parser.add_argument('--score_direction', help='+ if a higher score is better, - if a lower score is better', choices=['+', '-'])
parser.add_argument('--threshold', help='FDR/Q-value cutoff', type=float)
parser.add_argument('--psm_fdr_output')
parser.add_argument('--psm_q_output')
parser.add_argument('--peptide_fdr_output')
parser.add_argument('--peptide_q_output')
#if combined file input, the following options apply
parser.add_argument('--combined_input_file')
parser.add_argument('--label_column')
parser.add_argument('--target_label')
parser.add_argument('--decoy_label')
#if seperate files input, the following options apply
parser.add_argument('--target_file')
parser.add_argument('--decoy_file')

"""
A wrapper around the dict. This is so we can specify the label without modifying the row contents
"""
class Row:
    def __init__(self, row, label):
        self.row = row
        self.label = label

args = parser.parse_args()
assert(args.peptide_column)
assert(args.score_column)
assert(args.score_direction)
assert(args.threshold)
assert(args.psm_fdr_output)
assert(args.psm_q_output)
assert(args.peptide_fdr_output)
assert(args.peptide_q_output)
if args.combined_input_file:
    assert(args.label_column)
    assert(args.target_label)
    assert(args.decoy_label)
    assert(args.target_file is None)
    assert(args.decoy_file is None)
else:
    assert(args.label_column is None)
    assert(args.target_label is None)
    assert(args.decoy_label is None)
    assert(args.target_file)
    assert(args.decoy_file)

def read_tsv_file(input_path, arguments, file_type, fieldnames = None):
    rows = []
    assert(file_type is FileType.COMBINED or file_type is FileType.TARGET or file_type is FileType.DECOY)
    with open(input_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if fieldnames:
            assert(set(reader.fieldnames).issubset(fieldnames))
            assert(set(reader.fieldnames).issuperset(fieldnames))
        else:
            fieldnames = set(reader.fieldnames)
            assert(fieldnames)
        assert(arguments.peptide_column in fieldnames)
        assert(arguments.score_column in fieldnames)
        if file_type is FileType.COMBINED:
            assert(arguments.label_column in fieldnames)
        for row in reader:
            assert(arguments.peptide_column in row)
            assert(arguments.score_column in row)
            label = None
            if file_type is FileType.COMBINED:
                assert(arguments.label_column in row)
                assert(row[arguments.label_column] == arguments.target_label or row[arguments.label_column] == arguments.decoy_label)
                label = 1 if row[arguments.label_column] == arguments.target_label else -1
            elif file_type is FileType.TARGET:
                label = 1
            elif file_type is FileType.DECOY:
                label = -1
            assert(label == -1 or label == 1)
            rows.append(Row(row, label))
    return (rows, fieldnames)

fieldnames = None
rows = []

if args.combined_input_file:
    rows, fieldnames = read_tsv_file(args.combined_input_file, args, FileType.COMBINED, None)
else:
    target_rows, fieldnames = read_tsv_file(args.target_file, args, FileType.TARGET, None)
    decoy_rows, decoy_fieldnames = read_tsv_file(args.decoy_file, args, FileType.DECOY, fieldnames)
    rows = target_rows + decoy_rows
assert(rows)
assert(fieldnames)
parsed_peptide_rows = []
for row_object in rows:
    row = row_object.row
    label = row_object.label
    parsed_peptide_rows.append({'peptide': parse_peptide(row[args.peptide_column], peptide_regex, ptm_removal_regex), 'label': label, 'score': float(row[args.score_column])})

assert(parsed_peptide_rows)

psm_fdr_indices = fdr_cutoff(parsed_peptide_rows, args.threshold, args.score_direction, CutoffType.FDR, False)
psm_fdr_rows = [rows[i].row for i in psm_fdr_indices]


psm_q_value_indices = fdr_cutoff(parsed_peptide_rows, args.threshold, args.score_direction, CutoffType.Q_VALUE, False)
psm_q_value_rows = [rows[i].row for i in psm_q_value_indices]


peptide_fdr_indices = fdr_cutoff(parsed_peptide_rows, args.threshold, args.score_direction, CutoffType.FDR)
peptide_fdr_rows = [rows[i].row for i in peptide_fdr_indices]

peptide_q_value_indices = fdr_cutoff(parsed_peptide_rows, args.threshold, args.score_direction, CutoffType.Q_VALUE)
peptide_q_value_rows = [rows[i].row for i in peptide_q_value_indices]


def write_rows(rows, fieldnames, output_path):
    with open(output_path, 'w+') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
write_rows(psm_fdr_rows, fieldnames, args.psm_fdr_output)
write_rows(psm_q_value_rows, fieldnames, args.psm_q_output)

write_rows(peptide_fdr_rows, fieldnames, args.peptide_fdr_output)
write_rows(peptide_q_value_rows, fieldnames, args.peptide_q_output)

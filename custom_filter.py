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
                print('label needs to be 1 or -1')
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
"""
    num_targets = 0
    num_decoys = 1
    cutoff_index = -1
    print('length of unique_peptide_entries: %d', len(unique_peptide_entries))
    for i in range(0, len(unique_peptide_entries)):
        entry = unique_peptide_entries[i]
        if entry['label'] == -1:
            num_decoys += 1
        elif entry['label'] == 1:
            num_targets += 1
        if num_targets == 0:
            fdr = 1.0
        else:
            fdr = 1.0*num_decoys/num_targets
        if cutoff_type is CutoffType.FDR and fdr >= cutoff:
            print('breaking out')
            break
        if fdr < cutoff:
            cutoff_index = i
    print('num_targets: %d' % num_targets)
    print('num decoys: %d' % num_decoys)
    print('cutoff index: %d' % cutoff_index)
    if cutoff_index == -1:
        return []
    else:
        return [x['index']  for x in unique_peptide_entries[0:(cutoff_index + 1)] if x['label'] == 1]
"""
    
def parse_peptide(peptide, peptide_regex, ptm_removal_regex = None):
    match = peptide_regex.match(peptide)
    if match and match.group('peptide'):
        peptide = match.group('peptide')
    if ptm_removal_regex:
        return ptm_removal_regex.sub('', peptide)
    else:
        return peptide
print(sys.argv)
sys.exit()
peptide_regex = re.compile('^[A-Z\-]\.(?P<peptide>.*)\.[A-Z\-]$')
ptm_removal_regex = re.compile('\[[^\]]*\]')

parser = argparse.ArgumentParser(description='Given TSV file(s), a peptide column, score column, score direction, target/decoy column, FDR threshold, and output directory, apply the peptide level FDR threshold with slight variations. The first variant is whether or not to remove PTMs before removing peptide duplicates. The second variant is whether to use the first (FDR) or last threshold crossing (Q-value). Output files should be target rows from input TSV file(s) that pass the FDR/Q-value threshold')
parser.add_argument('peptide_column', help='Which column contains the peptide')
parser.add_argument('score_column', help='Which column contains the score')
parser.add_argument('score_direction', help='+ if a higher score is better, - if a lower score is better', choices=['+', '-'])
parser.add_argument('label_column', help='Column that indicates if the row is for a target match (is a 1), or a decoy match (is a -1)')
parser.add_argument('threshold', help='FDR/Q-value cutoff', type=float)
parser.add_argument('output_directory', help='Where to put the output files')
parser.add_argument('input_file', help='Combined targets/decoys file. Or, if seperated, specify a decoys file and treat this as targets')
parser.add_argument('--decoy_input_file')

args = parser.parse_args()
fieldnames = None
rows = []
input_files = [args.input_file]
need_label_column = True
on_decoy_file = False
if args.decoy_input_file:
    input_files.append(args.decoy_input_file)
    need_label_column = False
    args.label_column = 'label'

for input_file in input_files:
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if fieldnames:
            assert(set(reader.fieldnames).issubset(fieldnames))
            assert(set(reader.fieldnames).issuperset(fieldnames))
        else:
            fieldnames = set(reader.fieldnames)
            assert(fieldnames)
        for row in reader:
            assert(args.peptide_column in row)
            assert(args.score_column in row)
            row_copy = dict(row)
            if need_label_column:
                assert(args.label_column in row)
                assert(row[args.label_column] in ['1', '-1'])
            else:
                row_copy['label'] = -1 if on_decoy_file else 1
            rows.append(row_copy)
    if not need_label_column:
        on_decoy_file = True
assert(fieldnames)
print('fieldnames')
print(fieldnames)

if not os.path.isdir(args.output_directory):
    #make sure output_directory isn't a file
    assert(not os.path.exists(args.output_directory))
    os.mkdir(args.output_directory)
    assert(os.path.isdir(args.output_directory))



parsed_peptide_rows = []
rows_no_parsed_peptide = []
for row in rows:
    parsed_peptide_rows.append({'peptide': parse_peptide(row[args.peptide_column], peptide_regex, ptm_removal_regex), 'label': int(row[args.label_column]), 'score': float(row[args.score_column])})
    rows_no_parsed_peptide.append({'peptide': parse_peptide(row[args.peptide_column], peptide_regex), 'label': int(row[args.label_column]), 'score': float(row[args.score_column])})


#PSM level FDR
print('psm fdr')
psm_fdr_indices = fdr_cutoff(rows_no_parsed_peptide, args.threshold, args.score_direction, CutoffType.FDR, False)
psm_fdr_rows = [rows[i] for i in psm_fdr_indices]
print('psm q value')
psm_q_value_indices = fdr_cutoff(rows_no_parsed_peptide, args.threshold, args.score_direction, CutoffType.Q_VALUE, False)
psm_q_value_rows = [rows[i] for i in psm_q_value_indices]


#FDR with parsing
fdr_with_parsing_indices = fdr_cutoff(parsed_peptide_rows, args.threshold, args.score_direction, CutoffType.FDR)
fdr_with_parsing_rows = [rows[i] for i in fdr_with_parsing_indices]
#Q-value with parsing
q_value_with_parsing_indices = fdr_cutoff(parsed_peptide_rows, args.threshold, args.score_direction, CutoffType.Q_VALUE)
q_value_with_parsing_rows = [rows[i] for i in q_value_with_parsing_indices]
#FDR without parsing
fdr_no_parsing_indices = fdr_cutoff(rows_no_parsed_peptide, args.threshold, args.score_direction, CutoffType.FDR)
fdr_no_parsing_rows = [rows[i] for i in fdr_no_parsing_indices]
#Q-value without parsing
q_value_no_parsing_indices = fdr_cutoff(rows_no_parsed_peptide, args.threshold, args.score_direction, CutoffType.Q_VALUE)
q_value_no_parsing_rows = [rows[i] for i in q_value_no_parsing_indices]

def write_rows(rows, fieldnames, output_path):
    with open(output_path, 'w+') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
write_rows(fdr_with_parsing_rows, fieldnames, os.path.join(args.output_directory, 'fdr_with_parsing.txt'))
write_rows(psm_fdr_rows, fieldnames, os.path.join(args.output_directory, 'psm_fdr.txt'))
write_rows(psm_q_value_rows, fieldnames, os.path.join(args.output_directory, 'psm_q_value.txt'))
write_rows(q_value_with_parsing_rows, fieldnames, os.path.join(args.output_directory, 'q_value_with_parsing.txt'))
write_rows(fdr_no_parsing_rows, fieldnames, os.path.join(args.output_directory, 'fdr_no_parsing.txt'))
write_rows(q_value_no_parsing_rows, fieldnames, os.path.join(args.output_directory, 'q_value_no_parsing.txt'))

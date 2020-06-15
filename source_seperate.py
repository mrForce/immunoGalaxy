#!/usr/bin/python

import argparse
import csv
import re
parser = argparse.ArgumentParser()
parser.add_argument('input_file', type=str)
parser.add_argument('base_rows', type=str)
parser.add_argument('added_fasta_rows', type=str)
parser.add_argument('protein_column', type=str)

args = parser.parse_args()

with open(args.input_file, 'r') as f:
    base_output = open(args.base_rows, 'w')
    base_rows = []
    additional_fasta_output = open(args.added_fasta_rows, 'w')
    additional_fasta_rows = []
    reader = csv.DictReader(f, delimiter='\t')
    source_regex = re.compile('.*_source=(?P<source>0|1)')
    assert(args.protein_column in reader.fieldnames)
    for row in reader:
        protein = row[args.protein_column]
        match = source_regex.match(protein)
        if match:
            group = match.group('source')
            assert(group in ['0', '1'])
            if group == '0':
                base_rows.append(row)
            elif group == '1':
                additional_fasta_rows.append(row)
        else:
            print('No source')
            assert(match)


    base_writer = csv.DictWriter(base_output, reader.fieldnames, delimiter='\t')
    base_writer.writeheader()
    for x in base_rows:
        base_writer.writerow(x)
    additional_writer = csv.DictWriter(additional_fasta_output, reader.fieldnames, delimiter='\t')
    additional_fasta_output.writeheader()
    for x in additional_fasta_rows:
        additional_writer.writerow(x)
    

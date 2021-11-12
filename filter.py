#!/usr/bin/python

import argparse
import os
import zipfile
import fnmatch
import re
import csv
import subprocess
import sys

print('cwd: ' + os.getcwd())
command = ['python3', '/galaxy-prod/galaxy/tools/MSEpitope/custom_filter.py']
p = subprocess.Popen(command, stderr=sys.stdout.fileno())
assert(p.wait() == 0)
print('stuff')
parser = argparse.ArgumentParser()
#if percolator, then use the output of the search, which is from percolator. If 'raw', then unpack ZIP file, and compute (peptide level) Q-values ourself.
parser.add_argument('--type', choices=['percolator', 'raw'])
parser.add_argument('--qval_compute', choices=['msgf', 'custom'])
parser.add_argument('--output', type=str)
parser.add_argument('--input', type=str)
#peptide level Q-value cutoff
parser.add_argument('--threshold', type=float)

parser.add_argument('--raw_input_format', choices=['zip', 'pin'])

def fdr_cutoff(entries, cutoff):
    """
    entries should be a list of dictionaries, each of the form {'score':, 'label':, 'peptide':}

    label is -1 for decoy, 1 for target
    """
    #first, sort by score in descending order.
    sorted_entries = sorted(entries, key=lambda x: float(x['score']), reverse=True)
    #if a peptide has multiple entries, take the one with the best score
    peptides = []
    unique_peptide_entries = []
    for x in sorted_entries:
        if x['peptide'] not in peptides:
            peptides.append(x['peptide'])
            unique_peptide_entries.append(x)
    num_targets = 0
    num_decoys = 0
    cutoff_index = -1
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
        if fdr <= cutoff:
            cutoff_index = i
    if cutoff_index == -1:
        return []
    else:
        return list(filter(lambda x: x['label'] == 1, unique_peptide_entries[0:(cutoff_index + 1)]))
    
def parse_peptide(peptide, peptide_regex, ptm_removal_regex):
    match = peptide_regex.match(peptide)
    if match and match.group('peptide'):
        matched_peptide = match.group('peptide')
        return ptm_removal_regex.sub('', matched_peptide)
    else:
        return None
args = parser.parse_args()

if args.type == 'percolator':
    with open(args.input, 'r') as f:
        with open(args.output, 'w') as g:
            reader = csv.DictReader(f, delimiter='\t')
            fieldnames = list(reader.fieldnames)
            assert('percolator q-value' in reader.fieldnames)
            writer = csv.DictWriter(g, fieldnames, delimiter='\t')
            writer.writeheader()
            for row in reader:
                q_value = float(row['percolator q-value'])
                if q_value <= args.threshold:
                    writer.writerow(row)
            """
            header = f.readline()
            for line in f:
                parts = line.split()
                if len(parts) == 2:
                    peptide = parts[0]
                    q_value = float(parts[1])
                    if q_value <= args.threshold:
                        g.write(peptide + '\n')
            """
elif args.type == 'raw':
    pin_paths = []
    if args.raw_input_format == 'zip':
        with zipfile.ZipFile(args.input, 'r') as zip_object:
            zip_object.extractall()
        if args.qval_compute == 'msgf':
            tools_location = '/galaxy-prod/galaxy/tools-dependencies/bin/MSEpitope/tidePipeline'
            #tools_location = '/home/jforce/git/tidePipeline'
            project_directory = ''
            with zipfile.ZipFile(args.input, 'r') as zip_object:
                zip_object.extractall()
            for path, dirs, files in os.walk('.', topdown=True):            
                filtered = fnmatch.filter(files, 'database.db')
                if filtered:
                    project_directory = os.path.abspath(path)
                    break
            print('project directory')
            print(project_directory)
            
            p = subprocess.Popen(['python3', 'FilterQValue.py', project_directory, 'msgf', 'search', str(args.threshold), 'filtered_search', '--peptide_q_value'], cwd=tools_location, stderr=sys.stdout.fileno())
            assert(p.wait() == 0)
            print('output')
            print(os.path.abspath(args.output))
            p = subprocess.Popen(['python3', 'ExportPeptides.py', project_directory, 'FilteredSearchResult', 'filtered_search', os.path.abspath(args.output), '--overwrite', 'True'], cwd=tools_location, stderr=sys.stdout.fileno())
            assert(p.wait() == 0)
            sys.exit()
        elif args.qval_compute == 'custom':
            for path, dirs, files in os.walk('.', topdown=True):
                if fnmatch.fnmatch(path, '*msgfplus_search_results*'):
                    filtered = fnmatch.filter(files, 'search.mzid.pin')
                    if filtered:
                        filtered_with_path = [os.path.join(path, x) for x in filtered]
                        print('filtered')
                        print(filtered_with_path)
                        pin_paths.extend(filtered_with_path)
            #should only have 1 PIN
            print('pin paths')
            print(pin_paths)
            assert(len(pin_paths) == 1)
    elif args.raw_input_format == 'pin':
        pin_paths = [args.input]
    peptide_regex = re.compile('^[A-Z\-]\.(?P<peptide>.*)\.[A-Z\-]$')
    ptm_removal_regex = re.compile('\[[^\]]*\]')
    with open(pin_paths[0], 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        #in PIN file, row after the header row contains direction.
        reader.__next__()
        entries = []
        for row in reader:
            peptide = parse_peptide(row['Peptide'], peptide_regex, ptm_removal_regex)
            if peptide:
                entries.append({'score': float(row['lnEValue']), 'label': int(row['Label']), 'peptide': peptide})
            else:
                print('peptide None')
                print(row)
                assert(None)
        peptides = fdr_cutoff(entries, args.threshold)
        with open(args.output, 'w') as g:
            for x in peptides:
                g.write(x['peptide'] + '\n')

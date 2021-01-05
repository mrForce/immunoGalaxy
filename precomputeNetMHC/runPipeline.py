#!/usr/bin/python
import sys
import argparse
import subprocess
import glob
import tempfile
import shutil
import os
import itertools
import uuid

#locations of needed tools
NETMHC = ''
NETMHCPAN = ''
MSGFPLUS = ''
def getMemory():
    p = subprocess.run(['free'], stdout=subprocess.PIPE)
    s = p.stdout.decode()
    lines = s.split('\n')
    for line in lines[1::]:
        split_line = line.split()
        if split_line[0] == 'Mem:':
            return int(split_line[-1])
    return -1
    
parser = argparse.ArgumentParser()

TEST=False

parser.add_argument('baseDirectory')
parser.add_argument('baseFasta')
parser.add_argument('netmhcScoreDir')
parser.add_argument('netmhcPanScoreDir')
parser.add_argument('--allele', action='append')
parser.add_argument('--additional_proteome', action='append')
parser.add_argument('--mod', action='append')
parser.add_argument('--mode', choices=['filtered', 'unfiltered', 'netMHCPercolator'])
parser.add_argument('--pep_len', type=str)
parser.add_argument('--rank_filter', type=str)
parser.add_argument('--frag_method', type=str)
parser.add_argument('--instrument', type=str)
parser.add_argument('--precursor_tolerance', type=str)
parser.add_argument('--mgf', type=str)
parser.add_argument('--archive', type=str)
parser.add_argument('--msgf_unfiltered', type=str)
.parser.add_argument('--percolator_unfiltered', type=str)
parser.add_argument('--num_matches_per_spectrum', type=int)
parser.add_argument('--minLength', type=int, default=0)
parser.add_argument('--maxLength', type=int, default=0)


args = parser.parse_args()
print('args')
print(args)



memory = getMemory()
print('Memory: ' + str(memory))
assert(memory > 0)



allele_list = []
filtered = False
if args.mode == 'filtered':
    filtered = True
peptide_lengths = []

tempDir = tempfile.mkdtemp(prefix='pipelineRun')
outputPath = os.path.join(tempDir, 'search.mzid')
msgfCommand = ['java', '-Xmx' + str(int(memory/2048)), '-jar', MSGFPLUS,
               '-s', args.mgf,
               '-ignoreMetCleavage', '1',
               '-t', args.precursorTolerance,
               '-tda', '1',
               '-addFeatures', '1',
               '-n', str(args.num_matches_per_spectrum),
               '-m', args.frag_method,
               '-inst', args.instrument,
               '-o', outputPath
               ]
if args.mod:
    modPath = os.path.join(tempDir, 'mods.txt')
    with open(modPath, 'w') as f:
        for mod in args.mod:
            f.write(mod + '\n')
    msgfCommand.extend(['-mod', modPath])

if filtered:
    lengths = [int(x) for  x in args.pep_len.split(',')]
    msgfCommand.extend(['-e', '9'])
    msgfCommand.extend(['-minLength', str(min(lengths))])
    msgfCommand.extend(['-maxLength', str(max(lengths))])
    #TODO: Create chain and score table files for additional proteome, then create filtered FASTA
else:
    msgfCommand.extend(['-e', '0'])
    if args.minLength > 0:
        msgfCommand.extend(['-minLength', str(args.minLength)])
    if args.maxLength > 0:
        msgfCommand.extend(['-maxLength', str(args.maxLength)])
    fasta = args.baseFasta
    if args.additional_proteome:
        command = ['awk', '1', args.baseFasta, args.additional_proteome]
        fasta = os.path.join(tempDir, 'combined.fasta')
        f = open(fasta, 'w')
        proc = subprocess.Popen(command, stdout=f)
        outs, errors = proc.communicate()
    msgfCommand.extend(['-d', fasta])
length_to_allele_to_netmhc_map = {}
if filtered:
    for x in peptide_lengths:
        length_to_allele_to_netmhc_map[x] = {}
        for allele in args.allele:
            length_to_allele_to_netmhc_map[x][allele] = ['Proteome' + x + 'Mers_' + allele]
        

proteome = 'proteome'
fasta_link = os.path.join(project_directory, 'proteome.fasta')
if args.additional_proteome:
    i = 0
    additional_fasta_names = []
    for additional_fasta in args.additional_proteome:
        name = '%d_add' % i
        additional_fasta_names.append(name)
        command = ['python3', 'AddFASTA.py', project_directory, additional_fasta, name]
        print('going to call AddFASTA. Command: %s' % ' '.join(command))
        if TEST:
            print('skipping AddFASTA because of TEST')
        else:
            p = subprocess.Popen(command, cwd=tools_location, stderr=sys.stdout.fileno())
            assert(p.wait() == 0)
    if not filtered:
        """
        If we're doing NetMHC filtering, then there's no point in combining FASTA files.
        """
        command = ['python3', 'ConcatFASTA.py', project_directory, '--source', 'cproteome', 'proteome'] + additional_fasta_names
        print('going to call ConcatFASTA. Command: %s' % ' '.join(command))
        if TEST:
            print('skipping ConcatFASTA because of TEST')
        else:
            p = subprocess.Popen(command, cwd=tools_location, stderr=sys.stdout.fileno())
            assert(p.wait() == 0)
        proteome = 'cproteome'
    else:
        for x in peptide_lengths:
            x = x.strip()
            for fasta in additional_fasta_names:
                peptide_list_name = fasta + '_' + x
                command = ['python3', 'KChop.py', project_directory, fasta, x, peptide_list_name]
                print('going to call KChop. Command: %s' % ' '.join(command))
                if TEST:
                    print('skipping KChop because of TEST')
                else:
                    p = subprocess.Popen(command, cwd=tools_location, stderr=sys.stdout.fileno())
                    assert(p.wait() == 0)
                for allele in args.allele:
                    netmhc_command = ['python3', 'RunNetMHC.py', project_directory, peptide_list_name, allele]
                    print('going to call RunNetMHC. Command: %s' % ' '.join(netmhc_command))
                    if TEST:
                        print('skipping RunNetMHC because of TEST')
                    else:
                        p = subprocess.Popen(netmhc_command, cwd=tools_location, stderr= sys.stdout.fileno())
                        assert(p.wait() == 0)
                    netmhc_name = peptide_list_name + '_' + allele
                    length_to_allele_to_netmhc_map[x][allele].append(netmhc_name)


#use half the available memory. Divide by 1024 to get megabytes. Divide by 2 to cut in half. 
command = ['python3', 'RunMSGFPlusSearch.py', project_directory, 'mgf', 'index', 'search', '--modifications_name', 'mod', '--memory', str(int(memory/2048)), '--thread', '4', '--n', str(args.num_matches_per_spectrum), '--t', args.precursor_tolerance]
if args.minLength > 0:
    command.append('--minLength')
    command.append(str(args.minLength))
if args.maxLength > 0:
    command.append('--maxLength')
    command.append(str(args.maxLength))
print('going to call RunMSGFPlusSearch. Command: %s' % ' '.join(command))
if TEST:
    print('skipping RunMSGFPlusSearch because of TEST')
else:
    p = subprocess.Popen(command, cwd=tools_location, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)
print('ran msgfplus search')

command = ['python3', 'RunPercolator.py', project_directory, 'msgfplus', 'search', 'percolator', '--num_matches_per_spectrum', str(args.num_matches_per_spectrum)]
if args.mode == 'netMHCPercolator':
    command.append('--allele')
    for x in args.allele:
        command.append(x)
print('going to call RunPercolator. Command: %s' % ' '.join(command))
if TEST:
    print('skipping RunPercolator because of TEST')
else:
    p = subprocess.Popen(command, cwd=tools_location, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)

print('ran percolator')
#print('going to call ExportPeptidesWithQValues. Command: %s' % ' '.join(['python3', 'ExportPeptidesWithQValues.py', project_directory, 'percolator', os.path.abspath(args.peptide_output)]))
#p = subprocess.Popen(['python3', 'ExportPeptidesWithQValues.py', project_directory, 'percolator', os.path.abspath(args.peptide_output)], cwd=tools_location, stderr=sys.stdout.fileno())
#assert(p.wait() == 0)
#print('ran peptides with q values')

#print('going to call ExportPSMWithQValues. Command: %s' % ' '.join(['python3', 'ExportPSMWithQValues.py', project_directory, 'percolator', os.path.abspath(args.psm_output)]))
#p = subprocess.Popen(['python3', 'ExportPSMWithQValues.py', project_directory, 'percolator', os.path.abspath(args.psm_output)], cwd=tools_location, stderr=sys.stdout.fileno())
#assert(p.wait() == 0)
#print('got psms')

command = ['python3', 'ExportElliePIN.py', project_directory, 'search', '0.01', '/dev/null', args.msgf_unfiltered]
print('going to call ExportElliePIN. Command: %s' % ' '.join(command))
if TEST:
    print('skipping ExportElliePIN because of TEST')
else:
    p = subprocess.Popen(command, cwd=tools_location, stderr=sys.stdout.fileno())
    if p.wait() != 0:
        print('COULD NOT CREATE ELLIE PIN')



percolator_target_file = glob.glob(project_directory + '/**/percolator.target.psms.txt', recursive=True)
assert(len(percolator_target_file) == 1)
percolator_decoy_file = glob.glob(project_directory + '/**/percolator.decoy.psms.txt', recursive=True)
assert(len(percolator_decoy_file) == 1)

with tempfile.NamedTemporaryFile() as f:
    command = ['awk', 'BEGIN {OFS="\t"} NR==1 {print "Label", $0} NR>1 {print 1, $0}', percolator_target_file[0]]
    p = subprocess.Popen(command, stdout=f)
    assert(p.wait() == 0)
    with tempfile.NamedTemporaryFile() as g:
        command = ['awk', 'BEGIN {OFS="\t"} NR>1 {print -1, $0}', percolator_decoy_file[0]]
        p = subprocess.Popen(command, stdout=g)
        assert(p.wait() == 0)
        command = ['cat', f.name, g.name]
        with open(args.percolator_unfiltered, 'w') as h:
            p = subprocess.Popen(command, stdout=h)
            assert(p.wait() == 0)

"""
#use custom script for FDR calculations.
command = ['python3', '/galaxy-prod/galaxy/tools/MSEpitope/custom_filter.py', 'other', args.ellie_unknowns, '--peptide_column', 'Peptide', '--label_column', 'Label', '--score_column', 'lnEValue', '--target_label', '1', '--decoy_label', '-1', '--score_direction', '+', '--threshold', str(args.ellie_fdr), '--psm_q_output', args.ellie_positives, '--psm_fdr_output', '/dev/null', '--peptide_fdr_output', '/dev/null', '--peptide_q_output', '/dev/null']
print('going to call custom_filter. Command: %s' % ' '.join(command))
p = subprocess.Popen(command, stderr=sys.stdout.fileno())
assert(p.wait() == 0)
print('got ellie positive output')
"""

if args.archive:
    print('project directory: %s' % project_directory)
    print('archive: %s' % args.archive)
    subprocess.run(['zip', '-r', args.archive + '.zip', os.path.join(project_directory, 'database.db'), os.path.join(project_directory, 'percolator_results'), os.path.join(project_directory, 'msgfplus_search_results'), os.path.join(project_directory, 'msgfplus_indices'), os.path.join(project_directory, 'MGF'), os.path.join(project_directory, 'Modifications')])
    #subprocess.run(['zip', '-r', args.archive + '.zip', project_directory])
    shutil.move(args.archive + '.zip', args.archive)
    print('Zip file size: %d' % os.path.getsize(args.archive))

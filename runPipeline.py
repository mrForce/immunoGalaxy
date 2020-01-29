#!/usr/bin/python
import sys
import argparse
import subprocess
import shutil
import os
import uuid
parser = argparse.ArgumentParser()


parser.add_argument('--base_project', action='append')
parser.add_argument('--allele', action='append')
parser.add_argument('--additional_proteome', action='append')
parser.add_argument('--pep_len', type=str)
parser.add_argument('--rank_filter', type=str)
parser.add_argument('--frag_method', type=str)
parser.add_argument('--instrument', type=str)
parser.add_argument('--mgf', type=str)
parser.add_argument('--output', type=str)
parser.add_argument('--archive', type=str)
parser.add_argument('--num_matches_per_spectrum', type=int)

args = parser.parse_args()
print('allele')
print(args.allele)
print('additional proteome')
print(args.additional_proteome)
print('pep len')
print(args.pep_len)
print('output')
print(args.output)

#tools_location = '/galaxy-prod/galaxy/tools-dependencies/bin/MSEpitope/tidePipeline'
tools_location = '/home/jordan/github/galaxy/tools-dependencies/bin/MSEpitope/tidePipeline'

allele_list = []
for x in args.allele:
    allele_list.append('--allele')
    allele_list.append(x)
peptide_lengths = args.pep_len.split(',')
filtered = False
if len(arg*s.allele) > 0:
    assert(len(peptide_lengths) > 0)
    filtered = True
project_directory = os.path.join(os.getcwd(), 'project', 'project')
assert(len(set(args.base_project)) == 1)
base_project = args.base_project[0]
print('going to copy project')
p = subprocess.Popen(['python3', 'CopyProject.py', base_project, project_directory] + allele_list, cwd=tools_location, stderr=sys.stdout.fileno())
print('copying project')
assert(p.wait() == 0)

length_to_allele_to_netmhc_map = {}
if filtered:
    for x in peptide_lengths:
        for allele in args.allele:
            length_to_netmhc_map[x][allele] = ['proteome' + x + 'Mers_' + allele]
            

proteome = 'proteome'
fasta_link = os.path.join(project_directory, 'proteome.fasta')
if args.additional_proteome:
    i = 0
    additional_fasta_names = []
    for additional_fasta in args.additional_proteome:
        name = '%d_add' % i
        additional_fasta_names.append(name)
        print('going to call AddFASTA. Command: %s' % ' '.join(command))
        p = subprocess.Popen(command, cwd=tools_location, stderr=sys.stdout.fileno())
        assert(p.wait() == 0)
    if not filtered:
        """
        If we're doing NetMHC filtering, then there's no point in combining FASTA files.
        """
        command = ['python3', 'ConcatFASTA.py', project_directory, 'cproteome', 'proteome'] + additional_fasta_names
        print('going to call ConcatFASTA. Command: %s' % ' '.join(command))
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
                p = subprocess.Popen(command, cwd=tools_location, stderr=sys.stdout.fileno())
                assert(p.wait() == 0)
                for allele in args.allele:
                    netmhc_command = ['python3', 'RunNetMHC.py', project_directory, peptide_list_name, allele]
                    print('going to call RunNetMHC. Command: %s' % ' '.join(netmhc_command))
                    p = subprocess.Popen(netmhc_command, cwd=tools_location, stderr= sys.stdout.fileno())
                    assert(p.wait() == 0)
                    netmhc_name = peptide_list_name + '_' + allele
                    length_to_allele_to_netmhc_map[x][allele].append(netmhc_name)

print('current working: ' + os.getcwd())
print('listdir')
print(os.listdir('.'))
filtered_netmhc_names = []
if filtered:
    for x in peptide_lengths:
        for allele in args.allele:
            filtered_name = 'joined' + x + 'Mers_' + allele + '_filtered'
            #notice how I expect FilterNetMHC.py to work. Change this in tidePipeline
            filter_command = ['python3', 'FilterNetMHC.py', project_directory, args.rank_filter, filtered_name] + length_to_allele_to_netmhc_map[x][allele]
            print('going to call FilterNetMHC. Command: %s' % ' '.join(filter_command))
            p = subprocess.Popen(filter_command, cwd=tools_location, stderr=sys.stdout.fileno())
            assert(p.wait() == 0)
            filtered_netmhc_names.append(filtered_name)
add_tools = lambda x: os.path.join(tools_location, x)


#p = subprocess.Popen(['python3', 'Initialize.py', project_directory, '../pipeline_config.ini', '../unimod.xml'], cwd=tools_location, stderr=sys.stdout.fileno())

#assert(p.wait() == 0)

"""
for x in args.allele:
    p = subprocess.Popen(['python3', 'AddHLA.py', project_directory, x], cwd=tools_location, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)
"""

mgf_link = os.path.join(project_directory, 'thing.mgf')

os.symlink(args.mgf, mgf_link)

print('going to call AddMGF. Command: %s' % ' '.join(['python3', 'AddMGF.py', project_directory, mgf_link, 'mgf', '8', args.frag_method, args.instrument]))
p = subprocess.Popen(['python3', 'AddMGF.py', project_directory, mgf_link, 'mgf', '8', args.frag_method, args.instrument], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)


if filtered_netmhc_names:
    print('going to call CreateTargetSet. Command: %s' % ' '.join(['python3', 'CreateTargetSet.py', project_directory, 'thing'] + filtered_netmhc_names))
    p = subprocess.Popen(['python3', 'CreateTargetSet.py', project_directory, 'thing'] + filtered_netmhc_names, cwd=tools_location, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)
    print('created target set')

    print('going to call CreateMSGFPlusIndex. Command: %s' % ' '.join(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'TargetSet', 'thing', 'index']))
    p = subprocess.Popen(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'TargetSet', 'thing', 'index'], cwd=tools_location, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)
else:
    print('going to call CreateMSGFPlusIndex. Command: %s' % ' '.join(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'FASTA', 'proteome', 'index']))
    p = subprocess.Popen(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'FASTA', 'human', 'index'], cwd=tools_location, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)
    
print('created msgfplus index')
print('going to call RunMSGFPlusSearch. Command: %s' % ' '.join(['python3', 'RunMSGFPlusSearch.py', project_directory, 'mgf', 'index', 'search',  '--memory', '10000', '--thread', '4', '--n', str(args.num_matches_per_spectrum)]))
p = subprocess.Popen(['python3', 'RunMSGFPlusSearch.py', project_directory, 'mgf', 'index', 'search',  '--memory', '10000', '--thread', '4', '--n', str(args.num_matches_per_spectrum)], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)
print('ran msgfplus search')

print('going to call RunPercolator. Command: %s' % ' '.join(['python3', 'RunPercolator.py', project_directory, 'msgfplus', 'search', 'percolator']))
p = subprocess.Popen(['python3', 'RunPercolator.py', project_directory, 'msgfplus', 'search', 'percolator'], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)

print('ran percolator')
print('going to call ExportPeptidesWithQValues. Command: %s' % ' '.join(['python3', 'ExportPeptidesWithQValues.py', project_directory, 'percolator', os.path.abspath(args.output)]))
p = subprocess.Popen(['python3', 'ExportPeptidesWithQValues.py', project_directory, 'percolator', os.path.abspath(args.output)], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)
print('ran peptides with q values')


if args.archive:
    print('project directory: %s' % project_directory)
    print('archive: %s' % args.archive)
    subprocess.run(['zip', '-r', args.archive + '.zip', os.path.join(project_directory, 'percolator_results'), os.path.join(project_directory, 'msgfplus_search_results'), os.path.join(project_directory, 'msgfplus_indices'), os.path.join(project_directory, 'TargetSet'), os.path.join(project_directory, 'FilteredNetMHC')])
    shutil.move(args.archive + '.zip', args.archive)
    print('Zip file size: %d' % os.path.getsize(args.archive))

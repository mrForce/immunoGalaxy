#!/usr/bin/python
import sys
import argparse
import subprocess
import shutil
import os
import uuid
import tempfile
parser = argparse.ArgumentParser()

parser.add_argument('--base_project', action='append')
parser.add_argument('--allele', action='append')
parser.add_argument('--additional_proteome', action='append')
parser.add_argument('--pep_len', type=str)
parser.add_argument('--rank_filter', type=str)
parser.add_argument('--frag_method', type=str)
parser.add_argument('--instrument', type=str)
parser.add_argument('--mgf', type=str)
parser.add_argument('--peptides', type=str)
parser.add_argument('--fdr', type=float)
parser.add_argument('--archive', type=str)


args = parser.parse_args()
print('allele')
print(args.allele)
print('additional proteome')
print(args.additional_proteome)
print('pep len')
print(args.pep_len)
print('output')
print(args.output)

tools_location = '/galaxy-prod/galaxy/tools-dependencies/bin/MSEpitope/tidePipeline'
allele_list = []
for x in args.allele:
    allele_list.append('--allele')
    allele_list.append(x)
project_directory = os.path.join(os.getcwd(), 'project', 'project')
assert(len(set(args.base_project)) == 1)
base_project = args.base_project[0]
print('going to copy project')
p = subprocess.Popen(['python3', 'CopyProject.py', base_project, project_directory] + allele_list, cwd=tools_location, stderr=sys.stdout.fileno())
print('copying project')
fdr = str(args.fdr/100.0)
filtered_netmhc_names = []
assert(p.wait() == 0)
print('going to filter netMHC')
for allele in args.allele:
    for pep_len in args.pep_len.split(','):
        netmhc_name = 'Human' + pep_len + 'Mers_' + allele
        print('going to call FilterNetMHC. Command: %s' % ' '.join(['python3', 'FilterNetMHC.py', project_directory, netmhc_name, args.rank_filter, netmhc_name + '_filtered']))
        p = subprocess.Popen(['python3', 'FilterNetMHC.py', project_directory, netmhc_name, args.rank_filter, netmhc_name + '_filtered'], cwd=tools_location, stderr=sys.stdout.fileno())
        assert(p.wait() == 0)
        print('called FilterNetMHC.py')
        filtered_netmhc_names.append('--FilteredNetMHC')
        filtered_netmhc_names.append(netmhc_name + '_filtered')

print('filtered netMHC')
print('current working: ' + os.getcwd())
print('listdir')
print(os.listdir('.'))

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


fasta_link = os.path.join(project_directory, 'proteome.fasta')
if args.additional_proteome:
    with open(fasta_link, 'w') as f:
        for x in args.additional_proteome:
            with open(x, 'r') as g:
                shutil.copyfileobj(g, f)
    print('going to call AddFASTA. Command: %s' % ' '.join(['python3', 'AddFASTA.py', project_directory, fasta_link, 'proteome']))
    p = subprocess.Popen(['python3', 'AddFASTA.py', project_directory, fasta_link, 'proteome'], cwd=tools_location, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)
    for x in args.pep_len.split(','):
        x = x.strip()
        print('going to call KChop. Command: %s' % ' '.join(['python3', 'KChop.py', project_directory, 'proteome', x, 'proteome' + x]))
        p = subprocess.Popen(['python3', 'KChop.py', project_directory, 'proteome', x, 'proteome' + x], cwd=tools_location, stderr=sys.stdout.fileno())
        
        assert(p.wait() == 0)


    for x in args.allele:
        for y in args.pep_len.split(','):
            y = y.strip()
            netmhc_name = 'proteome' + y + '_' + x
            print('going to call RunNetMHC. Command: %s' % ' '.join(['python3', 'RunNetMHC.py', project_directory, 'proteome' + y, x, args.rank_filter]))
            p = subprocess.Popen(['python3', 'RunNetMHC.py', project_directory, 'proteome' + y, x, args.rank_filter], cwd=tools_location, stderr=sys.stdout.fileno())
            assert(p.wait() == 0)
            print('going to call FilterNetMHC. Command: %s' % ' '.join(['python3', 'FilterNetMHC.py', project_directory, netmhc_name, args.rank_filter, netmhc_name + '_filtered']))
            p = subprocess.Popen(['python3', 'FilterNetMHC.py', project_directory, netmhc_name, args.rank_filter, netmhc_name + '_filtered'], cwd=tools_location, stderr=sys.stdout.fileno())
            assert(p.wait() == 0)
            filtered_netmhc_name = netmhc_name + '_filtered'
            filtered_netmhc_names.append('--FilteredNetMHC')
            filtered_netmhc_names.append(filtered_netmhc_name)

print('filtered netmhc')
print(filtered_netmhc_names)

print('going to call CreateTargetSet. Command: %s' % ' '.join(['python3', 'CreateTargetSet.py', project_directory, 'thing'] + filtered_netmhc_names))
p = subprocess.Popen(['python3', 'CreateTargetSet.py', project_directory, 'thing'] + filtered_netmhc_names, cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)
print('created target set')

print('going to call CreateMSGFPlusIndex. Command: %s' % ' '.join(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'TargetSet', 'thing', 'filtered_index']))
p = subprocess.Popen(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'TargetSet', 'thing', 'filtered_index'], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)


print('going to call CreateMSGFPlusIndex. Command: %s' % ' '.join(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'FASTA', 'proteome', 'index']))
p = subprocess.Popen(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'FASTA', 'human', 'index'], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)

print('created msgfplus indices')
print('Going to add percolator parameter file')
with tempfile.TemporaryFile() as f:
    f.write('top_matches=1\n')
    p = subprocess.Popen(['python3', 'AddParamFile.py', project_directory, f.name, 'percolator'], cwd=tools_location, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)

print('going to call RunIterativeMSGFPlusSearch. Command: %s' % ' '.join(['python3', 'RunIterativeMSGFPlusSearch.py', project_directory, 'mgf', fdr, 'iterative_search', 'index', 'filtered_index', '--use_percolator', 'percolator',  '--memory', '10000']))
p = subprocess.Popen(['python3', 'RunIterativeMSGFPlusSearch.py', project_directory, 'mgf', fdr, 'iterative_search', 'index', 'filtered_index', '--use_percolator', 'percolator',  '--memory', '10000'], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)
print('ran iterative msgfplus search')

print('going to call ExportPeptides. Command: %s' % ' '.join(['python3', 'ExportPeptides.py', project_directory, 'MSGFPlusIterativeSearch', 'iterative_search', args.peptides]))
p = subprocess.Popen(['python3', 'ExportPeptides.py', project_directory, 'MSGFPlusIterativeSearch', 'iterative_search', args.peptides], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)


if args.archive:
    print('project directory: %s' % project_directory)
    print('archive: %s' % args.archive)
    subprocess.run(['zip', '-r', args.archive + '.zip', os.path.join(project_directory, 'percolator_results'), os.path.join(project_directory, 'msgfplus_search_results'), os.path.join(project_directory, 'msgfplus_indices'), os.path.join(project_directory, 'TargetSet'), os.path.join(project_directory, 'FilteredNetMHC'), os.path.join(project_directory, 'FASTA')])
    shutil.move(args.archive + '.zip', args.archive)
    print('Zip file size: %d' % os.path.getsize(args.archive))

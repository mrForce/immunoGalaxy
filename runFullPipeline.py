#!/usr/bin/python
import sys
import argparse
import subprocess
import shutil
import os
import uuid
parser = argparse.ArgumentParser()

parser.add_argument('--base_project', action='append')
parser.add_argument('--frag_method', type=str)
parser.add_argument('--instrument', type=str)
parser.add_argument('--mgf', type=str)
parser.add_argument('--output', type=str)



args = parser.parse_args()
print('output')
print(args.output)

tools_location = '/galaxy-prod/galaxy/tools-dependencies/bin/MSEpitope/tidePipeline'
project_directory = os.path.join(os.getcwd(), 'project', 'project')
assert(len(set(args.base_project)) == 1)
base_project = args.base_project[0]
print('going to copy project')
p = subprocess.Popen(['python3', 'CopyProject.py', base_project, project_directory], cwd=tools_location, stderr=sys.stdout.fileno())
print('copying project')



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



print('going to call CreateMSGFPlusIndex. Command: %s' % ' '.join(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'FASTA', 'proteome', 'index']))
p = subprocess.Popen(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'FASTA', 'proteome', 'index'], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)

print('created msgfplus index')
print('going to call RunMSGFPlusSearch. Command: %s' % ' '.join(['python3', 'RunMSGFPlusSearch.py', project_directory, 'mgf', 'index', 'search',  '--memory', '10000', '--thread', '4']))
p = subprocess.Popen(['python3', 'RunMSGFPlusSearch.py', project_directory, 'mgf', 'index', 'search',  '--memory', '10000', '--thread', '4'], cwd=tools_location, stderr=sys.stdout.fileno())
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


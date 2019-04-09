#!/usr/bin/python
import sys
import argparse
import subprocess
import os
import uuid
parser = argparse.ArgumentParser()

parser.add_argument('--proteome', type=str)
parser.add_argument('--allele', action='append')
parser.add_argument('--pep_len', type=str)
parser.add_argument('--rank_filter', type=str)
parser.add_argument('--frag_method', type=str)
parser.add_argument('--instrument', type=str)
parser.add_argument('--mgf', type=str)
parser.add_argument('--output', type=str)


args = parser.parse_args()
tools_location = '/galaxy-prod/galaxy/tools-dependencies/bin/MSEpitope/tidePipeline'


project_directory = os.path.join(os.getcwd(), 'project')

    
print('current working: ' + os.getcwd())
print('listdir')
print(os.listdir('.'))

add_tools = lambda x: os.path.join(tools_location, x)


p = subprocess.Popen(['python3', 'Initialize.py', project_directory, '../pipeline_config.ini', '../unimod.xml'], cwd=tools_location, stderr=sys.stdout.fileno())

assert(p.wait() == 0)


for x in args.allele:
    p = subprocess.Popen(['python3', 'AddHLA.py', project_directory, x], cwd=tools_location, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)


mgf_link = os.path.join(project_directory, 'thing.mgf')
os.link(args.mgf, mgf_link)


p = subprocess.Popen(['python3', 'AddMGF.py', project_directory, mgf_link, 'mgf', '8', args.frag_method, args.instrument], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)


fasta_link = os.path.join(project_directory, 'proteome.fasta')
os.link(args.proteome, fasta_link)
p = subprocess.Popen(['python3', 'AddFASTA.py', project_directory, fasta_link, 'proteome'], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)


for x in args.pep_len.split(','):
    x = x.strip()
    p = subprocess.Popen(['python3', 'KChop.py', project_directory, 'proteome', x, 'proteome' + x], cwd=tools_location, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)


for x in args.allele:
    for y in args.pep_len.split(','):
        y = y.strip()        
        p = subprocess.Popen(['python3', 'RunNetMHC.py', project_directory, 'proteome' + y, x, args.rank_filter], cwd=tools_location, stderr=sys.stdout.fileno())
        assert(p.wait() == 0)


p = subprocess.Popen(['python3', 'CreateTargetSetFromAllNetMHC.py', project_directory, 'thing'], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)

p = subprocess.Popen(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'TargetSet', 'thing', 'index'], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)

p = subprocess.Popen(['python3', 'RunMSGFPlusSearch.py', project_directory, 'mgf', 'index', 'search',  '--memory', '10000', '--thread', '4'], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)

p = subprocess.Popen(['python3', 'RunPercolator.py', project_directory, 'msgfplus', 'search', 'percolator'], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)

p = subprocess.Popen(['python3', 'ExportPeptidesWithQValues.py', project_directory, 'percolator', args.output], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)


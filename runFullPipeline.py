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
parser.add_argument('--archive', type=str)


args = parser.parse_args()
print('output')
print(args.output)
print('base project')
print(args.base_project)
tools_location = '/galaxy-prod/galaxy/tools-dependencies/bin/MSEpitope/tidePipeline'
#p = subprocess.Popen(['ls', '-R', os.getcwd()], stdout=sys.stdout.fileno())
project_directory = os.path.join(os.getcwd(), 'project', 'project')
assert(len(set(args.base_project)) == 1)
base_project = args.base_project[0]
print('going to copy project')
p = subprocess.Popen(['python3', 'CopyProject.py', base_project, project_directory], cwd=tools_location, stderr=sys.stdout.fileno())
assert(p.wait() == 0)
print('copying project')

assert os.path.exists(project_directory)


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



print('going to call CreateMSGFPlusIndex. Command: %s' % ' '.join(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'FASTA', 'proteome', 'index']))
p = subprocess.Popen(['python3', 'CreateMSGFPlusIndex.py', project_directory, 'FASTA', 'human', 'index'], cwd=tools_location, stderr=sys.stdout.fileno())
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

if args.archive:
    print('project directory: %s' % project_directory)
    print('archive: %s' % args.archive)
    #subprocess.run(['zip', '-r', args.archive + '.zip', os.path.join(project_directory, 'percolator_results'), os.path.join(project_directory, 'msgfplus_search_results'), os.path.join(project_directory, 'msgfplus_indices'), os.path.join(project_directory, 'TargetSet')])
    subprocess.run(['zip', '-r', args.archive + '.zip', project_directory])
    shutil.move(args.archive + '.zip', args.archive)
    print('Zip file size: %d' % os.path.getsize(args.archive))

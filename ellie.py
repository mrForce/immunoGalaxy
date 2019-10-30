#!/usr/bin/python

import argparse
import zipfile
import os
import sys
import fnmatch
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument('--archive', type=str)
parser.add_argument('--positive_fdr', type=float)
parser.add_argument('--positives_output', type=str)
parser.add_argument('--unknown_output', type=str)

args = parser.parse_args()

tools_location = '/galaxy-prod/galaxy/tools-dependencies/bin/MSEpitope/tidePipeline'

with zipfile.ZipFile(args.archive, 'r') as zip_object:
    zip_object.extractall()

project_path = None
for path, dirs, files in os.walk('.', topdown=True):
    if len(fnmatch.filter(files, 'database.db')) > 0:
        project_path = path
        break

    
assert(project_path)

p = subprocess.Popen(['python3', 'ExportElliePIN.py', os.path.abspath(project_path), 'search', str(args.positive_fdr), args.positives_output, args.unknown_output], cwd=tools_location, stderr=sys.stdout.fileno())        
assert(p.wait() == 0)

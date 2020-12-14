#!/usr/bin/python3.9
import shutil
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-a', nargs=1)
parser.add_argument('-p', action='store_true')
parser.add_argument('-xls', action='store_true')
parser.add_argument('-xlsfile', nargs=1)

parser.add_argument('-f', nargs=1)
args = parser.parse_args()
print(args.xlsfile)
shutil.copy('output.xls', args.xlsfile[0])
print('hello')
print('stderr', file=sys.stderr)
sys.exit(1)

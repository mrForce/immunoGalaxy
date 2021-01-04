import io
import argparse
import random
import os
import pickle
import sys
import functools
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable, ScoreCategory
from netMHCCalling import NetMHCScorer, NetMHCRunFailedError
from filterNetMHC import extractPeptidesAndScores

parser = argparse.ArgumentParser(description='List alleles contained in score table')
parser.add_argument('scoreTable')
args=  parser.parse_args()
f = open(args.scoreTable, 'rb')
st = ScoreTable.readExisting(f)
alleles = st.getAlleles()
for x in alleles:
    print(x)
f.close()

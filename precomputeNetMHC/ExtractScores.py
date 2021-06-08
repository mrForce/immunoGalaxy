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

parser = argparse.ArgumentParser(description='Extract a randomly selected sample of peptides and their score')
parser.add_argument('fasta')
parser.add_argument('chains')
parser.add_argument('scoreTable')
parser.add_argument('allele')
parser.add_argument('k', type=int)
args=  parser.parse_args()

assert(os.path.isfile(args.chains))
assert(os.path.isfile(args.fasta))
assert(os.path.isfile(args.scoreTable))
assert(args.k > 0)
chainCollection = None
with open(args.chains, 'rb') as f:
    chainCollection = pickle.load(f)

f = open(args.scoreTable, 'rb')
st = ScoreTable.readExisting(f)
print('alleles: ')
print(st.getAlleles())
assert(args.allele in st.getAlleles())

length = st.peptideLength
n = st.numPeptides
assert(args.k <= n)
print('num peptides: ' + str(n))
pepGen = peptideGenerator(chainCollection, args.fasta, length)
selection = sorted(random.sample(range(0, n), args.k))
scoreIter = st.scoreIter(args.allele)
pepAndScores = extractPeptidesAndScores(scoreIter, pepGen, selection)
for pep, score in pepAndScores:
    print(str(pep) + '\t' + str(score))

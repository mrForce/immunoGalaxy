import io
import argparse
import os
import pickle
import sys
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable
from netMHCCalling import NetMHCScorer, NetMHCRunFailedError

parser = argparse.ArgumentParser(description='Precompute NetMHC scores for peptides')
parser.add_argument('netMHC')
parser.add_argument('fasta')
parser.add_argument('chains')
parser.add_argument('scoreTable')

parser.add_argument('allele')
parser.add_argument('length', type=int)


args = parser.parse_args()
assert(args.length >= 8 and args.length <= 13)
assert(os.path.exists(args.netMHC))
assert(os.path.exists(args.fasta))
assert(os.path.exists(args.chains))

chainCollection = None
with open(args.chains, 'rb') as f:
    chainCollection = pickle.load(f)

scoreTableFile = None
scoreTable = None
if os.path.isfile(args.scoreTable):
    scoreTableFile = open(args.scoreTable, 'rb+')
    scoreTable = ScoreTable.readExisting(scoreTableFile)
else:
    scoreTableFile = open(args.scoreTable, 'wb+')
    scoreTable = ScoreTable.empty(scoreTableFile, 'H', 1, ' ')
scorer = NetMHCScorer(args.netMHC, 5000)

def getPeptideGen(chainCollection, fastaPath, pepLen):
    return map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, fastaPath, pepLen))
try:
    scoreTable.addAllele(scorer, args.allele, getPeptideGen(chainCollection, args.fasta, args.length))
except NetMHCRunFailedError as e:
    for x in e.runs.runs:
        print(x.success)
    print(e.message)
    sys.exit(1)



import io
import argparse
import os
import pickle
import sys
import functools
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable
from netMHCCalling import NetMHCScorer, NetMHCRunFailedError


parser = argparse.ArgumentParser(description='Precompute NetMHCPan scores for peptides')
parser.add_argument('path', help='Path to NetMHC executable')
parser.add_argument('fasta')
parser.add_argument('chains')
parser.add_argument('BAScoreTable')
parser.add_argument('ELScoreTable')

parser.add_argument('allele')
parser.add_argument('length', type=int)
parser.add_argument('threads', type=int)

args = parser.parse_args()
assert(args.length >= 8 and args.length <= 13)
assert(os.path.isfile(args.path))
assert(os.path.isfile(args.fasta))
assert(os.path.isfile(args.chains))
assert(args.threads >= 1)

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

def generateNetMHCCommand(netmhcPath, allele, inputFilePath, outputFilePath):
    return [netmhcPath, '-a', allele, '-p', '-xls', '-xlsfile', outputFilePath, '-f', inputFilePath, '-BA']

commandGen = functools.partial(generateNetMHCCommand, args.path, args.allele)
scorer = NetMHCScorer(5000, commandGen, args.threads)
scoreIter = scorer.scorePeptides(getPeptideGen(chainCollection, args.fasta, args.length), ['nM'])

def getPeptideGen(chainCollection, fastaPath, pepLen):
    return map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, fastaPath, pepLen))
try:
    scoreTable.addAllele(args.allele, scoreIter)
except NetMHCRunFailedError as e:
    for x in e.runs.runs:
        print(x.success)
    print(e.message)
    sys.exit(1)



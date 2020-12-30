import io
import argparse
import os
import pickle
import sys
import functools
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable, ScoreCategory, fileMD5
from netMHCCalling import NetMHCScorer, NetMHCRunFailedError


parser = argparse.ArgumentParser(description='Precompute NetMHC scores for peptides')
parser.add_argument('path', help='Path to NetMHC executable')
parser.add_argument('fasta')
parser.add_argument('chains')
parser.add_argument('scoreTable')

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

originalFASTAChecksum = chainCollection.md5
currentFASTAChecksum = fileMD5(args.fasta)
if currentFASTAChecksum != originalFASTAChecksum:
    print('MD5 of FASTA used to create chains file: ' + str(originalFASTAChecksum))
    print('MD5 of the given FASTA: '  + str(currentFASTAChecksum))
    print('These are not the same. Exiting')
    assert(False)

chainHash = chainCollection.checksum()

    
scoreTableFile = None
scoreTable = None
if os.path.isfile(args.scoreTable):
    scoreTableFile = open(args.scoreTable, 'rb+')
    scoreTable = ScoreTable.readExisting(scoreTableFile)
    if scoreTable.chainHash != chainHash:
        print('checksum for chains used to create score table: ' + str(scoreTable.chainHash))
        print('checksum for these chains: ' + str(chainHash))
        print('These are not the same. Exiting')
        assert(False)
else:
    scoreTableFile = open(args.scoreTable, 'wb+')
    scoreTable = ScoreTable.empty(scoreTableFile, ScoreCategory.AFFINITY, 'H', 1, ' ', args.length, chainHash)
    
assert(args.length == chainCollection.peptideLength)
assert(scoreTable.peptideLength == chainCollection.peptideLength)


if args.allele in scoreTable.getAlleles():
    print('Allele already in score table. Exiting')
    sys.exit(1)
def generateNetMHCCommand(netmhcPath, allele, inputFilePath):
    return [netmhcPath, '-a', allele, '-p', '-f', inputFilePath]
commandGen = functools.partial(generateNetMHCCommand, args.path, args.allele)
scorer = NetMHCScorer(5000, commandGen, args.threads)
def getPeptideGen(chainCollection, fastaPath, pepLen):
    return map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, fastaPath, pepLen))
scoreIter = map(lambda x: x[0], scorer.scorePeptides(getPeptideGen(chainCollection, args.fasta, args.length), ['Affinity(nM)']))
try:
    result = scoreTable.addAllele(args.allele, scoreIter)
except NetMHCRunFailedError as e:
    for x in e.runs.runs:
        print(x.success)
    print(e.message)
    sys.exit(1)



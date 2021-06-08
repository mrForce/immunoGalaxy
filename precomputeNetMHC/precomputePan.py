import io
import argparse
import os
import pickle
import sys
import functools
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable, ScoreCategory, fileMD5, ScoreTableGroup
from netMHCCalling import NetMHCScorer, NetMHCRunFailedError


parser = argparse.ArgumentParser(description='Precompute NetMHCPan scores for peptides')
parser.add_argument('path', help='Path to NetMHCPan executable')
parser.add_argument('fasta')
parser.add_argument('chains')
parser.add_argument('eluteScoreTable')
parser.add_argument('bindingScoreTable')

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

assert(args.length == chainCollection.peptideLength)

chainHash = chainCollection.checksum()

    
eluteScoreTableFile = None
eluteScoreTable = None
if os.path.isfile(args.eluteScoreTable):
    eluteScoreTableFile = open(args.eluteScoreTable, 'rb+')
    eluteScoreTable = ScoreTable.readExisting(eluteScoreTableFile)
    if eluteScoreTable.chainHash != chainHash:
        print('checksum for chains used to create elute score table: ' + eluteScoreTable.chainHash.hex())
        print('checksum for these chains: ' + chainHash.hex())
        print('These are not the same. Exiting')
        assert(False)
else:
    eluteScoreTableFile = open(args.eluteScoreTable, 'wb+')
    eluteScoreTable = ScoreTable.empty(eluteScoreTableFile, ScoreCategory.PAN_SCORE_EL, 'f', 1, ' ', args.length, chainHash)
    

assert(eluteScoreTable.peptideLength == chainCollection.peptideLength)

bindingScoreTableFile = None
bindingScoreTable = None
if os.path.isfile(args.bindingScoreTable):
    bindingScoreTableFile = open(args.bindingScoreTable, 'rb+')
    bindingScoreTable = ScoreTable.readExisting(bindingScoreTableFile)
    if bindingScoreTable.chainHash != chainHash:
        print('checksum for chains used to create binding score table: ' + bindingScoreTable.chainHash.hex())
        print('checksum for these chains: ' + chainHash.hex())
        print('These are not the same. Exiting')
        assert(False)
else:
    bindingScoreTableFile = open(args.bindingScoreTable, 'wb+')
    bindingScoreTable = ScoreTable.empty(bindingScoreTableFile, ScoreCategory.PAN_AFF, 'H', 1, ' ', args.length, chainHash)
    

assert(bindingScoreTable.peptideLength == chainCollection.peptideLength)



if args.allele in eluteScoreTable.getAlleles():
    print('Allele already in elute score table. Exiting')
    sys.exit(1)
if args.allele in bindingScoreTable.getAlleles():
    print('Allele already in binding score table. Exiting')
    sys.exit(1)
def generateNetMHCCommand(netmhcPath, allele, inputFilePath):
    return [netmhcPath, '-a', allele, '-p', '-f', inputFilePath, '-BA']
commandGen = functools.partial(generateNetMHCCommand, args.path, args.allele)
scorer = NetMHCScorer(5000, commandGen, args.threads)
def getPeptideGen(chainCollection, fastaPath, pepLen):
    return map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, fastaPath, pepLen))

#gives back an iterator
"""
This was the problem: how do we split up the iterators?
"""
scoreIter = scorer.scorePeptides(getPeptideGen(chainCollection, args.fasta, args.length), ['Score_EL', 'Aff(nM)'])
scoreTableGroup = ScoreTableGroup(eluteScoreTable, bindingScoreTable)

try:
    result = scoreTableGroup.addAllele(args.allele, scoreIter)
except NetMHCRunFailedError as e:
    for x in e.runs.runs:
        print(x.success)
    print(e.message)
    sys.exit(1)



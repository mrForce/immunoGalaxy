import io
import argparse
import os
import pickle
import numpy as np
import sys
import tempfile
import functools
import subprocess
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable, ScoreCategory, fileMD5, ScoreTableGroup

SCRIPT_LOCATION='/galaxy-prod/galaxy/tools-dependencies/bin/run_MHCFlurry_on_cd9.sh'
parser = argparse.ArgumentParser(description='Precompute MHCFlurry scores for peptides on remote server')
parser.add_argument('fasta')
parser.add_argument('chains')
parser.add_argument('bindingScoreTable')

parser.add_argument('allele')
parser.add_argument('length', type=int)

args = parser.parse_args()


assert(args.length >= 8 and args.length <= 13)
assert(os.path.isfile(args.fasta))
assert(os.path.isfile(args.chains))


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




if args.allele in bindingScoreTable.getAlleles():
    print('Allele already in binding score table. Exiting')
    sys.exit(1)

def getPeptideGen(chainCollection, fastaPath, pepLen):
    return map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, fastaPath, pepLen))



KNOWN_AMINOS=set('ACDEFGHIKLMNPQRSTVWY')
def unknownSub(peptides):
    pep = []
    for x in peptides:
        sub = ''.join([y if y in KNOWN_AMINOS else 'X' for y in x])
        pep.append(sub)
    return pep
peptides = [x.getPeptideSequence() for x in peptideGenerator(chainCollection, args.fasta, args.length)]
peptidesWithUnknown = unknownSub(peptides)

peptidesFile = tempfile.NamedTemporaryFile(mode='wt')
for x in peptidesWithUnknown:
    peptidesFile.write(x + '\n')
peptidesFile.flush()
scoresFile = tempfile.NamedTemporaryFile(mode='rb')
command = [SCRIPT_LOCATION, peptidesFile.name, scoresFile.name, args.allele]
proc = subprocess.Popen(command, stdout=subprocess.DEVNULL)
outs, errors = proc.communicate()
scores = list(np.load(scoresFile))


result = bindingScoreTable.addAllele(args.allele, iter(scores))



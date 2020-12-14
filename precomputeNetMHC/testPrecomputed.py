"""from precomputedNetMHCIndex import ChainCollection, NoDuplicatePeptideIterator, ScoreTable, DummyScorer"""
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable
from netMHCCalling import NetMHCScorer, NetMHCRunFailedError
import pickle
import sys
import itertools
import os
fastaPath = 'testBigger.fasta'
pepLen = 8
chains = ChainCollection(fastaPath, pepLen)
for chain in chains:
    print(chain)
with open('chains.pickle', 'wb') as f:
    pickle.dump(chains, f)

chainCollection = None
with open('chains.pickle', 'rb') as f:
    chainCollection = pickle.load(f)
for chain in chains:
    print(chain)
"""
pepIter = NoDuplicatePeptideIterator(chainCollection, 'test.fasta', 3)
while True:
    peptide = pepIter.getPeptideAndAdvance()
    if peptide:
        print(peptide)
    else:
        break
"""

"""pepIter = peptideGenerator(chainCollection, fastaPath, pepLen)

for pep in pepIter:
    print('>' + ' @ '.join(pep.getHeaders()))
    print(pep.getPeptideSequence())
"""
#os.remove('scores.bin')
scoreTable = ScoreTable('scores.bin')
"""first_scorer = DummyScorer(0)
second_scorer = DummyScorer(1)
third_scorer = DummyScorer(2)
scoreTable.addAllele(first_scorer, 'A',  map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, 'test.fasta', 3)))
scoreTable.addAllele(second_scorer, 'B', map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, 'test.fasta', 3)))
scoreTable.addAllele(third_scorer, 'C', map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, 'test.fasta', 3)))

for peptide, row in itertools.zip_longest(map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, 'test.fasta', 3)), scoreTable):
    print('peptide: ' + peptide)
    print('row')
    print(row) 
"""
scorer = NetMHCScorer('/home/jforce/Downloads/netMHC-4.0a.Linux/netMHC-4.0/netMHC', 10)
#scorer = NetMHCScorer('./dummyNetMHC.py', 10)
def getPeptideGen(chainCollection, fastaPath, pepLen):
    return map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, fastaPath, pepLen))
try:
    scoreTable.addAllele(scorer, 'HLA-A0101', getPeptideGen(chainCollection, fastaPath, pepLen))
except NetMHCRunFailedError as e:
    for x in e.runs.runs:
        print(x.success)
    print(e.message)
    sys.exit(1)
try:
    scoreTable.addAllele(scorer, 'HLA-A0201', getPeptideGen(chainCollection, fastaPath, pepLen))
except NetMHCRunFailedError as e:
    print(e.message)
    sys.exit(1)
num_peptides = 0
for peptide, row in itertools.zip_longest(getPeptideGen(chainCollection, fastaPath, pepLen), scoreTable):
    print('peptide: ' + peptide)
    print('row')
    print(row)
    num_peptides += 1
print('num peptides')
print(num_peptides)


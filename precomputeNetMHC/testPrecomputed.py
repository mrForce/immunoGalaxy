"""from precomputedNetMHCIndex import ChainCollection, NoDuplicatePeptideIterator, ScoreTable, DummyScorer"""
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable, DummyScorer, NetMHCScorer
import pickle
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
pepIter = peptideGenerator(chainCollection, fastaPath, pepLen)

for pep in pepIter:
    print('>' + ' @ '.join(pep.getHeaders()))
    print(pep.getPeptideSequence())
pepIter = map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, fastaPath, pepLen))
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
scoreTable.addAllele(scorer, 'HLA-A0101', map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, fastaPath, pepLen)))
scoreTable.addAllele(scorer, 'HLA-A0201', map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, fastaPath, pepLen)))

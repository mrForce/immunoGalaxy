"""from precomputedNetMHCIndex import ChainCollection, NoDuplicatePeptideIterator, ScoreTable, DummyScorer"""
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable, DummyScorer
import pickle
import itertools
import os
chains = ChainCollection('test.fasta', 3)
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
pepIter = peptideGenerator(chainCollection, 'test.fasta', 3)

for pep in pepIter:
    print('>' + ' @ '.join(pep.getHeaders()))
    print(pep.getPeptideSequence())
pepIter = map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, 'test.fasta', 3))
os.remove('scores.bin')
scoreTable = ScoreTable('scores.bin')
first_scorer = DummyScorer(0)
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
while True:
    pepAndHeader = pepIter.getPeptideWithHeadersAndAdvance()
    if pepAndHeader:
        peptide = pepAndHeader['peptide']
        header = ' @ '.join(pepAndHeader['headers'])
        print('>' + header)
        print(peptide)
    else:
        break
"""

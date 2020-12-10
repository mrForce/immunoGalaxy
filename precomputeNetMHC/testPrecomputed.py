"""from precomputedNetMHCIndex import ChainCollection, NoDuplicatePeptideIterator, ScoreTable, DummyScorer"""
from precomputedNetMHCIndex import ChainCollection, NoDuplicatePeptideIterator
import pickle
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
pepIter = NoDuplicatePeptideIterator(chainCollection, 'test.fasta', 3)
while True:
    pepAndHeader = pepIter.getPeptideWithHeadersAndAdvance()
    if pepAndHeader:
        peptide = pepAndHeader['peptide']
        header = ' @ '.join(pepAndHeader['headers'])
        print('>' + header)
        print(peptide)
    else:
        break


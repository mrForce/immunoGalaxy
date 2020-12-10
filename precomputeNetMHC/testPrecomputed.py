"""from precomputedNetMHCIndex import ChainCollection, NoDuplicatePeptideIterator, ScoreTable, DummyScorer"""
from precomputedNetMHCIndex import ChainCollection, NoDuplicatePeptideIterator
import pickle
chains = ChainCollection('test.fasta', 3)
for chain in chains:
    print(chain)
with open('chains.pickle', 'wb') as f:
    pickle.dump(chains, f)

with open('chains.pickle', 'rb') as f:
    chainCollection = pickle.load(f)
    for chain in chains:
        print(chain)


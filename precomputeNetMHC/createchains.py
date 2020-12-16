import argparse
from precomputedNetMHCIndex import ChainCollection
import pickle
import itertools
import os

parser = argparse.ArgumentParser()
parser.add_argument('fasta', help='path to FASTA file')
parser.add_argument('length', help='Peptide length', type = int)
parser.add_argument('output', help='Location to store chain collection pickle')

args = parser.parse_args()
assert(args.length >= 8)
assert(os.path.isfile(args.fasta))
assert(not os.path.isfile(args.output))

chains = ChainCollection(args.fasta, args.length)
with open(args.output, 'wb') as f:
    pickle.dump(chains, f)

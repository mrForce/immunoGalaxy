import argparse
import pickle
from Bio import SeqIO
from precomputedNetMHCIndex import peptideGenerator

parser = argparse.ArgumentParser(description='Extract peptides from chain file')
parser.add_argument('chains', help='path to chain file')
parser.add_argument('fasta', help='path to FASTA file')
parser.add_argument('output', help='Where to store peptides')

args=  parser.parse_args()

chainCollection = None
with open(args.chains, 'rb') as f:
    chainCollection = pickle.load(f)

pepGen = peptideGenerator(chainCollection, args.fasta, chainCollection.peptideLength)

g = open(args.output, 'w')
for peptide in pepGen:
    seq = peptide.getPeptideSequence()
    g.write(seq + '\n')
g.close()

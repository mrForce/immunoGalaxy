import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Extract peptides from FASTA file')
parser.add_argument('fasta', help='path to FASTA file')
parser.add_argument('length', help='Peptide length', type = int)
parser.add_argument('output', help='Where to store peptides')

args=  parser.parse_args()
f = open(args.fasta, 'r')
recordIterator = SeqIO.FastaIO.SimpleFastaParser(f)
g = open(args.output, 'w')

for header, sequence in recordIterator:
    pep = []
    for i in range(0, len(sequence) - args.length + 1):
        pep.append(sequence[i:(i + args.length)])
    if pep:
        pep.append('')
        g.write('\n'.join(pep))
f.close()
g.close()

from Bio import SeqIO
from precomputedNetMHCIndex import *
with open('test.fasta', 'r') as f:
    iterator = SeqIO.parse(f, 'fasta')
    occurences, proteinLengths = getPeptideOccurencesAndProteinLengths(iterator, 3)
    print('occurences')
    print(occurences)
    chains = createDuplicateChains(occurences, proteinLengths)
    sortDuplicateChains(chains)
    for i in range(0, len(chains)):
        print('protein: ' + str(i))
        for j in range(0, len(chains[i].chain)):
            print(chains[i].chain[j].__dict__)

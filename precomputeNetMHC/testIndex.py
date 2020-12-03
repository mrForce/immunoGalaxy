from Bio import SeqIO
from precomputedNetMHCIndex import *
with open('test.fasta', 'r') as f:
    iterator = SeqIO.parse(f, 'fasta')
    occurences, numProteins = getPeptideOccurences(iterator, 3)
    print('occurences')
    print(occurences)
    chains = createDuplicateChains(occurences, numProteins)
    sortDuplicateChains(chains)
    for i in range(0, len(chains)):
        print('protein: ' + str(i))
        for j in range(0, len(chains[i])):
            print(chains[i][j].__dict__)

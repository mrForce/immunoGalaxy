from Bio import SeqIO
from precomputedNetMHCIndex import duplicatePeptideChains

with open('test.fasta', 'r') as f:
    iterator = SeqIO.parse(f, 'fasta')
    table = duplicatePeptideChains(iterator, 3)
    for i in range(0, len(table)):
        print('protein: ' + str(i))
        for j in range(0, len(table[i])):
            print(table[i][j].__dict__)

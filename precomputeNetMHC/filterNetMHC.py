import io
import argparse
import os
import collections
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable

def getPeptideGen(chainCollection, fastaPath, pepLen):
    return map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, fastaPath, pepLen))

def scoreDistribution(table, allele, peptideGen, peptidesToExclude):
    #peptidesToExclude must be a set
    allelePosition = table.getAlleleList().index(allele)
    count = collections.Counter()
    for peptide, row in itertools.zip_longest(peptideGen, scoreTable):
        if peptide not in peptidesToExclude:
            count[row[allelePosition]] += 1
    return count

def filterNetMHC(allele, length, baseScoreTable, baseChainCollection, baseFasta, additionalScoreTable, additionalChainCollection, additionalFasta, k):
    baseScoreDist = collections.Count()    
    if baseScoreTable and baseChainCollection and baseFasta:
        exclusion = set()
        if additionalChainCollection:
            additionalGen = getPeptideGen(additionalChainCollection, additionalFasta, length)
            for peptide in additionalGen:                
                exclusion.add(peptide)
        baseGen = getPeptideGen(baseChainCollection, baseFasta, length)
        baseScoreDist = scoreDistribution(baseScoreTable, 
    pass
        

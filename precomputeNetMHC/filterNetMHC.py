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

def computeScoreThreshold(distribution, k):
    numPeptides = sum(distribution.values())
    maxVal = max(distribution.keys())
    numIncluded = 0
    for score,count in distribution:
        numIncluded += count
        #return the score for the bucket that pushes us over the edge
        if numIncluded >= k*numPeptides:
            return score
    return maxVal

def filterOnThreshold(table, allele, peptideGen, threshold):
    allelePosition = table.getAlleleList().index(allele)
    peptides = set()
    for peptide, row in itertools.zip_longest(peptideGen, scoreTable):
        if row[allelePosition] <= threshold:
            peptides.add(peptide)
    return peptides
    

def filterNetMHC(allele, length, baseScoreTable, baseChainCollection, baseFasta, additionalScoreTable, additionalChainCollection, additionalFasta, k):
    baseScoreDist = collections.Counter()    
    if baseScoreTable and baseChainCollection and baseFasta:
        exclusion = set()
        if additionalChainCollection:
            additionalGen = getPeptideGen(additionalChainCollection, additionalFasta, length)
            for peptide in additionalGen:                
                exclusion.add(peptide)
        baseGen = getPeptideGen(baseChainCollection, baseFasta, length)
        baseScoreDist = scoreDistribution(baseScoreTable, allele, baseGen, exclusion)
    additionalScoreDist = collections.Counter()
    if additionalScoreTable and additionalChainCollection and additionalFasta:
        additionalGen = getPeptideGen(additionalChainCollection, additionalFasta, length)
        additionalScoreDist = scoreDistribution(additionalScoreTable, allele, additionalGen, set())
    combinedScoreDist = additionalScoreDist + baseScoreDist
    threshold = computeScoreThreshold(combinedScoreDist, k)
    peptides = set()
    if baseScoreTable and baseChainCollection and baseFasta:
        peptides += filterOnThreshold(baseScoreTable, allele, getPeptideGen(baseChainCollection, baseFasta, length), threshold)
    if additionalScoreTable and additionalChainCollection and additionalFasta:
        peptides += filterOnThreshold(additionalScoreTable, allele, getPeptideGen(additionalChainCollection, additionalFasta, length), threshold)
    return peptides

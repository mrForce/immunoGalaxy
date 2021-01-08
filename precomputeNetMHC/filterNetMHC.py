import io
import itertools
import argparse
import os
import collections
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable

def getPeptideGen(chainCollection, fastaPath, pepLen):
    return map(lambda x: x.getPeptideSequence(), peptideGenerator(chainCollection, fastaPath, pepLen))

def scoreDistribution(scoreIter, peptideGen, peptidesToExclude):
    if peptidesToExclude:
        count = collections.Counter()
        for peptide, score in itertools.zip_longest(peptideGen, scoreIter):
            if peptide not in peptidesToExclude:
                count[score] += 1    
        return count
    else:
        return collections.Counter(scoreIter)


def computeScoreThreshold(distribution, k):
    numPeptides = sum(distribution.values())
    maxVal = max(distribution.keys())
    numIncluded = 0
    for score,count in distribution.items():
        numIncluded += count
        #return the score for the bucket that pushes us over the edge
        if numIncluded >= k*numPeptides:
            return score
    return maxVal

def extractPeptidesAndScores(scoreIter, peptideHolderGen, indices):
    i = 0
    j = 0
    for peptideHolder, score in itertools.zip_longest(peptideHolderGen, scoreIter):
        if indices[j] == i:
            yield (peptideHolder.getPeptideSequence(), score)
            j += 1
            if j == len(indices):
                break
        i += 1
            

def filterOnThreshold(scoreIter, peptideHolderGen, threshold):
    peptides = []
    for peptideHolder, score in itertools.zip_longest(peptideHolderGen, scoreIter):
        if score <= threshold:
            peptides.append(peptideHolder)
    return peptides

def mapPeptideToHeaders(peptideHolders):
    return {x.getPeptideSequence():x.getHeaders() for x in peptideHolders}


def filterNetMHC(allele, length, baseScoreTable, baseChainCollection, baseFasta, additionalScoreTable, additionalChainCollection, additionalFasta, k):
    baseScoreDist = collections.Counter()    
    if baseScoreTable and baseChainCollection and baseFasta:
        exclusion = None
        if additionalChainCollection:
            additionalGen = getPeptideGen(additionalChainCollection, additionalFasta, length)
            exclusion = set(additionalGen)
        baseGen = getPeptideGen(baseChainCollection, baseFasta, length)
        baseScoreDist = scoreDistribution(baseScoreTable.scoreIter(allele), baseGen, exclusion)
    additionalScoreDist = collections.Counter()
    if additionalScoreTable and additionalChainCollection and additionalFasta:
        additionalGen = getPeptideGen(additionalChainCollection, additionalFasta, length)
        additionalScoreDist = scoreDistribution(additionalScoreTable.scoreIter(allele), additionalGen, None)
    combinedScoreDist = additionalScoreDist + baseScoreDist
    threshold = computeScoreThreshold(combinedScoreDist, k)
    peptides = []
    pepToHeader = collections.defaultdict(list)
    if baseScoreTable and baseChainCollection and baseFasta:
        gen = peptideGenerator(baseChainCollection, baseFasta, length)
        pep = filterOnThreshold(baseScoreTable.scoreIter(allele), gen, threshold)
        for k, v in pep.items():
            pepToHeader[k].append(v)
    if additionalScoreTable and additionalChainCollection and additionalFasta:
        gen = peptideGenerator(additionalChainCollection, additionalFasta, length)        
        pep = filterOnThreshold(additionalScoreTable.scoreIter(allele), gen, threshold)
        for k,v in pep.items():
            pepToHeader[k].append(v)
    
    return pepToHeader

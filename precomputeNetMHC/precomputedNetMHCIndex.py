import io
import array
from Bio import SeqIO
from collections import Counter, namedtuple, defaultdict

class DuplicatePeptideMarker:
    def __init__(self, sequenceStart, nextProteinIndex, nextChainPosition, lastProteinIndex, lastChainPosition):
        self.sequenceStart = sequenceStart
        self.nextProteinIndex = nextProteinIndex
        self.nextChainPosition = nextChainPosition
        self.lastProteinIndex = lastProteinIndex
        self.lastChainPosition = lastChainPosition


class ProteinChain:
    def __init__(self, length):
        self.length = length
        self.chain = []
        
    
class TempPeptideTracker:
    def __init__(self, first, chainPosition, proteinIndex, sequenceStart):
        self.first = first
        self.chainPosition = chainPosition
        self.proteinIndex = proteinIndex
        self.sequenceStart = sequenceStart

class FASTAPeptideCollection:
    def __init__(self, fastaPeptides):
        self.fastaPeptides = fastaPeptides
    def writeToFile(self, delim, fastaPath):
        pass


    
def getPeptideOccurencesAndProteinLengths(recordIterator, peptideLength):
    peptides = defaultdict(list)
    proteinIndex = 0
    proteinLengths = []
    for record in recordIterator:
        for startIndex in range(0, len(record.seq) - peptideLength + 1):
            peptide = str(record.seq[startIndex:(startIndex + peptideLength)])
            peptides[peptide].append((proteinIndex, startIndex))
        proteinLengths.append(len(record.seq))
        proteinIndex += 1
        
    return (peptides, proteinLengths)

def createDuplicateChains(peptideOccurences, proteinLengths):
    chains = [ProteinChain(x) for x in proteinLengths]
    for peptide, occurences in dict(peptideOccurences).items():
        if len(occurences) > 1:
            print('occurences')
            print(occurences)
            lastProteinIndex = None
            lastChainPosition = None
            for i in range(0, len(occurences)):
                element = occurences[i]
                print('element')
                print(element)
                chains[element[0]].chain.append(DuplicatePeptideMarker(element[1], None, None, lastProteinIndex, lastChainPosition))
                print(chains[element[0]].chain)
                if i > 0:
                    assert(lastProteinIndex != None)
                    assert(lastChainPosition != None)
                    last = chains[lastProteinIndex].chain[lastChainPosition]
                    last.nextProteinIndex = element[0]
                    last.nextChainPosition = len(chains[element[0]].chain) - 1
                lastProteinIndex = element[0]
                lastChainPosition = len(chains[element[0]].chain) - 1
    return chains
def sortDuplicateChains(chains):
    for i in range(0, len(chains)):
        chains[i].chain.sort(key=lambda x: x.sequenceStart)
        for j in range(0, len(chains[i].chain)):
            element = chains[i].chain[j]
            if element.nextProteinIndex != None:
                chains[element.nextProteinIndex].chain[element.nextChainPosition].lastChainPosition = j
            if element.lastProteinIndex != None:
                chains[element.lastProteinIndex].chain[element.lastChainPosition].nextChainPosition = j



def duplicatePeptideChains(recordIterator, peptideLength):
    peptides = {}
    proteinIndex = 0
    chains = []
    for record in recordIterator:
        chains.append([])
        print('new protein')
        print(record.seq)
        for startIndex in range(0, len(record.seq) - peptideLength + 1):
            print('start')
            print(startIndex)
            peptide = str(record.seq[startIndex:(startIndex + peptideLength)])
            assert(len(peptide) == peptideLength)
            if peptide in peptides:
                tempTracker = peptides[peptide]
                if tempTracker.first:
                    #this is the second time we've seen this peptide
                    chains[tempTracker.proteinIndex].append(DuplicatePeptideMarker(tempTracker.sequenceStart, None, None, None, None))
                    tempTracker.chainPosition = len(chains[tempTracker.proteinIndex]) - 1
                chains[proteinIndex].append(DuplicatePeptideMarker(startIndex, None, None, tempTracker.proteinIndex, tempTracker.chainPosition))
                t = chains[tempTracker.proteinIndex][tempTracker.chainPosition]
                t.nextProteinIndex = proteinIndex
                t.nextChainPosition = len(chains[proteinIndex]) - 1
                peptides[peptide] = TempPeptideTracker(False, t.nextChainPosition, proteinIndex, startIndex)
            else:
                peptides[peptide] = TempPeptideTracker(True, None, proteinIndex, startIndex)
        proteinIndex += 1

                    
class ProteinPeptideTable:
    def __init__(self, fastaFile, peptideLength):
        self.duplicatePeptideMapper = []
            

    def peptideIterator(self):
        pass
    def peptidePositionIterator(self):
        pass
    def FASTAPeptideIterator(self, positions):
        pass
    def getNumPeptides(self):
        pass
class ScoreTable:
    #ppt is an instance of the ProteinPeptideTable class
    def __init__(self, ppt):
        self.ppt = ppt
        self.alleleScores = {}
    def addAllele(self, netmhcCaller, alleleName):
        assert(alleleName not in self.alleleScores)
        peptideIter = self.ppt.peptideIterator()
        scoreIter = netmhcCaller.scorePeptides(alleleName, peptideIter)
        scoreArray = array.array('I', scoreIter)
        assert(len(scoreArray) == self.ppt.getNumPeptides())
        self.alleleScores[alleleName] = scoreArray
    
    def getTopFraction(self, alleleName, fraction):
        #fraction should be 0.02 for 2%
        assert(alleleName in self.alleleScores)
        scores = self.alleleScores[alleleName]
        c = Counter(scores)
        numPeptides = len(scores)
        numPeptidesToInclude = int(numPeptides*fraction)
        threshold = -1
        numPeptidesIncluded = 0
        for score,count in c.items():
            threshold = score
            numPeptidesIncluded += count
            if numPeptidesIncluded >= numPeptidesToInclude:
                break
        assert(threshold > -1)
        scoreIndex = 0
        peptidePositionsToKeep = []
        for peptidePosition in self.ppt.peptidePositionIterator():
            if scores[scoreIndex] <= threshold:
                peptidePositionsToKeep.append(peptidePosition)
        return peptidePositionsToKeep    

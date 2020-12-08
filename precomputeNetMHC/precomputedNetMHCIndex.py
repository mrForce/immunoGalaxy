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
    def __init__(self, proteinLength):
        self.proteinLength = proteinLength
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

class NoDuplicatePeptideIterator:
    def __init__(self, chains, recordIterator, pepLen):
        self.chains = chains
        self.pepLen = pepLen
        self.recordIterator = recordIterator
        self.currentRecord = next(self.recordIterator)
        self.chainIndex = 0
        self.positionInProtein = 0
        self.positionInChain = 0
    def moveForward(self):
        link = None
        if self.positionInChain < len(self.chains[self.chainIndex].chain):
            currentChainLink = self.chains[self.chainIndex].chain[self.positionInChain]
            if currentChainLink.sequenceStart == self.positionInProtein:
                link = currentChainLink
                self.positionInChain += 1
        self.positionInProtein += 1
        return link
    def moveToNextProtein(self):
        self.chainIndex += 1
        self.positionInProtein = 0
        self.positionInChain = 0
        self.currentRecord = next(self.recordIterator)

    
    def moveToValid(self):
        """
        Returns False if we go past the chain collection. Returns True otherwise. 
        """
        if self.chainIndex >= len(self.chains):
            return False
        if self.positionInProtein > self.chains[self.chainIndex].proteinLength - self.pepLen:
            try:
                self.moveToNextProtein()
            except StopIteration:
                return False
            return self.moveToValid()
        else:
            #finish here
            link = self.moveForward()
            if link and link.lastProteinIndex != None:
                
                
            if self.positionInChain < len(self.chains[self.chainIndex]):
                if self.chains[self.chainIndex][self.positionInChain].sequenceStart == self.positionInProtein and self.chains[self.chainIndex][self.positionInChain].lastProteinIndex != None:
                    self.positionInProtein += 1
                    self.positionInChain += 1
                    return self.moveToValid()
            return True
            
    def skipPeptide(self):
        if self.positionInChain < len(self.chains[self.chainIndex]):
            if self.chains[self.chainIndex][self.positionInChain
        self.positionInProtein += 1
        self.advance()
    def getPeptide(self, headers=False):
        if 
        pass
        

                

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


class ScoreTable:
    def __init__(self, filename):
        self.filename = filename
        self.alleles = []
        #the number of bytes the list of alleles takes up. Basically, read this # of bytes, and you're add the start of the score table.
        self.headerSize = -1
        if os.path.isfile(filename):            
            reader = open(filename, 'rb')
            alleleListSize = int.from_bytes(reader.read(4), 'little')
            alleleListBytes = array.array('u', '')
            alleleListBytes.frombytes(reader.read(alleleListSize))
            alleleListString = alleleListBytes.tounicode()
            self.alleles = alleleListString.split(' ')
            reader.close()
            self.headerSize = 4 + alleleListSize
            z = array.array('I', [])
            self.rowSize = len(self.alleles)*z.itemsize
    def getTableReader(self):
        if self.headerSize < 0:
            return None
        else:
            reader = open(self.filename, 'rb')
            reader.read(self.headerSize)
            return reader
    def addAllele(self, netmhcCaller, alleleName, peptideIterator):
        assert(alleleName not in self.alleles)
        reader = self.getTableReader()
        self.alleles.append(alleleName)
        #create a temporary file to hold the output. 
        fd, path = temp.mkstemp()
        tempWriter = open(path, 'wb')
        alleleListString = ' '.join(self.alleles)
        alleleArray = array.array('u', alleleListString)
        alleleListBytes = alleleArray.tobytes()
        tempWriter.write(len(alleleListBytes).to_bytes(4, 'little'))
        tempWriter.write(alleleListBytes)
        scoreIter = netmhcCaller.scorePeptides(alleleName, peptideIter)
        z = array.array('I', [])        
        inputRowSize = z.itemsize*(len(self.alleles) - 1)
        outputRowSize = z.itemsize*len(self.alleles)
        for score in scoreIter:
            x = array.array('I', [])
            if inputRowSize > 0:
                inputRowBytes = self.reader.read(inputRowSize)
                x.extend(inputRowBytes)
            x.append(int(score))
            outputRowBytes = x.tobytes()
            assert(outputRowSize == len(outputRowBytes))
            tempWriter.write(outputRowBytes)
        tempWriter.flush()
      
        shutil.copyfile(path, self.filename)
        os.remove(path)
        self.__init__(self.filename)
    def getAlleleList(self):
        return self.alleles
    
    def __iter__(self):
        self.reader = self.getTableReader()
        assert(self.reader != None)
        return self
    def __next__(self):
        rowBytes = self.reader.read(self.rowSize)
        rowArray = array.array('I', rowBytes)
        return rowArray
        

                    
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
    
    def getScoreThreshold(self, alleleName, fraction):
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
        return threshold
    def __getstate__(self):
        
        pass
    def __setstate__(self):
        pass

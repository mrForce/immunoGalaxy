import io
import array
from Bio import SeqIO
from collections import Counter, namedtuple, defaultdict, UserList

class ChainLink:
    def __init__(self, sequenceStart, nextChainIndex, nextChainPosition, lastChainIndex, lastChainPosition):
        self.sequenceStart = sequenceStart
        self.nextChainIndex = nextChainIndex
        self.nextChainPosition = nextChainPosition
        self.lastChainIndex = lastChainIndex
        self.lastChainPosition = lastChainPosition


class Chain(UserList):
    def __init__(self, proteinLength, header):
        self.proteinLength = proteinLength
        self.header = header
        self.data = []

class ChainCollection(UserList):
    def __init__(self, fastaPath, pepLen):
        f = open(fastaPath, 'r')
        recordIterator = SeqIO.FastaIO.SimpleFastaParser(f)        
        self.data = []                    
        self.initChains(proteins)
        self.pruneChains()
    def initChains(self, recordIterator):
        peptides = {}
        protIndex = 0
        for header, sequence in recordIterator:
            chain = Chain(len(sequence), header)
            for i in range(0, len(sequence) - pepLen + 1):
                pep = sequence[i:(i + pepLen)]
                link = ChainLink(i, None, None, None, None)
                if pep in peptides:
                    chainIndex, chainPosition = peptides[pep]
                    self.data[chainIndex][chainPosition].nextProteinIndex = protIndex
                    self.data[chainIndex][chainPosition].nextChainPosition = i
                    link.lastProteinIndex = chainIndex
                    link.lastChainPosition = chainPosition
                peptides[pep] = (protIndex, i)                
            protIndex += 1
    def pruneChains(self):
        for i in range(0, len(self.data)):
            tempChain = Chain(self.data[i].proteinLength, self.data[i].header)
            for x in self.data[i]:
                if x.lastChainIndex != None or x.nextChainIndex != None:
                    if x.lastChainIndex != None:
                        self.data[x.lastChainIndex][x.lastChainPosition].nextChainPosition = len(tempChain)
                    if x.nextChainIndex != None:
                        self.data[x.nextChainIndex][x.nextChainPosition].lastChainPosition = len(tempChain)
                    tempChain.append(x)
            self.data[i] = tempChain
class NoDuplicatePeptideIterator:
    def __init__(self, chains, recordIterator, pepLen):
        self.chains = chains
        self.pepLen = pepLen
        self.recordIterator = recordIterator
        self.currentRecord = next(self.recordIterator)
        self.chainIndex = 0
        self.positionInProtein = 0
        self.positionInChain = 0
        self.moveToValid()
    def moveForward(self):
        if self.positionInChain < len(self.chains[self.chainIndex].chain):
            currentChainLink = self.chains[self.chainIndex].chain[self.positionInChain]
            if currentChainLink.sequenceStart == self.positionInProtein:
                self.positionInChain += 1
        self.positionInProtein += 1
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
            if self.positionInChain < len(self.chains[self.chainIndex].chain) and self.chains[self.chainIndex].chain[self.positionInChain].lastProteinIndex != None:
                self.moveForward()
                return self.moveToValid()
            return True
            
    def skipPeptide(self):
        self.moveForward()
        return self.moveToValid()
    def gatherHeaders(self):
        headers = []
        chainIndex = self.chainIndex
        chainPos = self.positionInChain
        if chainIndex >= len(self.chains) or chainPos >= len(self.chains[self.chainIndex].chain):
            return False
        while chainIndex != None and chainPos != None:
            headers.append(self.chains[chainIndex].header)
            tempChainPos = self.chains[chainIndex].chain[chainPos].nextChainPosition
            chainIndex = self.chains[chainIndex].chain[chainPos].nextProteinIndex
            chainPos = tempChainPos
        return headers
    def getPeptideAndAdvance(self):            
        if self.chainIndex >= len(self.chains) or self.positionInProtein > self.chains[self.chainIndex].proteinLength - self.pepLen:
            return False
        return self.currentRecord.seq[self.positionInProtein:(self.positionInProtein + self.pepLen)]
    def getPeptideWithHeadersAndAdvance(self):        
        headers = self.gatherHeaders()
        if headers:
            peptide = self.getPeptideAndAdvance()
            if peptide:            
                return {'peptide': peptide, 'headers': headers}
        return False
        

                


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

import io
import array
import os
import shutil
from Bio import SeqIO
import tempfile
from collections import Counter, namedtuple, defaultdict, UserList
import hashlib
from abc import ABC, abstractmethod


class ChainLink:
    def __init__(self, sequenceStart, nextChainIndex, nextChainPosition, lastChainIndex, lastChainPosition):
        self.sequenceStart = sequenceStart
        self.nextChainIndex = nextChainIndex
        self.nextChainPosition = nextChainPosition
        self.lastChainIndex = lastChainIndex
        self.lastChainPosition = lastChainPosition
    def __str__(self):
        return '({}, {}, {}, {}, {})'.format(self.sequenceStart, self.nextChainIndex, self.nextChainPosition, self.lastChainIndex, self.lastChainPosition)


class Chain(UserList):
    def __init__(self, proteinLength, header):
        self.proteinLength = proteinLength
        self.header = header
        self.data = []

    def __str__(self):
        return 'length: {}, header: {}, links: {}'.format(self.proteinLength, self.header, '->'.join([str(x) for x in self.data]))


def fileMD5(path):
    with open(path, 'rb') as f:
        fHash = hashlib.md5()
        chunk = f.read(8192)
        while chunk:
            fHash.update(chunk)
            chunk = f.read(8192)
        return fHash.hexdigest()
        
class ChainCollection(UserList):
    def __init__(self, fastaPath, pepLen):
        self.md5 = fileMD5(fastaPath)
        f = open(fastaPath, 'r')
        recordIterator = SeqIO.FastaIO.SimpleFastaParser(f)
        self.pepLen = pepLen
        self.data = []                    
        self.initChains(recordIterator)
        self.pruneChains()
        f.close()
    def __getstate__(self):
        state = {'md5': self.md5, 'data': self.data.copy(), 'pepLen': self.pepLen}
        return state
    def __setstate__(self, state):
        self.md5 = state['md5']
        self.data = state['data'].copy()
        self.pepLen = state['pepLen']
        
    def initChains(self, recordIterator):
        peptides = {}
        protIndex = 0
        for header, sequence in recordIterator:
            chain = Chain(len(sequence), header)
            self.data.append(chain)
            for i in range(0, len(sequence) - self.pepLen + 1):
                pep = sequence[i:(i + self.pepLen)]
                link = ChainLink(i, None, None, None, None)
                chain.append(link)
                if pep in peptides:
                    chainIndex, chainPosition = peptides[pep]
                    
                    self.data[chainIndex][chainPosition].nextChainIndex = protIndex
                    self.data[chainIndex][chainPosition].nextChainPosition = i
                    link.lastChainIndex = chainIndex
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
    def __init__(self, chains, fastaPath, pepLen):
        self.fastaPath = fastaPath
        self.openFASTA = open(fastaPath, 'r')
        self.recordIterator = SeqIO.FastaIO.SimpleFastaParser(self.openFASTA)
        self.chains = chains
        self.pepLen = pepLen
        _, self.currentRecord = next(self.recordIterator)
        self.chainIndex = 0
        self.positionInProtein = 0
        self.positionInChain = 0
        self.moveToValid()
    def terminate(self):
        self.openFASTA.close()
    
    def moveForward(self):
        if self.positionInChain < len(self.chains[self.chainIndex]):
            currentChainLink = self.chains[self.chainIndex][self.positionInChain]
            if currentChainLink.sequenceStart == self.positionInProtein:
                self.positionInChain += 1
        self.positionInProtein += 1
    def moveToNextProtein(self):
        self.chainIndex += 1
        self.positionInProtein = 0
        self.positionInChain = 0
        _, self.currentRecord = next(self.recordIterator)    
    def moveToValid(self):
        """
        Returns False if we go past the chain collection. Returns True otherwise. 
        """
        if self.chainIndex >= len(self.chains):
            self.terminate()
            return False
        if self.positionInProtein > self.chains[self.chainIndex].proteinLength - self.pepLen:
            try:
                self.moveToNextProtein()
            except StopIteration:
                self.terminate()
                return False
            return self.moveToValid()
        else:
            if self.positionInChain < len(self.chains[self.chainIndex]) and self.chains[self.chainIndex][self.positionInChain].sequenceStart == self.positionInProtein and self.chains[self.chainIndex][self.positionInChain].lastChainIndex != None:
                self.moveForward()
                return self.moveToValid()
            return True
            
    def skipPeptide(self):
        self.moveForward()
        return self.moveToValid()

    
    def gatherHeaders(self):
        headers = set()
        chainIndex = self.chainIndex
        chainPos = self.positionInChain
        
        if chainIndex >= len(self.chains):
            return False
        if chainPos >= len(self.chains[self.chainIndex]):
            return set([self.chains[chainIndex].header])
        while chainIndex != None and chainPos != None:            
            headers.add(self.chains[chainIndex].header)
            tempChainPos = self.chains[chainIndex][chainPos].nextChainPosition
            chainIndex = self.chains[chainIndex][chainPos].nextChainIndex
            chainPos = tempChainPos
        return headers
    def getPeptideAndAdvance(self):            
        if self.chainIndex >= len(self.chains) or self.positionInProtein > self.chains[self.chainIndex].proteinLength - self.pepLen:
            return False
        pep = self.currentRecord[self.positionInProtein:(self.positionInProtein + self.pepLen)]
        self.moveForward()
        self.moveToValid()
        return pep
    
    def getPeptideWithHeadersAndAdvance(self):        
        headers = self.gatherHeaders()
        if headers:
            peptide = self.getPeptideAndAdvance()
            if peptide:            
                return {'peptide': peptide, 'headers': headers}
        self.terminate()
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
        fd, path = tempfile.mkstemp()
        tempWriter = open(path, 'wb')
        alleleListString = ' '.join(self.alleles)
        alleleArray = array.array('u', alleleListString)
        alleleListBytes = alleleArray.tobytes()
        tempWriter.write(len(alleleListBytes).to_bytes(4, 'little'))
        tempWriter.write(alleleListBytes)
        scoreIter = netmhcCaller.scorePeptides(alleleName, peptideIterator)
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


class AbstractScorer(ABC):
    @abstractmethod
    def scorePeptides(self, alleleName, pepIter):
        pass

class DummyScorer(AbstractScorer):
    def scorePeptides(self, alleleName, pepIter):
        for x in pepIter:
            yield sum([ord(y) - ord('A') for y in x])

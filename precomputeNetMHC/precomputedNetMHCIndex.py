import io
import array
import os
import itertools
import mmap
from enum import Enum
import pdb
import csv
import subprocess
import struct
import shutil
import threading
from Bio import SeqIO
import tempfile
from collections import Counter, namedtuple, defaultdict, UserList
import queue
import hashlib
from abc import ABC, abstractmethod
MAX_SCORE = 50000

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

class PeptideHolder:
    def __init__(self, chains, proteinSequence, startPosition, chainIndex, chainPosition, pepLen):
        self.chains = chains
        self.proteinSequence = proteinSequence
        self.startPosition = startPosition
        self.chainIndex = chainIndex
        self.chainPosition = chainPosition
        self.pepLen = pepLen
    def getPeptideSequence(self):
        return self.proteinSequence[self.startPosition:(self.startPosition + self.pepLen)]
    def getHeaders(self):
        assert(self.chainIndex < len(self.chains))
        chainIndex = self.chainIndex
        chainPos = self.chainPosition
        headers = set([self.chains[chainIndex].header])
        while chainPos != None:
            headers.add(self.chains[chainIndex].header)
            tempChainPos = self.chains[chainIndex][chainPos].nextChainPosition
            chainIndex = self.chains[chainIndex][chainPos].nextChainIndex
            chainPos = tempChainPos
        return headers

def peptideGenerator(chains, fastaPath, pepLen):
    openFASTA = open(fastaPath, 'r')
    recordIter = SeqIO.FastaIO.SimpleFastaParser(openFASTA)
    chainIndex = 0
    chainPosition = 0
    proteinPosition = 0
    _, protein = next(recordIter)
    while chainIndex < len(chains):        
        if proteinPosition > chains[chainIndex].proteinLength - pepLen:
            chainIndex += 1
            proteinPosition = 0
            chainPosition = 0
            try:
                _, protein = next(recordIter)
            except StopIteration:
                break
        else:
            if chainPosition >= len(chains[chainIndex]) or chains[chainIndex][chainPosition].sequenceStart != proteinPosition:
                yield PeptideHolder(chains, protein, proteinPosition, chainIndex, None, pepLen)
            else:
                if chains[chainIndex][chainPosition].lastChainIndex == None:
                    yield PeptideHolder(chains, protein, proteinPosition, chainIndex, chainPosition, pepLen)                
                chainPosition += 1
            proteinPosition += 1
    openFASTA.close()

class ScoreCategory(Enum):
    PAN_SCORE_EL = 1
    PAN_RANK_EL = 2
    PAN_SCORE_BA = 3
    PAN_RANK_BA = 4
    PAN_AFF = 5
    #above were for NetMHCPan, below for regular NetMHC
    AFFINITY = 6
    RANK = 7

class ScoreTableMeta:
    """   
    In terms of bytes, this is structured as follows:

+------------+----------------+----------+--------------+--------------+--------------------------------------+---------+
| Multiplier | Score Category | Typecode | Num Peptides | Allele Delim | Size of following alleles (in bytes) | Alleles |
+------------+----------------+----------+--------------+--------------+--------------------------------------+---------+

The multiplier, score type, typecode, num Peptides, allele delim, and size of alleles are stored as a packed struct using Python's struct module, with the following format string: !HbQbL

Score category is one of the values for the ScoreCategory enum, indicating what kind of score it is. NetMHCPan gives 5 scores: elution score, elution rank, binding score, binding rank, predicted affinity. NetMHC gives 2 scores: Predicted affinity and rank.

Typecode is the typecode used in the array module for storing the NetMHC scores.

The list of allele names are stored as an ASCII string, with 'allele delim' seperating the different names. 

Note that both allele delim and typecode are stored as ascii characters (so a signed byte). But, pass them to the constructor as normal python strings. 
    """
    STRUCT_FORMAT='!HIcQcL'
    def __init__(self, multiplier, alleleSeperator, alleleList, numPeptides, scoreTypecode, scoreCategory):
        self.multiplier = multiplier        
        self.alleleSeperator = alleleSeperator.encode('ascii')
        assert(len(self.alleleSeperator) == 1)
        self.alleleList = [x.encode('ascii') for x in alleleList]
        self.numPeptides = numPeptides
        self.scoreTypecode = scoreTypecode.encode('ascii')
        assert(len(self.scoreTypecode) == 1)
        assert(isinstance(scoreCategory, ScoreCategory))
        self.scoreCategory = scoreCategory
    def getMultiplier(self):
        return self.multiplier
    def getAlleles(self):
        return [x.decode('ascii') for x in self.alleleList]
    def addAllele(self, allele):
        self.alleleList.append(allele.encode('ascii'))
    def getNumPeptides(self):
        return self.numPeptides
    def setNumPeptides(self, numPeptides):
        self.numPeptides = numPeptides
    def computeAlleleOffset(self, allele):
        print('list')
        print(self.alleleList)
        index = self.alleleList.index(allele.encode('ascii'))
        a = array.array(self.getScoreTypecode())
        return self.numPeptides*index*a.itemsize
    def getScoreTypecode(self):
        return self.scoreTypecode.decode('ascii')
    def getScoreCategory(self):
        return self.scoreCategory
    def computeScoreSectionSizeInBytes(self):
        a = array.array(self.getScoreTypecode())
        return self.numPeptides*len(self.alleleList)*a.itemsize
    @classmethod
    def fromBytes(cls, byteArray):
        print('byte array')
        print(byteArray)
        structSize = struct.calcsize(cls.STRUCT_FORMAT)
        multiplier, scoreCategoryInt, scoreTypecode, numPeptides, alleleSeperator, alleleArraySize = struct.unpack(cls.STRUCT_FORMAT, byteArray[0:structSize])        
        alleleArrayBytes = byteArray[structSize:(structSize + alleleArraySize)]
        alleleList = [x.decode('ascii') for x in alleleArrayBytes.split(alleleSeperator)]
        return cls(multiplier, alleleSeperator.decode('ascii'), alleleList, numPeptides, scoreTypecode.decode('ascii'), ScoreCategory(scoreCategoryInt))
    def toBytes(self):
        alleleArrayBytes = self.alleleSeperator.join(self.alleleList)
        metaHeader = struct.pack(self.STRUCT_FORMAT, self.multiplier, self.scoreCategory.value, self.scoreTypecode, self.numPeptides, self.alleleSeperator, len(alleleArrayBytes))
        return metaHeader + alleleArrayBytes


class ScoreTable:
    """
    The ScoreTable file should start with the MAGIC 4 bytes (4344382b in hex. It's just the ascii string 'CD8+' converted to a hex number). 
    """
    MAGIC=0X4344382b
    def __init__(self, fileObj, tableMeta):
        self.fileObj = fileObj
        self.tableMeta = tableMeta
    def getAlleles(self):
        return self.tableMeta.getAlleles()
    @classmethod
    def empty(cls, fileObj, scoreCategory, scoreTypecode, multiplier, alleleSeperator):
        meta = ScoreTableMeta(multiplier, alleleSeperator, [], 0, scoreTypecode, scoreCategory)
        return cls(fileObj, meta)
    @classmethod
    def readExisting(cls, fileObj):
        size = struct.calcsize('!IL')        
        magic, metaOffsetFromEnd = struct.unpack('!IL', fileObj.read(size))
        print('meta offset')
        print(metaOffsetFromEnd)
        assert(magic == cls.MAGIC)
        fileObj.seek(-1*metaOffsetFromEnd, os.SEEK_END)
        metaBytes = fileObj.read()
        meta = ScoreTableMeta.fromBytes(metaBytes)
        return cls(fileObj, meta)

    def addAllele(self, allele, scoreIter):
        if allele in self.tableMeta.getAlleles():
            return False
        else:
            location = struct.calcsize('!IL') + self.tableMeta.computeScoreSectionSizeInBytes()
            numScoresWritten = self.appendScores(scoreIter, location)
            if len(self.tableMeta.alleleList) == 0:
                self.tableMeta.setNumPeptides(numScoresWritten)
            if numScoresWritten == self.tableMeta.getNumPeptides():                    
                self.tableMeta.addAllele(allele)                
            self.writeMeta()
            self.fileObj.flush()
            return True
    def appendScores(self, scoreIter, location):
        chunkSize = 2**10
        self.fileObj.seek(location)
        numScores = 0
        while True:
            chunk = [int(x) for x in itertools.islice(scoreIter, chunkSize)]
            print('chunk')
            print(chunk)
            numScores += len(chunk)
            if len(chunk) == 0:
                break
            chunkArray = array.array(self.tableMeta.getScoreTypecode(), chunk)
            chunkBytes = chunkArray.tobytes()
            self.fileObj.write(chunkBytes)
        return numScores
    def writeMeta(self):
        beginningOffset = struct.calcsize('!IL') + self.tableMeta.computeScoreSectionSizeInBytes()
        self.fileObj.seek(beginningOffset)
        metaBytes = self.tableMeta.toBytes()
        self.fileObj.write(metaBytes)        
        self.fileObj.seek(0, os.SEEK_END)
        endLoc = self.fileObj.tell()
        self.fileObj.seek(0)
        print('end loc')
        print(endLoc)
        print('beginning offset')
        print(beginningOffset)
        self.fileObj.write(struct.pack('!IL', self.MAGIC, endLoc - beginningOffset))
        
    def scoreIter(self, allele):
        if allele in self.tableMeta.getAlleles():
            loc = struct.calcsize('!IL') + self.tableMeta.computeAlleleOffset(allele)
            self.fileObj.seek(loc)
            typecode = self.tableMeta.getScoreTypecode()            
            size = struct.calcsize(typecode)
            for i in range(0, self.tableMeta.getNumPeptides()):
                score,  = struct.unpack(typecode, self.fileObj.read(size))
                yield score

"""
def readScoreTableHeader(fileObj):
    size = struct.calcsize('!L')
    headerOffsetFromEnd, = struct.unpack('!L', fileObj.read(size))
    fileObj.seek(headerOffsetFromEnd, os.SEEK_END)
    headerBytes = fileObj.read()
    return ScoreTableHeader(headerBytes)
"""
""""
def setScoreTableHeader(fileObj, header, startLocation):
    headerBytes = header.toBytes()
    fileObj.seek(startLocation)
    
fileObj.write(struct.pack('!L', offset))
    fileObj.seek(0, os.SEEK_END)
    mmapObject.write(headerBytes)
    mmapObject.seek(currPos)
        


class ScoreTableReader:
    def __init__(self, filename):
        assert(os.path.isfile(filename))
        fileObject = open(filename, 'rb')
        self.memMap = mmap.mmap(fileObject.fileno(), 0)
        header = readScoreTableHeader(self.memMap)
        self.alleleList = header.getAlleleList()
        self.multiplier = header.getMultiplier()
        self.numPeptides = header.getNumPeptides()
    def scoreGenerator(self, allele):
        alleleLoc = self.alleleList.index(allele)
        alleleStartByteLocation = struct.calcsize('!L') + struct.calcsize('!H')*alleleLoc*self.numPeptides
        byteIter = itertools.islice(self.memMap, alleleStartByteLocation, alleleStartByteLocation + struct.calcsize('!H')*self.numPeptides)
        for score, in struct.iter_unpack('!H', byteIter):
            yield score
        
    def getAlleleList(self):
        return self.alleleList
    def getMultiplier(self):
        return self.multiplier
    def getNumPeptides(self):
        return self.numPeptides

def createEmptyScoreTable(filename, numPeptides, multiplier):
    assert(not os.path.isfile(filename))
    fileObject = open(filename, 'wb')
    memMap = mmap.mmap(fileObject.fileno(), 0)
    header = ScoreTableHeader(multiplier, '\t', [], numPeptides)
    setScoreTableHeader(memMap, header)
    memMap.flush()
    memMap.close()
    fileObject.close()


def addAlleleToScoreTable(fileObj, allele, netmhcCaller, peptideIterator):
    header = readScoreTableHeader(fileObj)
    assert(allele not in header.getAlleleList())
    
    with tempfile.TemporaryFile() as f:
        memMap = mmap.mmap(f.fileno())
    
"""
    
"""    
class ScoreTable:
    def __init__(self, filename):
        self.typecode = 'H'
        try:
            t = array.array(self.typecode, [50000])
        except OverflowError:
            self.typecode = 'I'
        t = array.array(self.typecode, [50000])
        print('typecode: ' + self.typecode)
        self.filename = filename
        self.alleles = []
        #the number of bytes the list of alleles takes up. Basically, read this # of bytes, and you're add the start of the score table.
        self.headerSize = -1
        if os.path.isfile(filename):            
            reader = open(filename, 'rb')
            alleleListSize = int.from_bytes(reader.read(4), 'little')
            print('allele list size: ' + str(alleleListSize))
            alleleListBytes = array.array('u', '')
            alleleListBytes.frombytes(reader.read(alleleListSize))
            alleleListString = alleleListBytes.tounicode()
            self.alleles = alleleListString.split(' ')
            reader.close()
            self.headerSize = 4 + alleleListSize
            z = array.array(self.typecode, [])
            self.rowSize = len(self.alleles)*z.itemsize
    def getTableReader(self):
        if self.headerSize < 0:
            return None
        else:
            reader = open(self.filename, 'rb')
            reader.read(self.headerSize)
            return reader
    def addAllele(self, netmhcCaller, alleleName, peptideIterator):
        print('starting to add allele')
        assert(alleleName not in self.alleles)
        reader = self.getTableReader()
        self.alleles.append(alleleName)
        #create a temporary file to hold the output. 
        fd, path = tempfile.mkstemp()
        tempWriter = open(path, 'wb')
        alleleListString = ' '.join(self.alleles)
        print('allele list string: ' + alleleListString)
        alleleArray = array.array('u', alleleListString)
        alleleListBytes = alleleArray.tobytes()
        print('allele list bytes')
        print(alleleListBytes)
        tempWriter.write(len(alleleListBytes).to_bytes(4, 'little'))
        tempWriter.write(alleleListBytes)
        scoreIter = netmhcCaller.scorePeptides(alleleName, peptideIterator)
        z = array.array(self.typecode, [])
        print(self.alleles)
        inputRowSize = z.itemsize*(len(self.alleles) - 1)
        print('input row size: ' + str(inputRowSize))
        outputRowSize = z.itemsize + inputRowSize
        print('output row size: ' + str(outputRowSize))
        for score in scoreIter:
            x = array.array(self.typecode, [])
            if inputRowSize > 0:
                inputRowBytes = array.array(self.typecode, reader.read(inputRowSize))
                print('input row bytes')
                print(inputRowBytes)
                print(x)
                x.extend(inputRowBytes)
            print(x)
            x.append(int(score))
            print(x)
            outputRowBytes = x.tobytes()
            print('output row size: ' + str(outputRowSize))
            print('output row bytes: ' + str(outputRowBytes))
            assert(outputRowSize == len(outputRowBytes))
            tempWriter.write(outputRowBytes)
        tempWriter.flush()
        if reader:
            #none left
            assert(len(reader.read(1)) == 0)
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
        if len(rowBytes) < self.rowSize:
            raise StopIteration
        else:
            rowArray = array.array(self.typecode, rowBytes)
            return rowArray
"""

class AbstractScorer(ABC):
    @abstractmethod
    def scorePeptides(self, alleleName, pepIter):
        pass

class DummyScorer(AbstractScorer):
    def __init__(self, offset):
        self.offset = offset
    def scorePeptides(self, alleleName, pepIter):
        for x in pepIter:
            yield sum([ord(y) - ord('A') for y in x]) + self.offset                             

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

def replaceNoneZero(l):
    return [0 if x == None else x for x in l]
    
class ChainLink:
    def __init__(self, sequenceStart, nextChainIndex, nextChainPosition, lastChainIndex, lastChainPosition):
        self.sequenceStart = sequenceStart
        self.nextChainIndex = nextChainIndex
        self.nextChainPosition = nextChainPosition
        self.lastChainIndex = lastChainIndex
        self.lastChainPosition = lastChainPosition
    def __str__(self):
        return '({}, {}, {}, {}, {})'.format(self.sequenceStart, self.nextChainIndex, self.nextChainPosition, self.lastChainIndex, self.lastChainPosition)
    @property
    def hashAttributeList(self):
        #must return a list of integers
        return replaceNoneZero([self.sequenceStart, self.nextChainIndex, self.nextChainPosition, self.lastChainIndex, self.lastChainPosition])


class Chain(UserList):
    def __init__(self, proteinLength, header):
        self.proteinLength = proteinLength
        self.header = header
        self.data = []

    def __str__(self):
        return 'length: {}, header: {}, links: {}'.format(self.proteinLength, self.header, '->'.join([str(x) for x in self.data]))
    @property
    def hashAttributeList(self):
        #must return a list of integers
        return [self.proteinLength, len(self.data)] + list(itertools.chain.from_iterable([x.hashAttributeList for x in self.data]))

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
    def checksum(self):
        hasher = hashlib.md5()
        attributeList = list(itertools.chain.from_iterable([x.hashAttributeList for x in self.data]))
        a = array.array('L', attributeList)
        hasher.update(a.tobytes())
        #16 byte hash
        return hasher.digest()
    @property
    def peptideLength(self):
        return self.pepLen
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
                    assert(protIndex is not None)
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
                        """
                        If x.lastChainIndex is i, then we need to reference the tempChain, because we haven't replaced self.data[i] with tempChain yet.
                        """
                        if x.lastChainIndex == i:
                            assert(tempChain[x.lastChainPosition].nextChainIndex is not None)
                            tempChain[x.lastChainPosition].nextChainPosition = len(tempChain)
                        else:
                            assert(self.data[x.lastChainIndex][x.lastChainPosition].nextChainIndex is not None)
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
        assert(self.chainIndex is not None)
        assert(self.chainIndex < len(self.chains))
        chainIndex = self.chainIndex
        chainPos = self.chainPosition
        headers = set([self.chains[chainIndex].header])
        while chainPos != None:
            if chainIndex is None:
                print('chain index is none')
                print('chain pos')
                print(chainPos)
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

+------------+----------------+----------------+--------------+--------------+--------------------------------------+---------+
| Chain Hash | Multiplier | Score Category | Peptide Length | Typecode | Num Peptides | Allele Delim | Size of following alleles (in bytes) | Alleles |
+------------+----------------+----------------+--------------+--------------+--------------------------------------+---------+

The chain hash is a 16 byte (128 bit) md5 checksum of the chains used to create the score table. The idea here is that if the FASTA is the same, then the peptides are defined by their length and starting positions within the proteins. 

The multiplier, score category, peptide length, typecode, num Peptides, allele delim, and size of alleles are stored as a packed struct using Python's struct module.

Score category is one of the values for the ScoreCategory enum, indicating what kind of score it is. NetMHCPan gives 5 scores: elution score, elution rank, binding score, binding rank, predicted affinity. NetMHC gives 2 scores: Predicted affinity and rank.

Typecode is the typecode used in the array module for storing the NetMHC scores.

The list of allele names are stored as an ASCII string, with 'allele delim' seperating the different names. 

Note that both allele delim and typecode are stored as ascii characters (so a signed byte). But, pass them to the constructor as normal python strings. 
    """
    STRUCT_FORMAT='!16sHIIcQcL'
    def __init__(self, multiplier, alleleSeperator, alleleList, numPeptides, scoreTypecode, scoreCategory, peptideLength, chainHash):
        self.multiplier = multiplier        
        self.alleleSeperator = alleleSeperator.encode('ascii')
        assert(len(self.alleleSeperator) == 1)
        self.alleleList = [x.encode('ascii') for x in alleleList]
        self.numPeptides = numPeptides
        self._peptideLength = peptideLength
        self.scoreTypecode = scoreTypecode.encode('ascii')
        assert(len(self.scoreTypecode) == 1)
        assert(isinstance(scoreCategory, ScoreCategory))
        self.scoreCategory = scoreCategory
        self._chainHash = chainHash
    @property
    def chainHash(self):
        return self._chainHash
    @property
    def peptideLength(self):
        return self._peptideLength
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
        structSize = struct.calcsize(cls.STRUCT_FORMAT)
        chainHash, multiplier, scoreCategoryInt, peptideLength, scoreTypecode, numPeptides, alleleSeperator, alleleArraySize = struct.unpack(cls.STRUCT_FORMAT, byteArray[0:structSize])        
        alleleArrayBytes = byteArray[structSize:(structSize + alleleArraySize)]
        alleleList = [x.decode('ascii') for x in alleleArrayBytes.split(alleleSeperator)]
        return cls(multiplier, alleleSeperator.decode('ascii'), alleleList, numPeptides, scoreTypecode.decode('ascii'), ScoreCategory(scoreCategoryInt), peptideLength, chainHash)
    def toBytes(self):
        alleleArrayBytes = self.alleleSeperator.join(self.alleleList)
        metaHeader = struct.pack(self.STRUCT_FORMAT, self.chainHash, self.multiplier, self.scoreCategory.value, self.peptideLength, self.scoreTypecode, self.numPeptides, self.alleleSeperator, len(alleleArrayBytes))
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
    @property
    def numPeptides(self):
        return self.tableMeta.getNumPeptides()
    @property
    def chainHash(self):
        return self.tableMeta.chainHash
    @property
    def peptideLength(self):
        return self.tableMeta.peptideLength
    @classmethod
    def empty(cls, fileObj, scoreCategory, scoreTypecode, multiplier, alleleSeperator, peptideLength, chainHash):
        #chain hash should be a bytes object
        meta = ScoreTableMeta(multiplier, alleleSeperator, [], 0, scoreTypecode, scoreCategory, peptideLength, chainHash)
        return cls(fileObj, meta)
    @classmethod
    def readExisting(cls, fileObj):
        size = struct.calcsize('!IL')        
        magic, metaOffsetFromEnd = struct.unpack('!IL', fileObj.read(size))
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
            assert(numScoresWritten == self.tableMeta.getNumPeptides())
            self.tableMeta.addAllele(allele)                
            self.writeMeta()
            self.fileObj.flush()
            return True
    def appendScores(self, scoreIter, location):
        chunkSize = 2**10
        self.fileObj.seek(location)
        numScores = 0
        multiplier = self.tableMeta.getMultiplier()
        while True:
            chunk = [int(x*multiplier) for x in itertools.islice(scoreIter, chunkSize)]
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
        self.fileObj.write(struct.pack('!IL', self.MAGIC, endLoc - beginningOffset))
        
    def scoreIter(self, allele):
        if allele in self.tableMeta.getAlleles():
            loc = struct.calcsize('!IL') + self.tableMeta.computeAlleleOffset(allele)
            self.fileObj.seek(loc)
            typecode = self.tableMeta.getScoreTypecode()            
            size = struct.calcsize(typecode)
            for i in range(0, self.tableMeta.getNumPeptides()):
                score,  = struct.unpack(typecode, self.fileObj.read(size))
                yield float(score)/self.tableMeta.getMultiplier()

class ScoreTableGroup:
    def __init__(self, *scoreTables):
        self.scoreTables = scoreTables
        self.numScoreTables = len(scoreTables)
    def addAllele(self, allele, scoreIter):
        if any([allele in table.tableMeta.getAlleles() for table in self.scoreTables]):
            return False
        else:
            locations = [struct.calcsize('!IL') + table.tableMeta.computeScoreSectionSizeInBytes() for table in self.scoreTables]
            numScoresWritten = self.appendScores(scoreIter, locations)
            for table in self.scoreTables:
                if len(table.tableMeta.alleleList) == 0:
                    table.tableMeta.setNumPeptides(numScoresWritten)
                assert(numScoresWritten == table.tableMeta.getNumPeptides())
                table.tableMeta.addAllele(allele)
                table.writeMeta()
                table.fileObj.flush()
            return True
    def appendScores(self, scoreIter, locations):
        chunkSize = 2**10
        for i in range(0, self.numScoreTables):
            self.scoreTables[i].fileObj.seek(locations[i])
        numScores = 0
        while True:
            chunk = [list(x) for x in itertools.islice(scoreIter, chunkSize)]
            numScores += len(chunk)
            if len(chunk) == 0:
                break
            for i in range(0, self.numScoreTables):
                table = self.scoreTables[i]
                multiplier = table.tableMeta.getMultiplier()
                chunkForTable = [int(x[i]*multiplier) for x in chunk]
                chunkArray = array.array(table.tableMeta.getScoreTypecode(), chunkForTable)
                chunkBytes = chunkArray.tobytes()
                table.fileObj.write(chunkBytes)
        return numScores
        


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

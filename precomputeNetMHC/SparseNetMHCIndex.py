import io
import argparse
import os
import pickle
from bitarray import bitarray
import typing
import sys
import functools
import struct
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable, ScoreCategory, fileMD5, ScoreTableGroup

"""
TODO: 
1) Complete sparsify in SparseScoreTable -- done
2) Complete scoreIter in SparseScoreTable -- done
3) Finish SparseScoresForAllele
4) Update computeAlleleOffset and computeScoreSectionSizeInBytes in SparseScoreTableMeta
"""

class SparseScoreTableMeta:
    """   
    In terms of bytes, this is structured as follows:

+------------+----------------+----------------+--------------+--------------+--------------------------------------+---------+
| Chain Hash | Multiplier | Score Category | Peptide Length | Typecode | Num Peptides | Num sparse peptides | Allele Delim | Size of following alleles (in bytes) | Alleles |
+------------+----------------+----------------+--------------+--------------+--------------------------------------+---------+

The chain hash is a 16 byte (128 bit) md5 checksum of the chains used to create the score table. The idea here is that if the FASTA is the same, then the peptides are defined by their length and starting positions within the proteins. 

The multiplier, score category, peptide length, typecode, num Peptides, num sparse peptides, allele delim, and size of alleles are stored as a packed struct using Python's struct module.

Score category is one of the values for the ScoreCategory enum, indicating what kind of score it is. NetMHCPan gives 5 scores: elution score, elution rank, binding score, binding rank, predicted affinity. NetMHC gives 2 scores: Predicted affinity and rank.

Typecode is the typecode used in the array module for storing the NetMHC scores.

The list of allele names are stored as an ASCII string, with 'allele delim' seperating the different names. 

Note that both allele delim and typecode are stored as ascii characters (so a signed byte). But, pass them to the constructor as normal python strings. 
    """
    STRUCT_FORMAT='!16sHIIcQQcL'
    def __init__(self, multiplier, alleleSeperator, alleleList, numPeptides, numSparsePeptides, scoreTypecode, scoreCategory, peptideLength, chainHash):
        self.multiplier = multiplier        
        self.alleleSeperator = alleleSeperator.encode('ascii')
        self._numSparsePeptides = numSparsePeptides
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
    @property
    def numSparsePeptides(self):
        return self._numSparsePeptides
    @classmethod
    def sparsifyScoreTableMeta(cls, scoreTableMeta, numSparsePeptides):
        return cls(scoreTableMeta.multiplier, scoreTableMeta.alleleSeperator, scoreTableMeta.alleleList, scoreTableMeta.numPeptides, numSparsePeptides, scoreTableMeta.scoreTypecode, scoreTableMeta.scoreCategory, scoreTableMeta.peptideLength, scoreTableMeta.chainHash)
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
        chainHash, multiplier, scoreCategoryInt, peptideLength, scoreTypecode, numPeptides, numSparsePeptides, alleleSeperator, alleleArraySize = struct.unpack(cls.STRUCT_FORMAT, byteArray[0:structSize])        
        alleleArrayBytes = byteArray[structSize:(structSize + alleleArraySize)]
        alleleList = [x.decode('ascii') for x in alleleArrayBytes.split(alleleSeperator)]
        return cls(multiplier, alleleSeperator.decode('ascii'), alleleList, numPeptides, numSparsePeptides, scoreTypecode.decode('ascii'), ScoreCategory(scoreCategoryInt), peptideLength, chainHash)
    def toBytes(self):
        alleleArrayBytes = self.alleleSeperator.join(self.alleleList)
        metaHeader = struct.pack(self.STRUCT_FORMAT, self.chainHash, self.multiplier, self.scoreCategory.value, self.peptideLength, self.scoreTypecode, self.numPeptides, self.numSparsePeptides, self.alleleSeperator, len(alleleArrayBytes))
        return metaHeader + alleleArrayBytes

class ScoreAndIndex:
    def __init__(self, score, index, direction):
        self.score = score
        self.index = index
        self.direction = direction
    def __lt__(self, other):
        if self.direction > 0:
            return self.score >= other.score
        if self.direction < 0:
            return self.score < other.score
    
class BestScoreStore:
    def __init__(self, approxSize, direction):
        #direction is positive if lower scores are better, or negative if higher scores are better
        self.heap = []
        self.approxSize = approxSize
        self.worstValueInHeap = None
        self.worstValueCounter = 0
        self.direction = direction
    def add(self, x, index):
        """
        This is so complicated because we have to deal with ties at the boundary. 
        """
        s = ScoreAndIndex(x, index, self.direction)
        if self.heapMax is None and self.heapMin is None:
            assert(len(self.heap) == 0)
            self.worstValueInHeap = x
            self.worstValueCounter = 1
            heapq.heappush(self.heap, s)
        elif len(self.heap) < self.approxSize:
            if s < ScoreAndIndex(self.worstValueInHeap, 0, self.direction):
                self.worstValueInHeap = x
                self.worstValueCounter = 1
            elif x == self.worstValueInHeap:
                self.worstValueCounter += 1
            heapq.heappush(self.heap, s)
        elif x == self.worstValueInHeap:
            self.worstValueCounter += 1
            heapq.heappush(self.heap, s)
        elif ScoreAndIndex(self.worstValueInHeap, 0, self.direction) < s:
            heapq.heappush(self.heap, s)
            if len(self.heap) - self.worstValueCounter >= self.approxSize:                
                for i in range(0, self.worstValueCounter):
                    worst = heapq.heappop(self.heap)
                    assert(worst.score == self.worstValueInHeap)
                self.worstValueCounter = 0
                self.worstValueInHeap = None
    def getHeap(self):
        return self.heap
            
    

class SparseScoresForAllele:
    @property
    def bitArray(self):
        return self._bitArray
    @property
    def scoreVector(self):
        return self._scoreVector
    @classmethod
    def sparsify(cls, scoreIter, approxNumSparse, direction):
        #numSparse is the number of peptides to store scores for.
        #direction is positive if lower scores are better (like affinity), or negative if higher scores are better.
        store = BestScoreStore(approxNumSparse, direction)
        for score, index in zip(scoreIter, itertools.count(0)):
            store.add(score, index)
        #sort based on index
        sortedScoresAndIndices = sorted(store.getHeap(), key=lambda x: x.index)
        self._scoreVector = [x.index for x in sortedScoresAndIndices]
        #Constuct the bit array
        self._bitArray = bitarray(endian='big')
        lastIndex = 0
        for scoreAndIndex in sortedScoresAndIndices:
            index = scoreAndIndex.index
            self._bitArray.extend([0]*(index - lastIndex))
            self._bitArray.append(1)
            lastIndex = index + 1
    @classmethod
    def fromBytes(cls, fileObj):
        #finish fromBytes
        size = int(fileObj.read(struct.calcsize('Q!')))
        pass
    def toBytes(self):
        bitArrayBytes = self.bitArray.tobytes()
        size = len(bitArrayBytes)
        return struct.pack('Q!', size) + bitArrayBytes
    def scoreIter(self):
        pass
        
class SparseScoreTable:
    """
    The ScoreTable file should start with the MAGIC 4 bytes (4344342b in hex. It's just the ascii string 'CD4+' converted to a hex number). 
    """
    MAGIC=0X4344342b
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
    def numSparsePeptides(self):
        return self.tableMeta.numSparsePeptides
    @property
    def peptideLength(self):
        return self.tableMeta.peptideLength
    @classmethod
    def sparsify(cls, scoreTableObject, sparseFraction, fileObject):
        #fileObject is where to store the SparseScoreTable
        startSize = struct.calcsize('!IL')
        fileObject.write(struct.pack('!IL', self.MAGIC, 0))
        numSparsePeptides = int(sparseFraction*scoreTableObject.tableMeta.numPeptides)
        sparseMeta = SparseScoreTableMeta.sparsifyScoreTableMeta(scoreTableObject.tableMeta, numSparsePeptides)
        for allele in self.tableMeta.getAlleles():
            scoreIter = scoreTableObject.scoreIter(allele)
            sparseScores = SparseScoresForAllele.sparsify(scoreIter, numSparsePeptides)
            sparseScoresBytes = sparseScores.toBytes()
            fileObject.seek(0, os.SEEK_END)
            fileObject.write(sparseScoreBytes)
        
        metaBytes = self.tableMeta.toBytes()
        fileObject.seek(0, os.SEEK_END)
        fileObject.write(metaBytes)
        fileObject.seeky(0, os.SEEK_END)
        endLoc = fileObject.tell()
        fileObject.seek(0)
        fileObject.write(struct.pack('!IL', self.MAGIC, len(metaBytes)))
        
    @classmethod
    def readExisting(cls, fileObj):
        size = struct.calcsize('!IL')        
        magic, metaOffsetFromEnd = struct.unpack('!IL', fileObj.read(size))
        assert(magic == cls.MAGIC)
        fileObj.seek(-1*metaOffsetFromEnd, os.SEEK_END)
        metaBytes = fileObj.read()
        meta = SparseScoreTableMeta.fromBytes(metaBytes)
        return cls(fileObj, meta)
        
    def scoreIter(self, allele):
        if allele in self.tableMeta.getAlleles():
            loc = struct.calcsize('!IL') + self.tableMeta.computeAlleleOffset(allele)
            self.fileObj.seek(loc)
            typecode = self.tableMeta.getScoreTypecode()
            sparseScoresObject = SparseScoresForAllele.fromBytes(self.fileObj, self.tableMeta.getNumPeptides(), self.tableMeta.numSparsePeptides, typecode)
            for x in sparseScoresObject.scoreIter():
                if x is None:
                    yield None
                else:
                    yield float(x)/self.tableMeta.getMultiplier()

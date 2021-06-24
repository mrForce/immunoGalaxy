import io
import argparse
import os
import pickle
import typing
import sys
import functools
from precomputedNetMHCIndex import ChainCollection, peptideGenerator, ScoreTable, ScoreCategory, fileMD5, ScoreTableGroup
from netMHCCalling import NetMHCScorer, NetMHCRunFailedError


"""
TODO: 
1) Complete sparsify in SparseScoreTable
2) Complete scoreIter in SparseScoreTable
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


class SparseScoresForAllele:
    @property
    def bitVector(self):
        return self._bitVector
    @property
    def scoreVector(self):
        return self._scoreVector
    @classmethod
    def sparsify(cls, scoreIter, numSparse):
        pass
    @classmethod
    def fromBytes(cls, byteArray):
        pass
    def toBytes(self):
        pass
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
        return self._numSparsePeptides
    @property
    def peptideLength(self):
        return self.tableMeta.peptideLength
    @classmethod
    def sparsify(cls, scoreTableObject, sparseFraction, fileObject):
        #fileObject is where to store the SparseScoreTable
        return cls(scoreTableObject, sparseFraction)
    
    @classmethod
    def readExisting(cls, fileObj):
        size = struct.calcsize('!IL')        
        magic, metaOffsetFromEnd = struct.unpack('!IL', fileObj.read(size))
        assert(magic == cls.MAGIC)
        fileObj.seek(-1*metaOffsetFromEnd, os.SEEK_END)
        metaBytes = fileObj.read()
        meta = ScoreTableMeta.fromBytes(metaBytes)
        return cls(fileObj, meta)

    def appendScores(self, scoreIter, location):
        chunkSize = 2**10
        self.fileObj.seek(location)
        numScores = 0
        multiplier = self.tableMeta.getMultiplier()
        while True:
            if self.tableMeta.getScoreTypecode() in ['f', 'd']:
                chunk = [float(x*multiplier) for x in itertools.islice(scoreIter, chunkSize)]
            else:
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


parser = argparse.ArgumentParser(description='Sparsify index. In other words, just store the scores for the top X percent of scorers.')
parser.add_argument('inputIndex')
parser.add_argument('sparseOutput')
parser.add_argument('percent', type=float, help='Store scores for top percent percent')
args = parser.parse_args()
percent = args.percent
assert(percent > 0)
assert(percent <= 100.0)
assert(os.path.isfile(args.inputIndex))
assert(not os.path.isfile(args.sparseOutput))
scoreTableFile = open(args.inputIndex, 'rb+')
scoreTable = ScoreTable.readExisting(scoreTableFile)


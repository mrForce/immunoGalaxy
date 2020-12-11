import io
import array
import os
import itertools
import pdb
import csv
import subprocess
import shutil
import threading
from Bio import SeqIO
import tempfile
from collections import Counter, namedtuple, defaultdict, UserList
import queue
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



def getMemoryUsage(self, pids):
    #call smem to determine the USS for each pid
    pass
            
def writePeptidesToFile(peptides, filePath):
    with open(filePath, 'w') as f:
        for x in peptides:
            f.write(x + '\n')
def runNetMHC(peptides, allele, netmhcPath):
    with tempfile.NamedTemporaryFile(mode='w') as inputFile:
        for x in peptides:
            inputFile.write(x + '\n')
            inputFile.flush()
        fd, outputFilePath = tempfile.mkstemp(suffix='.xls')
        command = [netmhcPath, '-a', allele, '-p', '-xls', '-xlsfile', outputFilePath, '-f', inputFile.name]
        print('command: ' + ' '.join(command))
        subprocess.call(command)
        #return a list of tuple (peptide, score)
        with open(outputFilePath, 'r') as outputFile:
            outputFile.readline()
            results = []
            reader = csv.DictReader(outputFile, delimiter='\t')
            for line in reader:
                results.append((line['Peptide'], line['nM']))                
            return results

class NetMHCRunnerThread(threading.Thread):
    def __init__(self, netmhcPath, allele, inputQ, outputQ, completedFlag):
        threading.Thread.__init__(self)
        self.netmhcPath = netmhcPath
        self.allele = allele
        self.inputQ = inputQ
        self.outputQ = outputQ
        self.completedFlag = completedFlag

    def run(self):
        while True:
            batch = None
            pdb.set_trace()
            try:
                batch = self.inputQ.get(timeout=1)
            except queue.Empty:
                if self.completedFlag.is_set():
                    return
            else:                
                assert(batch)
                results = runNetMHC(batch, self.allele, self.netmhcPath)
                self.outputQ.put(results)
                self.inputQ.task_done()

class QueueInserterThread(threading.Thread):
    def __init__(self, q, iterator, batchSize, completedFlag):
        threading.Thread.__init__(self)
        self.q = q
        self.iterator = iterator
        self.batchSize = batchSize
        self.completedFlag = completedFlag
    def run(self):
        while True:
            batch = list(itertools.islice(self.iterator, 0, self.batchSize))
            if batch:
                self.q.put(batch)
            else:
                self.completedFlag.set()
                return
class NetMHCScorer(AbstractScorer):
    def __init__(self, path, batchSize):
        self.path = path
        self.batchSize = batchSize
    def scorePeptides(self, alleleName, pepIter):
        num_threads = 1
        inputQ = queue.Queue(num_threads)
        outputQ = queue.Queue(num_threads)
        completeFlag = threading.Event()
        threads = []
        qThread = QueueInserterThread(inputQ, pepIter, self.batchSize, completeFlag)
        qThread.start()
        for t in range(0, num_threads):
            thread = NetMHCRunnerThread(self.path, alleleName, inputQ, outputQ, completeFlag)
            threads.append(thread)
        for t in threads:
            t.start()
        while True:
            resultBatch = None
            try:
                resultBatch = outputQ.get(timeout=30)
            except queue.Empty:
                if not any([t.is_alive() for t in threads]):
                    break
            else:
                assert(resultBatch)
                for x in resultBatch:
                    yield x
        for t in threads:
            #they should all be gone by now, but do this for safety
            t.join()
            
            

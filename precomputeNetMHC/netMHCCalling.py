import threading
import tempfile
import subprocess
import csv
import queue
import os
import pdb
import itertools
from precomputedNetMHCIndex import AbstractScorer

class NetMHCError(Exception):
    pass
class NetMHCRunFailedError(NetMHCError):
    def __init__(self, command):
        self.message = 'NetMHC Run failed. This is a terminal error; none of the NetMHC results will be saved for this allele-length. Command: ' + command
    

def writePeptidesToFile(peptides, filePath):
    with open(filePath, 'w') as f:
        for x in peptides:
            f.write(x + '\n')
def runNetMHC(peptides, allele, netmhcPath):
    fdIn, inputFilePath  = tempfile.mkstemp()
    numPeptides = 0
    with open(inputFilePath, 'w') as f:
        for x in peptides:
            f.write(x + '\n')
            numPeptides += 1
    fd, outputFilePath = tempfile.mkstemp(suffix='.xls')
    command = [netmhcPath, '-a', allele, '-p', '-xls', '-xlsfile', outputFilePath, '-f', inputFilePath]
    rc = subprocess.call(command)
    assert(rc == 0)
    results = []    
    #return a list of tuple (peptide, score)
    with open(outputFilePath, 'r') as outputFile:
        outputFile.readline()
        reader = csv.DictReader(outputFile, delimiter='\t')
        for line in reader:
            results.append((line['Peptide'], float(line['nM'])))
    assert(len(results) == numPeptides)
    os.remove(inputFilePath)
    os.remove(outputFilePath)
    return (' '.join(command), results)

class CoordinatedOutput:
    def __init__(self, startCount):
        self.counter = startCount
        self.data = None
        self.dataSet = False
        self.cv = threading.Condition()
        self.terminate = False
    def waitToSet(self, countTarget, data):
        with self.cv:
            self.cv.wait_for(lambda: self.counter == countTarget)
            self.data = data
            self.dataSet = True
            self.cv.notify_all()
    def getDataAndIncrement(self):
        with self.cv:
            self.cv.wait_for(lambda: self.dataSet)
            self.dataSet = False
            self.counter += 1
            data = self.data
            self.cv.notify_all()
            return data

class NetMHCRunResults:
    def __init__(self, success, command, data):
        self.success = success
        self.command = command
        self.data = data
class NetMHCRunnerThread(threading.Thread):
    def __init__(self, netmhcPath, allele, inputQ, coordOutput):
        threading.Thread.__init__(self)
        self.netmhcPath = netmhcPath
        self.allele = allele
        self.inputQ = inputQ
        self.coordOutput = coordOutput

    def run(self):
        while True:
            batch = None
            i = -1
            try:
                i, batch = self.inputQ.get(timeout=1)
            except queue.Empty:
                return
            else:
                if batch is None:
                    self.inputQ.task_done()
                    self.coordOutput.waitToSet(i, None)
                    return
                else:
                    try:
                        command, results = runNetMHC(batch, self.allele, self.netmhcPath)
                    except:
                        self.coordOutput.waitToSet(i, NetMHCRunResults(False, ' '.join(command), None))
                        return
                    else:
                        self.inputQ.task_done()
                        resultsDict = dict(results)
                        self.coordOutput.waitToSet(i, NetMHCRunResults(True, ' '.join(command), [resultsDict[x] for x in batch]))

class QueueInserterThread(threading.Thread):
    def __init__(self, q, iterator, batchSize, numThreads):
        threading.Thread.__init__(self)
        self.q = q
        self.iterator = iterator
        self.batchSize = batchSize
        self.numThreads = numThreads
        
    def run(self):
        i = 0
        while True:
            batch = list(itertools.islice(self.iterator, 0, self.batchSize))
            if batch:
                self.q.put((i, batch))
                i += 1
            else:
                for x in range(0, self.numThreads):                    
                    self.q.put((x + i, None))
                return
class NetMHCScorer(AbstractScorer):
    def __init__(self, path, batchSize):
        self.path = path
        self.batchSize = batchSize
    def scorePeptides(self, alleleName, pepIter):
        num_threads = 2
        inputQ = queue.PriorityQueue(num_threads)
        coordOutput = CoordinatedOutput(0)
        threads = []
        qThread = QueueInserterThread(inputQ, pepIter, self.batchSize, num_threads)
        qThread.start()
        for t in range(0, num_threads):
            thread = NetMHCRunnerThread(self.path, alleleName, inputQ, coordOutput)
            threads.append(thread)
        for t in threads:
            t.start()
        deadThreads = 0
        while deadThreads < num_threads:
            resultBatch = coordOutput.getDataAndIncrement()
            if resultBatch is None:
                deadThreads += 1
            else:
                success = resultBatch.success
                if success:
                    assert(resultBatch)
                    for x in resultBatch.data:
                        yield x
                else:
                    raise NetMHCRunFailedError(resultBatch.command)
                        
        for t in threads:
            #they should all be gone by now, but do this for safety
            t.join()

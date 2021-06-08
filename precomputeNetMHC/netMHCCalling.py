import threading
import tempfile
import subprocess
import csv
import queue
import os
import pdb
import itertools
from precomputedNetMHCIndex import AbstractScorer
#retry 2 times if necessary, for a total of 3 runs for a single batch of peptides.
NUM_RETRY=2
class NetMHCError(Exception):
    pass
class NetMHCRunFailedError(NetMHCError):
    def __init__(self, runs):
        self.runs = runs
        messageLines = []
        messageLines.append('NetMHC Run failed for all retries')
        messageLines.append('Input File (containing peptides): ' + runs.peptidesFilePath)
        messageLines.append('This was ran ' + str(len(runs.runs)) + ' times')
        for x in runs.runs:
            messageLines.append('command: ' + x.command)
            messageLines.append('return code: ' + str(x.returnCode))
            messageLines.append('stdout file: ' + x.stdoutPath)
            messageLines.append('stderr file: ' + x.stderrPath)
        self.message = '\n'.join(messageLines)




class NetMHCRun:
    def __init__(self, stdoutPath, stderrPath, returnCode, command):
        self.stdoutPath = stdoutPath
        self.stderrPath = stderrPath
        self.returnCode = returnCode
        self.command = command
        self.success = True
        self.results = None
    def setResults(self, results):
        self.results = results
    def getResults(self):
        return self.results
    def setFailed(self):
        self.success = False
    def __del__(self):
        if self.success:
            print('removing files')
            print(self.stdoutPath)
            print(self.stderrPath)
            os.remove(self.stdoutPath)
            os.remove(self.stderrPath)
        
class NetMHCRuns:
    def __init__(self, peptidesFilePath):
        self.runs = []
        self.peptidesFilePath = peptidesFilePath
        self.success = True
    def __del__(self):
        if self.success:
            os.remove(self.peptidesFilePath)
    def setFailed(self):
        self.success = False
    def addRun(self, netmhcRun):
        self.runs.append(netmhcRun)
    def getRuns(self):
        return self.runs
    def getLastRun(self):
        return self.runs[-1]
    def getPeptidesFilePath(self):
        return self.peptidesFilePath

def extractScores(outputFilePath, fieldsToExtract):
    results = []
    with open(outputFilePath, 'r') as outputFile:
        lines = outputFile.readlines()
        fields = []
        i = 0
        while i < len(lines):
            line = lines[i]
            if fields:
                if line.strip().startswith('-'):
                    fields = []
                else:
                    parts = line.split()
                    lineDict = dict(zip(fields, parts))
                    if 'peptide' in lineDict and 'Peptide' not in lineDict:
                        lineDict['Peptide'] = lineDict['peptide']
                        
                    results.append((lineDict['Peptide'], [float(lineDict[x]) for x in fieldsToExtract]))
                i += 1
            elif line.strip().lower().startswith('pos') and i >= 1 and lines[i - 1].startswith('-') and i + 1 < len(lines) and lines[i + 1].startswith('-'):
                fields = [x.strip() for x in line.split()]
                i += 2
            else:
                i += 1
    return results
def writePeptidesToFile(peptides, filePath):
    with open(filePath, 'w') as f:
        for x in peptides:
            f.write(x + '\n')

KNOWN_AMINOS=set('ACDEFGHIKLMNPQRSTVWY')
def unknownSub(peptides):
    pep = []
    for x in peptides:
        sub = ''.join([y if y in KNOWN_AMINOS else 'X' for y in x])
        pep.append(sub)
    return pep
    
def runNetMHC(peptides, commandGenerator, fieldsToExtract):
    #commandGenerator is a function that takes in the input file path and output file path, and returns the command to run
    fdIn, inputFilePath  = tempfile.mkstemp()
    with open(inputFilePath, 'w') as f:
        for x in peptides:
            f.write(x + '\n')
    os.close(fdIn)
    runs = NetMHCRuns(inputFilePath)
    tries = 0
    while tries <= NUM_RETRY:
        stdoutFD, stdoutPath = tempfile.mkstemp()
        command = commandGenerator(inputFilePath)
        stderrFD, stderrPath = tempfile.mkstemp()
        print('stdout: ' + stdoutPath + ', stderr: ' + stderrPath)
        rc = subprocess.call(command, stdout=stdoutFD, stderr=stderrFD)
        os.close(stdoutFD)
        os.close(stderrFD)
        run = NetMHCRun(stdoutPath, stderrPath, rc, ' '.join(command))
        runs.addRun(run)
        tries += 1
        if rc == 0:
            scores = extractScores(stdoutPath, fieldsToExtract)
            scoredPeptides = set([x[0] for x in scores])
            pepSub = unknownSub(peptides)
            peptideSet = set(pepSub)
            if scoredPeptides <= peptideSet and peptideSet <= scoredPeptides:
                scoreDict = dict(scores)
                run.setResults([scoreDict[x] for x in pepSub])
                return runs
        run.setFailed()
    runs.setFailed()
    return runs

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
            status = self.cv.wait_for(lambda: self.dataSet)

            self.dataSet = False
            self.counter += 1
            data = self.data
            self.cv.notify_all()
            return data

class NetMHCRunnerThread(threading.Thread):
    def __init__(self, commandGenerator,  fieldsToExtract, inputQ, coordOutput):
        threading.Thread.__init__(self)
        self.commandGenerator = commandGenerator
        self.fieldsToExtract = fieldsToExtract
        self.inputQ = inputQ
        self.coordOutput = coordOutput

    def run(self):
        while True:
            batch = None
            i = -1
            try:
                i, batch = self.inputQ.get()
            except queue.Empty:                      
                return
            else:
                if batch is None:
                    self.inputQ.task_done()
                    self.coordOutput.waitToSet(i, None)
                    return
                else:
                    runs = runNetMHC(batch, self.commandGenerator, self.fieldsToExtract)
                    success = runs.success
                    self.coordOutput.waitToSet(i, runs)
                    if not success:
                        return

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
    def __init__(self, batchSize, commandGenerator, numThreads):
        self.commandGenerator = commandGenerator
        self.batchSize = batchSize
        self.numThreads = numThreads
    def scorePeptides(self, pepIter, fieldsToExtract):
        inputQ = queue.PriorityQueue(self.numThreads)
        coordOutput = CoordinatedOutput(0)
        threads = []
        qThread = QueueInserterThread(inputQ, pepIter, self.batchSize, self.numThreads)
        qThread.start()
        for t in range(0, self.numThreads):
            thread = NetMHCRunnerThread(self.commandGenerator, fieldsToExtract, inputQ, coordOutput)
            threads.append(thread)
        for t in threads:
            t.start()
        deadThreads = 0
        while deadThreads < self.numThreads:
            runs = coordOutput.getDataAndIncrement()
            if runs is None:
                deadThreads += 1
            else:
                success = runs.success
                if success:
                    for x in runs.getLastRun().getResults():
                        yield x
                else:
                    raise NetMHCRunFailedError(runs)
                        
        #for t in threads:
        #they should all be gone by now, but do this for safety
        #    t.join()

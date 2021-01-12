#!/usr/bin/python
import sys
import argparse
import pickle
import re
import functools
import subprocess
from PINFile import PINFile
from precomputedNetMHCIndex import ChainCollection,  ScoreTable
from netMHCCalling import NetMHCScorer
from filterNetMHC import filterNetMHC
import glob
import tempfile
from collections import defaultdict
import shutil
import os
import itertools
import uuid

#locations of needed tools
NETMHC = '/home/jordan/netMHC-4.0/netMHC'
NETMHCPAN = ''
MSGFPLUS = '/home/jordan/MSGFPlus/MSGFPlus.jar'
CRUX='/home/jordan/crux'
MSGF2PIN='/home/jordan/msgf2pin'
THREADS=16
PRECOMPUTE_SCRIPTS='/home/jordan/github/immunoGalaxy/precomputeNetMHC'

def generateNetMHCCommand(netmhcPath, allele, inputFilePath):
    return [netmhcPath, '-a', allele, '-p', '-f', inputFilePath]



def writePeptideHeaderMapToFasta(peptideToHeaders, fastaPath):
    #pass in a dictionary mapping each peptide to a set of headers
    with open(fastaPath, 'w') as f:
        for k,v in peptideToHeaders.items():
            header = '>' + ' @@ '.join(list(v))
            f.write(header + '\n')
            f.write(k + '\n')
def addRevcat(fastaPath):
    root, ext = os.path.splitext(fastaPath)
    return root + '.revCat' + ext
def getMemory():
    p = subprocess.run(['free'], stdout=subprocess.PIPE)
    s = p.stdout.decode()
    lines = s.split('\n')
    for line in lines[1::]:
        split_line = line.split()
        if split_line[0] == 'Mem:':
            return int(split_line[-1])
    return -1
    
parser = argparse.ArgumentParser()

TEST=False

parser.add_argument('baseDirectory')
parser.add_argument('baseFasta')
parser.add_argument('netmhcScoreDir')
parser.add_argument('netmhcPanScoreDir')
parser.add_argument('--allele', action='append')
parser.add_argument('--additional_proteome', action='append')
parser.add_argument('--mod', action='append')
parser.add_argument('--mode', choices=['filtered', 'unfiltered', 'netMHCPercolator'])
parser.add_argument('--pep_len', type=str)
parser.add_argument('--rank_filter', type=float)
parser.add_argument('--frag_method', type=str)
parser.add_argument('--instrument', type=str)
parser.add_argument('--precursor_tolerance', type=str)
parser.add_argument('--mgf', type=str)
parser.add_argument('--archive', type=str)
parser.add_argument('--msgf_unfiltered', type=str)
parser.add_argument('--percolator_unfiltered', type=str)
parser.add_argument('--num_matches_per_spectrum', type=int)
parser.add_argument('--minLength', type=int, default=0)
parser.add_argument('--maxLength', type=int, default=0)


args = parser.parse_args()
print('args')
print(args)



memory = getMemory()
print('Memory: ' + str(memory))
assert(memory > 0)



allele_list = []
filtered = False
if args.mode == 'filtered':
    filtered = True
peptide_lengths = []

tempDir = tempfile.mkdtemp(prefix='pipelineRun')

outputPath = os.path.join(tempDir, 'search.mzid')
revCatFastaPath = None
pinOutputPath = os.path.join(tempDir, 'search.mzid.pin')
mgf = os.path.join(tempDir, 'spectra.mgf')
shutil.copyfile(args.mgf, mgf)
msgfCommand = ['java', '-Xmx' + str(int(memory/2048)) + 'M', '-jar', MSGFPLUS,
               '-s', mgf,
               '-ignoreMetCleavage', '1',
               '-t', args.precursor_tolerance,
               '-tda', '1',
               '-addFeatures', '1',
               '-n', str(args.num_matches_per_spectrum),
               '-m', args.frag_method,
               '-inst', args.instrument,
               '-o', outputPath
               ]

modPath = os.path.join(tempDir, 'mods.txt')
with open(modPath, 'w') as f:
    if args.mod:
        for mod in args.mod:
            f.write(mod + '\n')
msgfCommand.extend(['-mod', modPath])

if filtered:
    lengths = [int(x) for  x in args.pep_len.split(',')]
    msgfCommand.extend(['-e', '9'])
    msgfCommand.extend(['-minLength', str(min(lengths))])
    msgfCommand.extend(['-maxLength', str(max(lengths))])
    if args.additional_proteome:
        for x in lengths:
            chains = ChainCollection(args.additional_proteome, x)
            chainPath = os.path.join(tempDir, str(x) + 'additional.chains')
            with open(chainPath, 'wb') as f:
                pickle.dump(chains, f)
            scoreTablePath = os.path.join(tempDir, str(x) + 'additional.scores')
            for allele in args.allele:
                precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precompute.py'), NETMHC, args.additional_proteome, chainPath, scoreTablePath, allele, str(x), str(THREADS)]
                proc = subprocess.Popen(precomputeCommand, stdout=subprocess.DEVNULL)
                outs, errors = proc.communicate()
    pepToHeaders = defaultdict(set)
    for x in lengths:
        additionalChains = None
        additionalScoreTable = None
        if args.additional_proteome:
            with open(os.path.join(tempDir, str(x) + 'additional.chains'), 'rb') as f:
                additionalChains = pickle.load(f)
            with open(os.path.join(tempDir, str(x) + 'additional.scores'), 'rb') as f:
                additionalScoreTable = ScoreTable.readExisting(f)
        baseChains = None
        with open(os.path.join(args.baseDirectory, str(x) + '.chains'), 'rb') as f:
            baseChains = pickle.load(f)
        baseScoreTable = None
        scoresFile = open(os.path.join(args.baseDirectory, args.netmhcScoreDir, str(x) + '.scores'), 'rb')
        
        baseScoreTable = ScoreTable.readExisting(scoresFile)
        additionalFasta = args.additional_proteome if args.additional_proteome else None
        for allele in args.allele:
            pepToHeader = filterNetMHC(allele, x, baseScoreTable, baseChains, os.path.join(args.baseDirectory, args.baseFasta), additionalScoreTable, additionalChains, additionalFasta, args.rank_filter/100.0)
            for k,v in pepToHeader.items():
                pepToHeaders[k].update(v)
        scoresFile.close()
    fasta = os.path.join(tempDir, 'peptides.fasta')
    revCatFastaPath  = addRevcat(fasta)
    writePeptideHeaderMapToFasta(pepToHeaders, fasta)
    msgfCommand.extend(['-d', fasta])
    
else:
    msgfCommand.extend(['-e', '0'])
    if args.minLength > 0:
        msgfCommand.extend(['-minLength', str(args.minLength)])
    if args.maxLength > 0:
        msgfCommand.extend(['-maxLength', str(args.maxLength)])
    fasta = os.path.join(args.baseDirectory, args.baseFasta)
    if args.additional_proteome:
        command = ['awk', '1', os.path.join(args.baseDirectory, args.baseFasta), args.additional_proteome]
        fasta = os.path.join(tempDir, 'combined.fasta')
        f = open(fasta, 'w')
        proc = subprocess.Popen(command, stdout=f)
        assert(proc.wait() == 0)
    else:
        fasta = os.path.join(tempDir, 'base.fasta')
        shutil.copyfile(os.path.join(args.baseDirectory, args.baseFasta), fasta)
        
    revCatFastaPath  = addRevcat(fasta)
    msgfCommand.extend(['-d', fasta])
print('going to call RunMSGFPlusSearch. Command: %s' % ' '.join(msgfCommand))
if TEST:
    print('skipping RunMSGFPlusSearch because of TEST')
else:
    p = subprocess.Popen(msgfCommand, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)
print('ran msgfplus search')

"""
Run msgf2pin
"""

msgf2pinCommand = [MSGF2PIN, outputPath, '-o', pinOutputPath, '-e', 'no_enzyme', '-P', 'XXX_', '-m', str(args.num_matches_per_spectrum), '-F', revCatFastaPath]
p = subprocess.Popen(msgf2pinCommand, stderr=sys.stdout.fileno())
assert(p.wait() == 0)

#msgf unfiltered
if args.msgf_unfiltered:
    with open(pinOutputPath, 'r') as f:
        with open(args.msgf_unfiltered, 'w') as g:
            lines = f.readlines()
            #remove the direction row
            lines.pop(1)
            for line in lines:
                g.write(line)

if args.mode == 'netMHCPercolator':
    print('NetMHC + Percolator')
    print('Pin output path: ' + pinOutputPath)
    pin = PINFile(pinOutputPath)
    peptides = list(pin.peptides)
    print('number of peptides: ' + str(len(peptides)))
    scoreDict = defaultdict(list)
    for allele in args.allele:
        commandGen = functools.partial(generateNetMHCCommand, NETMHC, allele)
        scorer = NetMHCScorer(5000, commandGen, 1)
        scores = [x[0] for x in scorer.scorePeptides(peptides, ['Affinity(nM)'])]
        assert(len(scores) == len(peptides))
        for i in range(0, len(scores)):
            scoreDict[peptides[i]].append(scores[i])
    singleScoreDict = {k: min(v) for k,v in scoreDict.items()}
    pin.addScores(singleScoreDict, 'NetMHC', '-1')
    
    
percolatorOutputDir = os.path.join(tempDir, 'percolatorOutput')
percolatorCommand = [CRUX, 'percolator', '--top-match', str(args.num_matches_per_spectrum), '--output-dir', percolatorOutputDir, pinOutputPath]
p = subprocess.Popen(percolatorCommand, stderr=sys.stdout.fileno())
assert(p.wait() == 0)

print('ran percolator')
percolator_target_path = os.path.join(percolatorOutputDir, 'percolator.target.psms.txt')
percolator_decoy_path = os.path.join(percolatorOutputDir, 'percolator.decoy.psms.txt')
assert(os.path.isfile(percolator_target_path))
assert(os.path.isfile(percolator_decoy_path))


with tempfile.NamedTemporaryFile() as f:
    command = ['awk', 'BEGIN {OFS="\t"} NR==1 {print "Label", $0} NR>1 {print 1, $0}', percolator_target_path]
    p = subprocess.Popen(command, stdout=f)
    assert(p.wait() == 0)
    with tempfile.NamedTemporaryFile() as g:
        command = ['awk', 'BEGIN {OFS="\t"} NR>1 {print -1, $0}', percolator_decoy_path]
        p = subprocess.Popen(command, stdout=g)
        assert(p.wait() == 0)
        command = ['cat', f.name, g.name]
        with open(args.percolator_unfiltered, 'w') as h:
            p = subprocess.Popen(command, stdout=h)
            assert(p.wait() == 0)

if args.archive:
    subprocess.run(['zip', '-r', args.archive + '.zip', tempDir])
    shutil.move(args.archive + '.zip', args.archive)
    print('Zip file size: %d' % os.path.getsize(args.archive))

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
NETMHCPAN = '/data5/jordan/netMHCpan-4.1/netMHCpan'
MSGFPLUS = '/home/jordan/MSGFPlus/MSGFPlus.jar'
CRUX='/home/jordan/crux'
MSGF2PIN='/home/jordan/msgf2pin'
THREADS=16
PRECOMPUTE_SCRIPTS='/home/jordan/github/immunoGalaxy/precomputeNetMHC'
MAX_HEADER_LENGTH=5000



def generateNetMHCCommand(netmhcPath, allele, inputFilePath):
    return [netmhcPath, '-a', allele, '-p', '-f', inputFilePath]



def writePeptideHeaderMapToFasta(peptideToHeaders, fastaPath):
    #pass in a dictionary mapping each peptide to a set of headers
    with open(fastaPath, 'w') as f:
        for k,v in peptideToHeaders.items():
            header = '>' + ' @@ '.join(list(v))
            if len(header) > MAX_HEADER_LENGTH:
                header = header[0:MAX_HEADER_LENGTH]
            f.write(header + '\n')
            f.write(k + '\n')
def addRevcat(fastaPath):
    root, ext = os.path.splitext(fastaPath)
    return root + '.revCat' + ext

def create_decoy_file(targetFile, outputPath):
    with tempfile.TemporaryDirectory() as tempDir:
        shutil.copy(targetFile, tempDir)
        buildSACommand = ['java', '-Xmx10000M', '-cp', MSGFPLUS, 'edu.ucsd.msjava.msdbsearch.BuildSA', '-d', os.path.join(tempDir, targetFile), '-tda', '1']
        print('SA command: ' + ' '.join(buildSACommand), flush=True)
        p = subprocess.Popen(buildSACommand, stderr=sys.stdout.fileno())
        assert(p.wait() == 0)
        revcatFile = os.path.join(tempDir, addRevcat(targetFile))
        
        assert(os.path.isfile(revcatFile))
        f = open(outputPath, 'w')
        g = open(revcatFile, 'r')
        copyLines = False
        for line in g:
            if (not copyLines) and line.startswith('>XXX_'):
                copyLines = True
            if copyLines:
                f.write(line)
        f.close()
        g.close()
    
parser = argparse.ArgumentParser()

TEST=False

parser.add_argument('baseDirectory')
parser.add_argument('baseFasta')
parser.add_argument('netmhcScoreDir')
parser.add_argument('netmhcPanScoreDir')
parser.add_argument('mhcFlurryScoreDir')
parser.add_argument('--allele', action='append', nargs=2)
parser.add_argument('--panAllele', action='append', nargs=3)
parser.add_argument('--additional_proteome', type=str)
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
parser.add_argument('--mix_max', action='store_true')
parser.add_argument('--minLength', type=int, default=0)
parser.add_argument('--maxLength', type=int, default=0)

targetFasta = None
decoyFasta = None
args = parser.parse_args()
print('args')
print(args)
if args.allele:
    for x in args.allele:
        assert(x[0] in ['netmhcPrecompute', 'netmhcOnFly', 'MHCFlurryPrecompute', 'MHCFlurryOnFly'])
else:
    args.allele = []
if args.panAllele:
    for x in args.panAllele:
        assert(x[0] in ['netmhcPanPrecompute', 'netmhcPanOnFly'])
        assert(x[2] in ['ba', 'elute'])
else:
    args.panAllele = []
usingBase = True
if args.baseDirectory == 'None':
    assert(args.additional_proteome)
    usingBase = False


allele_list = []
filtered = False
if args.mode == 'filtered':
    filtered = True
peptide_lengths = []

tempDir = tempfile.mkdtemp(prefix='pipelineRun')

target_outputPath = os.path.join(tempDir, 'search.mzid')
decoy_outputPath = os.path.join(tempDir, 'search_decoy.mzid')

target_pinOutputPath = os.path.join(tempDir, 'search.mzid.pin')
decoy_pinOutputPath = os.path.join(tempDir, 'search_decoy.mzid.pin')
combined_pinOutputPath = os.path.join(tempDir, 'search_combined.mzid.pin')
mgf = os.path.join(tempDir, 'spectra.mgf')
shutil.copyfile(args.mgf, mgf)
mixMax = True#args.mix_max


    
msgfCommand = ['java', '-Xmx10000M', '-jar', MSGFPLUS,
               '-s', mgf,
               '-ignoreMetCleavage', '1',
               '-t', args.precursor_tolerance,
               '-addFeatures', '1',
               '-n', str(args.num_matches_per_spectrum),
               '-m', args.frag_method,
               '-tda', '0',
               '-inst', args.instrument]
target_msgfCommand = []
decoy_msgfCommand = []
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
    pepToHeaders = defaultdict(set)
    for x in lengths:
        additionalChains = None
        additionalChainPath = None
        additionalPanBAScoreTablePath = None
        additionalPanEluteScoreTablePath = None
        additionalNetMHCscoreTablePath = None
        additionalFlurryScoreTablePath = None
        if args.additional_proteome:
            additionalChains = ChainCollection(args.additional_proteome, x)
            additionalChainPath = os.path.join(tempDir, str(x) + 'additional.chains')
            with open(additionalChainPath, 'wb') as f:
                pickle.dump(additionalChains, f)
            additionalPanBAScoreTablePath = os.path.join(tempDir, str(x) + 'additional_pan_ba.scores')
            additionalPanEluteScoreTablePath = os.path.join(tempDir, str(x) + 'additional_pan_elute.scores')
            additionalNetMHCscoreTablePath = os.path.join(tempDir, str(x) + 'additional.scores')
            additionalFlurryScoreTablePath = os.path.join(tempDir, str(x) + 'additional_flurry.scores')
        for allele in args.panAllele:
            additionalScoreTable = None
            if args.additional_proteome:
                precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precomputePan.py'), NETMHCPAN, args.additional_proteome, additionalChainPath, additionalPanEluteScoreTablePath, additionalPanBAScoreTablePath, allele[1].strip(), str(x), str(THREADS)]
                proc = subprocess.Popen(precomputeCommand, stdout=subprocess.DEVNULL)
                outs, errors = proc.communicate()
                additionalScoreFile = None
                if allele[2] == 'elute':
                    additionalScoreFile = open(additionalPanEluteScoreTablePath, 'rb')
                elif allele[2] == 'ba':
                    additionalScoreFile = open(additionalPanBAScoreTablePath, 'rb')
                else:
                    assert(0)
                additionalScoreTable = ScoreTable.readExisting(additionalScoreFile)
            baseScoreTable = None
            baseChains = None
            baseFasta = None
            if usingBase:
                baseChainPath = os.path.join(args.baseDirectory, str(x) + '.chains')
                with open(baseChainPath, 'rb') as f:
                    baseChains = pickle.load(f)
                baseScoreFile = None
                baseFasta = os.path.join(args.baseDirectory, args.baseFasta)
                if allele[0] == 'netmhcPanPrecompute':
                    if allele[2] == 'elute':
                        baseScoreFile = open(os.path.join(args.baseDirectory, args.netmhcPanScoreDir, str(x) + '_elute.scores'), 'rb')
                    elif allele[2] == 'ba':
                        baseScoreFile = open(os.path.join(args.baseDirectory, args.netmhcPanScoreDir, str(x) + '_ba.scores'), 'rb')
                    else:
                        assert(0)
                elif allele[0] == 'netmhcPanOnFly':
                    eluteScoresFileName = os.path.join(tempDir, allele[1].strip() + '_' + str(x) + '_elute_netmhcpan_fly.scores')
                    baScoresFileName = os.path.join(tempDir, allele[1].strip() + '_' + str(x) + '_ba_netmhcpan_fly.scores')
                    precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precomputePan.py'), NETMHCPAN, baseFasta, baseChainPath, eluteScoresFileName, baScoresFileName, allele[1].strip(), str(x), str(THREADS)]
                    proc = subprocess.Popen(precomputeCommand, stdout=subprocess.DEVNULL)
                    outs, errors = proc.communicate()
                    if allele[2] == 'elute':
                        baseScoreFile = open(eluteScoresFileName, 'rb')
                    elif allele[2] == 'ba':
                        baseScoreFile = open(baScoresFileName, 'rb')
                    else:
                        assert(0)
                baseScoreTable = ScoreTable.readExisting(baseScoreFile)
            reverse = True if allele[2] == 'elute' else False
            pepToHeader = filterNetMHC(allele[1].strip(), x, baseScoreTable, baseChains, baseFasta, additionalScoreTable, additionalChains, args.additional_proteome,  args.rank_filter/100.0, reverse)
            for k,v in pepToHeader.items():
                pepToHeaders[k].update(v)
        for allele in args.allele:
            additionalScoreTable = None
            if args.additional_proteome:
                additionalScoreTablePath = None
                precomputeCommand = None
                if allele[0] in ['netmhcPrecompute', 'netmhcOnFly']:                    
                    precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precompute.py'), NETMHC, args.additional_proteome, additionalChainPath, additionalNetMHCscoreTablePath, allele[1].strip(), str(x), str(THREADS)]
                    additionalScoreTablePath = additionalNetMHCscoreTablePath
                    print('additional score table path: ' + additionalScoreTablePath)
                elif allele[0] in ['MHCFlurryPrecompute', 'MHCFlurryOnFly']:
                    precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precomputeRemoteFlurry.py'), args.additional_proteome, additionalChainPath, additionalFlurryScoreTablePath, allele[1].strip(), str(x)]
                    additionalScoreTablePath = additionalFlurryScoreTablePath                    
                proc = subprocess.Popen(precomputeCommand, stdout=subprocess.DEVNULL)
                outs, errors = proc.communicate()
                additionalScoreFile = open(additionalScoreTablePath, 'rb')
                additionalScoreTable = ScoreTable.readExisting(additionalScoreFile)
            baseScoreTable = None
            baseChains = None
            baseFasta = None
            if usingBase:
                baseChainPath = os.path.join(args.baseDirectory, str(x) + '.chains')
                with open(baseChainPath, 'rb') as f:
                    baseChains = pickle.load(f)
                baseScoreFile = None
                baseFasta = os.path.join(args.baseDirectory, args.baseFasta)
                if allele[0] == 'netmhcPrecompute':
                    filePath = os.path.join(args.baseDirectory, args.netmhcScoreDir, str(x) + '.scores')
                    baseScoreFile = open(filePath, 'rb')
                elif allele[0] == 'MHCFlurryPrecompute':
                    baseScoreFile = open(os.path.join(args.baseDirectory, args.mhcFlurryScoreDir, str(x) + '_top_alleles.scores'), 'rb')
                elif allele[0] == 'netmhcOnFly':
                    baseScoreFileName = os.path.join(tempDir, allele[1].strip() + '_' + str(x) + '_netmhc_fly.scores')
                    precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precompute.py'), NETMHC, baseFasta, baseChainPath, baseScoreFileName, allele[1].strip(), str(x), str(THREADS)]
                    proc = subprocess.Popen(precomputeCommand, stdout=subprocess.DEVNULL)
                    outs, errors = proc.communicate()
                    baseScoreFile = open(baseScoreFileName, 'rb')
                elif allele[0] == 'MHCFlurryOnFly':
                    baseScoreFileName = os.path.join(tempDir, allele[1].strip() + '_' + str(x) + '_mhcflurry_fly.scores')
                    precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precomputeRemoteFlurry.py'), baseFasta, baseChainPath, baseScoreFileName, allele[1].strip(), str(x)]
                    proc = subprocess.Popen(precomputeCommand, stdout=subprocess.DEVNULL)
                    outs, errors = proc.communicate()
                    baseScoreFile = open(baseScoreFileName, 'rb')
                else:
                    assert(0)
                baseScoreTable = ScoreTable.readExisting(baseScoreFile)
            reverse = True if allele[0] in ['MHCFlurryPrecompute', 'MHCFlurryOnFly'] else False
            pepToHeader = filterNetMHC(allele[1].strip(), x, baseScoreTable, baseChains, baseFasta, additionalScoreTable, additionalChains, args.additional_proteome,  args.rank_filter/100.0, reverse)
            for k,v in pepToHeader.items():
                pepToHeaders[k].update(v)
            

                
    targetFasta = os.path.join(tempDir, 'peptides.fasta')
    decoyFasta = os.path.join(tempDir, 'peptides_decoys.fasta')
    writePeptideHeaderMapToFasta(pepToHeaders, targetFasta)
    create_decoy_file(targetFasta, decoyFasta)
    target_msgfCommand.extend(msgfCommand)
    decoy_msgfCommand.extend(msgfCommand)
    target_msgfCommand.extend(['-d', targetFasta, '-o', target_outputPath])
    decoy_msgfCommand.extend(['-d', decoyFasta, '-o', decoy_outputPath ])
    
else:
    msgfCommand.extend(['-e', '0'])
    if args.minLength > 0:
        msgfCommand.extend(['-minLength', str(args.minLength)])
    if args.maxLength > 0:
        msgfCommand.extend(['-maxLength', str(args.maxLength)])
    if usingBase:
        targetFasta = os.path.join(args.baseDirectory, args.baseFasta)
        if args.additional_proteome:
            command = ['awk', '1', os.path.join(args.baseDirectory, args.baseFasta), args.additional_proteome]
            targetFasta = os.path.join(tempDir, 'combined.fasta')
            f = open(targetFasta, 'w')
            proc = subprocess.Popen(command, stdout=f)
            assert(proc.wait() == 0)
        else:
            targetFasta = os.path.join(tempDir, 'base.fasta')
            shutil.copyfile(os.path.join(args.baseDirectory, args.baseFasta), targetFasta)
    else:
        targetFasta = os.path.join(tempDir, 'additional.fasta')
        shutil.copyfile(args.additional_proteome, targetFasta)
    assert(targetFasta)
    assert(os.path.isfile(targetFasta))
    decoyFasta = os.path.join(tempDir, 'decoy.fasta')
    create_decoy_file(targetFasta, decoyFasta)
    target_msgfCommand.extend(msgfCommand)
    decoy_msgfCommand.extend(msgfCommand)    
    target_msgfCommand.extend(['-d', targetFasta, '-o', target_outputPath])
    decoy_msgfCommand.extend(['-d', decoyFasta, '-o', decoy_outputPath])
    

if TEST:
    print('skipping RunMSGFPlusSearch because of TEST')
else:
    print('going to call RunMSGFPlusSearch. Command: %s' % ' '.join(target_msgfCommand), flush=True)
    p = subprocess.Popen(target_msgfCommand, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)
    print('going to call RunMSGFPlusSearch. Command: %s' % ' '.join(decoy_msgfCommand), flush=True)
    p = subprocess.Popen(decoy_msgfCommand, stderr=sys.stdout.fileno())
    assert(p.wait() == 0)
print('ran msgfplus search')

"""
Run msgf2pin
"""

target_msgf2pinCommand = [MSGF2PIN, target_outputPath, '-o', target_pinOutputPath, '-e', 'no_enzyme', '-m', str(args.num_matches_per_spectrum), '-F', targetFasta]
print('Target msgf2pin command: ')
print(target_msgf2pinCommand)
p = subprocess.Popen(target_msgf2pinCommand, stderr=sys.stdout.fileno())
assert(p.wait() == 0)

decoy_msgf2pinCommand = [MSGF2PIN, decoy_outputPath, '-o', decoy_pinOutputPath, '-e', 'no_enzyme', '-m', str(args.num_matches_per_spectrum), '-F', decoyFasta]
p = subprocess.Popen(decoy_msgf2pinCommand, stderr=sys.stdout.fileno())
assert(p.wait() == 0)

shutil.copyfile(target_pinOutputPath, combined_pinOutputPath)
combinedPIN = PINFile(combined_pinOutputPath)
combinedPIN.addPin(decoy_pinOutputPath, decoy=True)

#msgf unfiltered
if args.msgf_unfiltered:
    with open(combined_pinOutputPath, 'r') as f:
        with open(args.msgf_unfiltered, 'w') as g:
            lines = f.readlines()
            #remove the direction row
            lines.pop(1)
            for line in lines:
                g.write(line)

if args.mode == 'netMHCPercolator':
    print('NetMHC + Percolator')
    print('Combined pin output path: ' + combined_pinOutputPath)
    combinedPin = PINFile(combined_pinOutputPath)
    peptides = list(combinedPin.peptides)
    print('number of  peptides: ' + str(len(peptides)))
    scoreDict = defaultdict(list)
    for allele in args.allele:
        commandGen = functools.partial(generateNetMHCCommand, NETMHC, allele.strip())
        scorer = NetMHCScorer(5000, commandGen, 1)
        scores = [x[0] for x in scorer.scorePeptides(iter(peptides), ['Affinity(nM)'])]
        print('number of scores: ' + str(len(scores)))
        print('scores: ')
        print(scores)
        assert(len(scores) == len(peptides))
        for i in range(0, len(scores)):
            scoreDict[peptides[i]].append(scores[i])
    print('score dict')
    print(scoreDict)
    singleScoreDict = {k: min(v) for k,v in scoreDict.items()}
    print('single score dict')
    print(singleScoreDict)
    singleScoreDictPath = os.path.join(os.path.split(combined_pinOutputPath)[0], 'singleScoreDict.pickle')
    print('single score dict being saved to disk with path: ' + singleScoreDictPath)
    print('size of score dict: ' + str(len(singleScoreDict.items())))
    with open(singleScoreDictPath, 'wb+') as f:
        pickle.dump(singleScoreDict, f)
    combinedPin.addScores(singleScoreDict, 'NetMHC', '-1')


percolatorOutputDir = os.path.join(tempDir, 'percolatorOutput')
percolatorCommand = [CRUX, 'percolator', '--top-match', str(args.num_matches_per_spectrum), '--tdc', 'F', '--output-dir', percolatorOutputDir, combined_pinOutputPath]
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
        command = ['awk', 'BEGIN {OFS="\t"} NR>1 {print -1, $0}', percolator_target_path]
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

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
parser.add_argument('--minLength', type=int, default=0)
parser.add_argument('--maxLength', type=int, default=0)


args = parser.parse_args()
print('args')
print(args)
if args.allele:
    for x in args.allele:
        assert(x[0] in ['netmhcPrecompute', 'netmhcOnFly', 'MHCFlurryPrecompute', 'MHCFlurryOnFly'])
else:
    args.allele = []
if args.panAlleles:
    for x in args.panAllele:
        assert(x[0] in ['netmhcPanPrecompute', 'netmhcPanOnFly'])
        assert(x[2] ain ['ba', 'elute'])
else:
    args.panAlleles = []
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

outputPath = os.path.join(tempDir, 'search.mzid')
revCatFastaPath = None
pinOutputPath = os.path.join(tempDir, 'search.mzid.pin')
mgf = os.path.join(tempDir, 'spectra.mgf')
shutil.copyfile(args.mgf, mgf)
msgfCommand = ['java', '-Xmx1000M', '-jar', MSGFPLUS,
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
                precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precomputePan.py'), NETMHCPAN, args.additional_proteome, additionalChainPath, additionalPanEluteScoreTablePath, additionalPanBAScoreTablePath, allele[1], str(x), str(THREADS)]
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
                        baseScoreFile = open(os.path.join(args.baseDirectory, args.netmhcPanScoreDir, str(x) + '_top_alleles_elute.scores'), 'rb')
                    elif allele[2] == 'ba':
                        baseScoreFile = open(os.path.join(args.baseDirectory, args.netmhcPanScoreDir, str(x) + '_top_alleles_ba.scores'), 'rb')
                    else:
                        assert(0)
                elif allele[0] == 'netmhcPanOnFly':
                    eluteScoresFile = tempfile.NamedTemporaryFile()
                    baScoresFile = tempfile.NamedTemporaryFile()
                    precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precomputePan.py'), NETMHCPAN, baseFasta, baseChainPath, eluteScoresFile.name, baScoresFile.name, allele[1], str(x), str(THREADS)]
                    proc = subprocess.Popen(precomputeCommand, stdout=subprocess.DEVNULL)
                    outs, errors = proc.communicate()
                    if allele[2] == 'elute':
                        baseScoreFile = eluteScoresFile
                        baScoresFile.close()
                    elif allele[2] == 'ba':
                        baseScoreFile = baScoresFile
                        eluteScoresFile.close()
                    else:
                        assert(0)
                baseScoreTable = ScoreTable.readExisting(baseScoreFile)
            reverse = True if allele[2] == 'elute' else False
            pepToHeader = filterNetMHC(allele[1], x, baseScoreTable, baseChains, baseFasta, additionalScoreTable, additionalChains, args.additional_proteome,  args.rank_filter/100.0, reverse)
            for k,v in pepToHeader.items():
                pepToHeaders[k].update(v)
        for allele in args.allele:
            additionalScoreTable = None
            if args.additional_proteome:
                additionalScoreTablePath = None
                precomputeCommand = None
                if allele[0] in ['netmhcPrecompute', 'netmhcOnFly']:                    
                    precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precompute.py'), NETMHC, args.additional_proteome, additionalChainPath, additionalNetMHCscoreTablePath, allele[1], str(x), str(THREADS)]
                    additionalScoreTablePath = additionalNetMHCscoreTablePath
                elif allele[0] in ['MHCFlurryPrecompute', 'MHCFlurryOnFly']:
                    precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precomputeRemoteFlurry.py'), args.additional_proteome, additionalChainPath, additionalFlurryScoreTablePath, allele[1], str(x)]
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
                    baseScoreFile = open(os.path.join(args.baseDirectory, args.netmhcScoreDir, str(x) + '_top_alleles.scores'), 'rb')
                elif allele[0] == 'MHCFlurryPrecompute':
                    baseScoreFile = open(os.path.join(args.baseDirectory, args.mhcFlurryScoreDir, str(x) + '_top_alleles.scores'), 'rb')
                elif allele[0] == 'netmhcOnFly':
                    baseScoreFile = tempfile.NamedTemporaryFile()
                    precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precompute.py'), NETMHC, baseFasta, baseChainPath, baseScoreFile.name, allele[1], str(x), str(THREADS)]
                    proc = subprocess.Popen(precomputeCommand, stdout=subprocess.DEVNULL)
                    outs, errors = proc.communicate()
                elif allele[0] == 'MHCFlurryOnFly':
                    baseScoreFile = tempfile.NamedTemporaryFile()
                    precomputeCommand = ['python3', os.path.join(PRECOMPUTE_SCRIPTS, 'precomputeRemoteFlurry.py'), baseFasta, baseChainPath, baseScoreFile.name, allele[1], str(x)]
                    proc = subprocess.Popen(precomputeCommand, stdout=subprocess.DEVNULL)
                    outs, errors = proc.communicate()
                else:
                    assert(0)
                baseScoreTable = ScoreTable.readExisting(baseScoreFile)
            reverse = True if allele[0] in ['MHCFlurryPrecompute', 'MHCFlurryOnFly'] else False
            pepToHeader = filterNetMHC(allele[1], x, baseScoreTable, baseChains, baseFasta, additionalScoreTable, additionalChains, args.additional_proteome,  args.rank_filter/100.0, reverse)
            for k,v in pepToHeader.items():
                pepToHeaders[k].update(v)
            

                
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
    fasta = None
    if usingBase:
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
    else:
        fasta = os.path.join(tempDir, 'additional.fasta')
        shutil.copyfile(args.additional_proteome, fasta)
    assert(fasta)
    assert(os.path.isfile(fasta))
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
        scores = [x[0] for x in scorer.scorePeptides(iter(peptides), ['Affinity(nM)'])]
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

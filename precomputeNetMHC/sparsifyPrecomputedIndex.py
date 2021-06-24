import io
import argparse
import os
import pickle
import typing
import sys
import functools

from netMHCCalling import NetMHCScorer, NetMHCRunFailedError




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


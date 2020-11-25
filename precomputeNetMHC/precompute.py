import io
import argparse
import os

parser = argparse.ArgumentParser(description='Precompute NetMHC scores for peptides')
parser.add_argument('fastaInput')
parser.add_argument('output')
parser.add_argument('netMHC')
parser.add_argument('minLength', type=int)
parser.add_argument('maxLength', type=int)

args = parser.parse_args()
assert(args.minLength >= 8)
assert(args.maxLength >= 8)
assert(os.path.exists(args.netMHC))
assert(os.path.exists(args.fastaInput))

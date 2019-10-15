#!/usr/bin/python
import sys
import argparse
import subprocess
import shutil
import csv
import io
import os
import uuid
import zipfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
parser = argparse.ArgumentParser()

parser.add_argument('--unfiltered_archive', type=str)
parser.add_argument('--filtered_archive', type=str)
parser.add_argument('--fdr', type=float)
parser.add_argument('--switch_plot', type=str)
parser.add_argument('--spectra_overlap_plot', type=str)


args = parser.parse_args()
assert(args.unfiltered_archive)
assert(args.filtered_archive)
assert(args.fdr)
assert(args.switch_plot)
assert(args.spectra_overlap_plot)



def get_psms(zip_path, fdr_cutoff):
    with zipfile.ZipFile(zip_path, 'r') as f:
        locations = []
        for x in f.namelist():
            if x.endswith('percolator.target.psms.txt'):
                locations.append(x)
        assert(len(locations) == 1)
        psms = []
        with f.open(locations[0], 'r') as psms_binary_file:
            psms_text_file = io.TextIOWrapper(psms_binary_file)
            psms_reader = csv.DictReader(psms_text_file, delimiter='\t')
            for row in psms_reader:
                scan = row['scan']
                q_val = float(row['percolator q-value'])
                peptide = row['sequence']
                if q_val <= fdr_cutoff:
                    psms.append((scan, peptide))

        return psms

unfiltered_archive_path = args.unfiltered_archive
filtered_archive_path = args.filtered_archive
fdr = args.fdr
switch_plot = args.switch_plot
spectra_overlap_plot = args.spectra_overlap_plot
unfiltered_psms = get_psms(unfiltered_archive_path, fdr)
filtered_psms = get_psms(filtered_archive_path, fdr)

def get_spectra(psms):
    return set([x[0] for x in psms])

#First, how many spectra are shared between the filtered and unfiltered search at the FDR cutoff?
venn2([get_spectra(unfiltered_psms), get_spectra(filtered_psms)], ['Unfiltered Search Spectra', 'Filtered Search Spectra'])
matplotlib.pyplot.savefig(spectra_overlap_plot, format='png')
plt.clf()

#Now, of those spectra shared between the two searches, how many are matched to the same peptide, and how many are matched to a different peptide?
shared_spectra = get_spectra(unfiltered_psms).intersection(get_spectra(filtered_psms))
shared_spectra_same_peptide_count = 0
shared_spectra_different_peptide_count = 0
unfiltered_psm_dict = dict(unfiltered_psms)
filtered_psm_dict = dict(filtered_psms)
for x in list(shared_spectra):
    if unfiltered_psm_dict[x] == filtered_psm_dict[x]:
        shared_spectra_same_peptide_count += 1
    else:
        shared_spectra_different_peptide_count += 1
plt.pie([shared_spectra_same_peptide_count, shared_spectra_different_peptide_count], labels=['Same Peptide', 'Different Peptide'], colors=['gold', 'lightskyblue'])
matplotlib.pyplot.savefig(switch_plot, format='png')



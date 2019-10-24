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
from matplotlib_venn import _math
from matplotlib_venn import venn2_circles
from matplotlib.patches import ConnectionPatch, Circle
from matplotlib.text import Text
import numpy as np
parser = argparse.ArgumentParser()

parser.add_argument('--unfiltered_archive', type=str)
parser.add_argument('--filtered_archive', type=str)
parser.add_argument('--fdr', type=float)
parser.add_argument('--plot', type=str)


args = parser.parse_args()
assert(args.unfiltered_archive)
assert(args.filtered_archive)
assert(args.fdr)
assert(args.plot)

def plot_venn(search_one, search_two, output_location):
    """
    search_one and search_two should both be lists of the form [(spectra, peptide)...]
    """
    spectra_one = set([x[0] for x in search_one])
    assert(len(spectra_one) == len(search_one))
    spectra_two = set([x[0] for x in search_two])
    assert(len(spectra_two) == len(search_two))
    search_one_dict = dict(search_one)
    print('search one dict')
    print(search_one_dict)
    search_two_dict = dict(search_two)
    print('search two dict')
    print(search_two_dict)
    common_spectra = spectra_one.intersection(spectra_two)
    print('common spectra')
    print(common_spectra)
    #the number of spectra shared between the two searches that match against different peptides
    discordant_spectra = 0
    #the number of spectra shared between the two searches that match against the same peptide 
    concordant_spectra = 0
    for spectra in common_spectra:
        if search_one_dict[spectra] == search_two_dict[spectra]:
            concordant_spectra += 1
        else:
            discordant_spectra += 1
            
    circles = venn2_circles([spectra_one, spectra_two])
    sorted_circles = sorted(circles, key=lambda x: x.center[0])
    bigger_circle = max(circles, key=lambda x: x.radius)
    bigger_radius = bigger_circle.radius
    left_point = np.array([sorted_circles[0].center[0] - sorted_circles[0].radius, sorted_circles[0].center[1]])
    right_point = np.array([sorted_circles[1].center[0] + sorted_circles[1].radius, sorted_circles[1].center[1]])
    left_intersection = max(_math.circle_line_intersection(sorted_circles[0].center, sorted_circles[0].radius, left_point, right_point), key=lambda x: x[0])
    right_intersection = min(_math.circle_line_intersection(sorted_circles[1].center, sorted_circles[1].radius, left_point, right_point), key=lambda x: x[0])
    line = ConnectionPatch(left_intersection, right_intersection, 'data', 'data')
    plt.gca().add_patch(line)
    print(sorted_circles[0].center)
    print(sorted_circles[1].center)
    circle_intersections = _math.circle_circle_intersection(sorted_circles[0].center, sorted_circles[0].radius, sorted_circles[1].center, sorted_circles[1].radius)
    upper_circle_intersection = max(circle_intersections, key=lambda x: x[1])
    #take the centroid
    upper_text_location = (left_intersection + right_intersection + upper_circle_intersection)/3.0
    #plt.rc('text', usetex=True)
    plt.text(upper_text_location[0], upper_text_location[1], str(concordant_spectra) + '\n' + r'$p_i = p_j$')
    lower_circle_intersection = min(circle_intersections, key=lambda x: x[1])
    lower_text_location = (left_intersection + right_intersection + lower_circle_intersection)/3.0
    plt.text(lower_text_location[0], lower_text_location[1], str(discordant_spectra) + '\n' + r'$p_i \neq p_j$')
    venn_diagram = venn2([spectra_one, spectra_two], ['Unfiltered', 'Filtered'])
    venn_diagram.get_label_by_id('11').set_text('')
    matplotlib.pyplot.savefig(output_location, format='png')



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
plot = args.plot
unfiltered_psms = get_psms(unfiltered_archive_path, fdr)
filtered_psms = get_psms(filtered_archive_path, fdr)

def get_spectra(psms):
    return set([x[0] for x in psms])


plot_venn(unfiltered_psms, filtered_psms, plot)



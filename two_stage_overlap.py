#!/usr/bin/python

import argparse
import matplotlib
#matplotlib.use('Agg')
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import _math
from matplotlib_venn import venn2_circles
from matplotlib.patches import ConnectionPatch, Circle
from matplotlib.text import Text
parser = argparse.ArgumentParser()
parser.add_argument('--search_one', type=str)
parser.add_argument('--search_two', type=str)
parser.add_argument('--output', type=str)


def formatter(s):
    print('s: %d' %s)
    return ''

def plot_venn(search_one, search_two):
    """
    search_one and search_two should both be lists of the form [(spectra, peptide)...]
    """
    spectra_one = set([x[0] for x in search_one])
    assert(len(spectra_one) == len(search_one))
    spectra_two = set([x[0] for x in search_two])
    assert(len(spectra_two) == len(search_two))
    search_one_dict = dict(search_one)
    search_two_dict = dict(search_two)
    common_spectra = spectra_one.intersection(spectra_two)
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
    print('circles')
    print(circles[0].center[0])
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
    venn_diagram = venn2([spectra_one, spectra_two], ['one', 'two'])
    venn_diagram.get_label_by_id('11').set_text('')
    plt.show()
    return True
"""    line_middle = np.array([(left_intersection[0] + right_intersection[0])/2, left_intersection[1]])
    #consider a vertical line positioned at the center of the divider line. See where it intersects the two circles.
    vertical_line_upper_point = np.array([line_middle[0], line_middle[1] + bigger_radius])
    vertical_line_lower_point = np.array([line_middle[0], line_middle[1] - bigger_radius])
    left_circle_vertical_intersections = _math.circle_line_intersection(sorted_circles[0].center, sorted_circles[0].radius, v0
    vertical_line_patch = ConnectionPatch(vertical_line_upper_point, vertical_line_lower_point, 'data', 'data')
    venn2([spectra_one, spectra_two], ['one', 'two'])
    plt.gca().add_patch(line)
    plt.gca().add_patch(vertical_line_patch)
    plt.show()
    return True"""
    #upper_intersection = min([_math.circle_line_intersection(x.center, x.radius, np.array(
    #circle_intersections = _math.circle_circle_intersection(sorted_circles[0].center, sorted_circles[1].center

plot_venn([(1, 'A'), (2, 'B'), (3, 'C'), (4, 'D')], [(4, 'D'), (5, 'E'), (6, 'F')])
sys.exit()


args = parser.parse_args()
one = set()
two = set()
with open(args.sample_one, 'r') as f:
    for x in f:
        line = x.strip()
        if len(line) > 0:
            one.add(line)
with open(args.sample_two, 'r') as f:
    for x in f:
        line = x.strip()
        if len(line) > 0:
            two.add(line)

venn2([one, two], [args.sample_one_name, args.sample_two_name])
matplotlib.pyplot.savefig(args.output, format='png')

with open(args.union, 'w') as f:
    for x in list(one.union(two)):
        f.write(x + '\n')
with open(args.intersection, 'w') as f:
    for x in list(one.intersection(two)):
        f.write(x + '\n')
with open(args.sample_one_unique, 'w') as f:
    for x in list(one):
        if x not in two:
            f.write(x + '\n')
with open(args.sample_two_unique, 'w') as f:
    for x in list(two):
        if x not in one:
            f.write(x + '\n')

import Histogram

import numpy as np
import matplotlib
# Don't require an X-Server
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import csv

def plot_histogram(plot_file, histogram, width, tree_name, d, t, l, max_x=None, max_y=None, title=True):
    # Set the max limits
    if max_y is not None:
        plt.ylim(top=float(max_y))
    if max_x is not None:
        plt.xlim(right=float(max_x))
    plt.bar(histogram.keys(), histogram.values(), width)
    # Force y-axis to use scientific notation
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    # Find the exponent in order to put it in the ylabel
    # Force offset text to update using draw
    #TODO: there MUST be a better way to do this
    plt.draw()
    ax = plt.gca()
    # matplotlib sure is intuitive and easy to use!
    # Get the exponent text from the y-axis and format it into latex
    exponent_text = ax.get_yaxis().get_offset_text().get_text()
    exponent = float(exponent_text.split("e")[-1])
    latex_exponent = r"x$10^{%d}$" % exponent
    # Don't display it because we're going to use it in the y-axis label
    ax.yaxis.offsetText.set_visible(False)
    # Set the labels
    plt.xlabel("Distance", fontsize=18)
    plt.ylabel("Number of MPR Pairs {}".format(latex_exponent), fontsize=18)
    # y=1.08 is a hack to make the title display above
    if title:
        plt.title("{} with costs D:{}, T:{}, L:{}".format(tree_name, d, t, l), y=1.08, fontsize=18)
    plt.savefig(plot_file, bbox_inches='tight')
    plt.clf()

def csv_histogram(csv_file, histogram):
    with open(csv_file, 'w') as csv_handle:
        writer = csv.writer(csv_handle)
        for key, value in histogram.items():
            writer.writerow([key, value])

def normalize_xvals(histogram):
    max_xval = float(max(histogram.keys()))
    if max_xval == 0:
        return histogram
    new_hist = { k/max_xval : v for k,v in histogram.items() }
    return new_hist

def normalize_yvals(histogram):
    total_yval = float(sum(histogram.values()))
    if total_yval == 0:
        return histogram
    new_hist = { k : v/total_yval for k,v in histogram.items() }
    return new_hist

def cumulative(histogram):
    total = 0
    new_hist = {}
    for k, v in histogram.items():
        new_hist[k] = v + total
        total += v
    return new_hist

def omit_zeros(histogram):
    return { k : v for k,v in histogram.items() if k != 0 }

def compute_stats(histogram):
    # Diameter of MPR-space
    # Average distance between MPRs and standard deviation
    diameter = max(histogram.keys())
    # Re-convert to a Histogram to calculate stats
    h = Histogram.Histogram(histogram)
    mean = h.mean()
    std = h.standard_deviation()
    return diameter, mean, std


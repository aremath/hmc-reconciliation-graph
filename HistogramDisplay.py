import numpy as np
import matplotlib
# Don't require an X-Server
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import csv

def plot_histogram(plot_file, histogram, width, tree_name, d, t, l):
    plt.bar(histogram.keys(), histogram.values(), width)
    # Force y-axis to use scientific notation
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    # Find the exponent in order to put it in the ylabel
    # Force offset text to update using draw
    #TODO: there MUST be a better way to do this
    plt.draw()
    ax = plt.gca()
    # matplotlib sure is intuitive and easy to use!
    exponent_text = ax.get_yaxis().get_offset_text().get_text()
    exponent = float(exponent_text.split("e")[-1])
    latex_exponent = r"x$10^{%d}$" % exponent
    # Don't display it because we're going to use it in the y-axis
    ax.yaxis.offsetText.set_visible(False)
    # Set the labels
    plt.xlabel("Distance", fontsize=18)
    plt.ylabel("Number of MPR Pairs {}".format(latex_exponent), fontsize=18)
    # y=1.08 is a hack to make the title display above 
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
    #TODO: use the histogram methods from histogram.py
    # Flattening is a bad idea because this array might be exponentially large.
    flat_h = flatten(histogram)
    flat_h = np.asarray(flat_h)
    mean = flat_h.mean()
    std = flat_h.std()
    return diameter, mean, std

def flatten(histogram):
    l = []
    for k, v in histogram.keys():
        for i in range(v):
            l.append(v)
    return l


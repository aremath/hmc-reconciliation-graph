import numpy as np
import matplotlib
# Don't require an X-Server
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv

def plot_histogram(plot_file, histogram, width, tree_name):
    plt.bar(histogram.keys(), histogram.values(), width)
    plt.xlabel("Distance")
    plt.ylabel("Number of MPR Pairs")
    # y=1.08 is a hack to make the title display above 
    plt.title("Pairwise Distance Vector for {}".format(tree_name), y=1.08)
    plt.savefig(plot_file)
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


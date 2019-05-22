import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot_histogram(plot_file, histogram, x_norm, y_norm):
    if x_norm:
        max_xval = float(max(histogram.keys()))
        width = 1/max_xval
        histogram = normalize_xvals(histogram)
    else:
        width = 1
    if y_norm:
        histogram = normalize_yvals(histogram)
    plt.bar(histogram.keys(), histogram.values(), width)
    # plt.xlabel
    # plt.ylabel
    # plt.title
    plt.savefig(plot_file)
    plt.clf()

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


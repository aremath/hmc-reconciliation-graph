import ClusterUtil
import ReconciliationVisualization as RV
import HistogramDisplay

import argparse
import signal
from pathlib import Path
import collections
import csv
import pickle
import itertools
import time

import numpy as np
import scipy.stats
import matplotlib
# Don't require an X-Server
matplotlib.use("Agg")
import matplotlib.pyplot as plt

#TODO: if args.input is specified, then no need for d,t,l,input,etc.
def process_args():
    # Required arguments - input file, D T L costs
    parser = argparse.ArgumentParser("")
    parser.add_argument("--input", metavar="<filename>", required=True,
        help="The path to a folder of .newick files.")
    parser.add_argument("-d", type=int, metavar="<duplication_cost>", required=True,
        help="The relative cost of a duplication.")
    parser.add_argument("-t", type=int, metavar="<transfer_cost>", required=True,
        help="The relative cost of a transfer.")
    parser.add_argument("-l", type=int, metavar="<loss_cost>", required=True,
        help="The relative cost of a loss.")
    parser.add_argument("-k", type=int, metavar="<number_of_clusters>", required=True,
        help="How many clusters to create.")
    parser.add_argument("--output", metavar="<filename>", required=False,
        help="The path to a file which will store the scores.")
    parser.add_argument("--absolute", action="store_true", required=False,
        help="Use absolute improvement or relative improvement")
    depth_or_n = parser.add_mutually_exclusive_group(required=True)
    depth_or_n.add_argument("--depth", type=int, metavar="<tree_depth>",
        help="How far down the graph to consider event splits.")
    depth_or_n.add_argument("--nmprs", type=int, metavar="<tree_depth>",
        help="How many MPRs to consider")
    which_plot = parser.add_mutually_exclusive_group(required=True)
    which_plot.add_argument("--s1s2", action="store_true",
        help="Correlate the improvements across multiple scores.")
    which_plot.add_argument("--ki", action="store_true",
        help="Correlate the improvement with the number of clusters.")
    which_plot.add_argument("--ni", action="store_true",
        help="Correlate the improvement with the number of splits used.")
    which_plot.add_argument("--mpri", action="store_true",
        help="Correlate the improvement with the number of MPRs.")
    which_plot.add_argument("--time", action="store_true",
        help="Correlate the running time with the number of MPRs.")
    score = parser.add_mutually_exclusive_group(required=True)
    score.add_argument("--pdv", action="store_true",
        help="Use the weighted average distance to evaluate clusters.")
    score.add_argument("--support", action="store_true",
        help="Use the weighted average event support to evaluate clusters.")
    args = parser.parse_args()
    # Raw improvement -- no option since you probably don't want it.
    args.raw = False
    return args

class TimeoutError(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutError

def choose_score(args):
    if args.pdv:
        return ClusterUtil.mk_pdv_score
    elif args.support:
        return ClusterUtil.mk_support_score
    else:
        assert False

def choose_name(args):
    if args.pdv:
        return "WAD"
    elif args.support:
        return "WAS"
    else:
        assert False

def mk_title(title, args):
    #return title + " ({})".format(choose_name(args))
    if args.pdv:
        return "Weighted Average Distance"
    elif args.support:
        return "Weighted Average Support"
    else:
        assert False

def mk_savename(savename, args):
    return savename + "_{}.pdf".format(choose_name(args))

def choose_imp(args):
    if args.pdv:
        return ClusterUtil.calc_improvement_pdv
    elif args.support:
        return ClusterUtil.calc_improvement_support
    else:
        assert False

# (s0, s1), (s1, s2), ...
def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)

def get_tree_files(pathstr):
    p = Path(pathstr)
    all_files = [f for f in p.glob("**/*") if f.is_file()]
    tree_files = [f for f in all_files if f.suffix==".newick"]
    return tree_files

# Do the clustering with a timeout
def timeout_cluster(recon_g, gene_root, score, mpr_count, args, timeout, max_splits):
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(timeout)
    try:
        c = cluster(recon_g, gene_root, score, mpr_count, args, max_splits)
    except TimeoutError:
        return None
    signal.alarm(0)
    return c

# Actually perform the clustering
def cluster(recon_g, gene_root, score, mpr_count, args, max_splits):
    if args.depth is not None:
        return ClusterUtil.cluster_graph(recon_g, gene_root, score, args.depth, args.k, max_splits)
    elif args.nmprs is not None:
        return ClusterUtil.cluster_graph_n(recon_g, gene_root, score, args.nmprs, mpr_count, args.k, max_splits)
    else:
        assert False

# # # #
# Get timing info
# # # #

def get_times(tree_files, mk_score, args, timeout=1200, min_mprs=1000, max_splits=200):
    nmprs = []
    times = []
    ftrees = []
    n_timed_out = 0
    n = len(tree_files)
    for (i, f) in enumerate(tree_files):
        print("{}: {}/{}".format(f, i+1, n))
        # Get the recon graph + other info
        gene_tree, species_tree, gene_root, recon_g, mpr_count = \
            ClusterUtil.get_tree_info(str(f), args.d,args.t,args.l)
        # Only care about trees with a certain number of MPRs
        if mpr_count < min_mprs:
            continue
        score = mk_score(species_tree, gene_tree, gene_root)
        start = time.time()
        t = timeout_cluster(recon_g, gene_root, score, mpr_count, args, timeout, max_splits)
        end = time.time()
        elapsed = end - start
        if t is None:
            print("{} timed out".format(f))
            n_timed_out += 1
            continue
        nmprs.append(mpr_count)
        times.append(elapsed)
        ftrees.append(f)
    return (nmprs, times, n_timed_out), ftrees

def analyze_times(data, args):
    nmprs, times, n_timed_out = data
    time_mean = np.mean(times)
    time_std = np.std(times)
    print("Time:")
    print("Mean: {}".format(time_mean))
    print("Standard Deviation: {}".format(time_std))
    print("Timeouts: {}".format(n_timed_out))

    # Also plot it with number of MPRs
    plt.scatter(nmprs, times, alpha=0.5)
    plt.ylim((0, 1200))
    title = "Time Elapsed by Number of MPRs"
    plt.xlabel("Number of MPRs")
    plt.xscale("log")
    plt.ylabel("Running Time (seconds)")
    plt.title(mk_title(title, args))
    plt.savefig(mk_savename("nmpr_time_plot", args), bbox_inches="tight")
    plt.clf()

def get_time_data(trees, args):
    assert args.k == 1
    mk_score = choose_score(args)
    return get_times(trees, mk_score, args)

# # # #
# Correlate improvement with number of mprs
# # # #

def get_scores_mprs(tree_files, mk_score, args, timeout=1200, min_mprs=1000, max_splits=200):
    nmprs = []
    scores = []
    ftrees = []
    n = len(tree_files)
    for (i, f) in enumerate(tree_files):
        print("{}: {}/{}".format(f, i+1, n))
        # Get the recon graph + other info
        gene_tree, species_tree, gene_root, recon_g, mpr_count = \
            ClusterUtil.get_tree_info(str(f), args.d,args.t,args.l)
        # Only care about trees with a certain number of MPRs
        if mpr_count < min_mprs:
            continue
        score = mk_score(species_tree, gene_tree, gene_root)
        t = timeout_cluster(recon_g, gene_root, score, mpr_count, args, timeout, max_splits)
        if t is None:
            print("{} timed out".format(f))
            continue
        _,tree_scores = t
        nmprs.append(mpr_count)
        scores.append((tree_scores[0], tree_scores[1]))
        ftrees.append(f)
    return (nmprs, scores), ftrees

def plot_mprs_improvement(data, args):
    nmprs, scores = data
    imp = choose_imp(args)
    improvements = [imp(s[1], s[0]) for s in scores]
    # Fixes error when nmprs is python L (since numpy has no 128 bit intergers)
    np_nmprs = np.array(nmprs, dtype="float64")
    slope, intercept, r, p, stderr = scipy.stats.linregress(np_nmprs, improvements)
    print("R: {}, P: {}, E: {}, N: {}".format(r, p, stderr, len(nmprs)))
    plt.scatter(nmprs, improvements, alpha=0.5)
    plt.ylim((0.95, 1.6))
    title = "Absolute Improvement by Number of MPRs"
    plt.xlabel("Number of MPRs")
    plt.xscale("log")
    plt.ylabel("Improvement")
    plt.title(mk_title(title, args))
    plt.savefig(mk_savename("nmpr_improvement_plot", args), bbox_inches="tight")
    plt.clf()

def get_mpri_data(trees, args):
    assert args.k == 1
    mk_score = choose_score(args)
    return get_scores_mprs(trees, mk_score, args)

# # # #
# Correlate improvement with number of clusters
# # # #

def get_scores(tree_files, mk_score, args, timeout=1200, min_mprs=1000, max_splits=200):
    scores = []
    ftrees = []
    n = len(tree_files)
    for (i, f) in enumerate(tree_files):
        print("{}: {}/{}".format(f, i+1, n))
        # Get the recon graph + other info
        gene_tree, species_tree, gene_root, recon_g, mpr_count = \
            ClusterUtil.get_tree_info(str(f), args.d,args.t,args.l)
        # Only care about trees with a certain number of MPRs
        if mpr_count < min_mprs:
            continue
        score = mk_score(species_tree, gene_tree, gene_root)
        t = timeout_cluster(recon_g, gene_root, score, mpr_count, args, timeout, max_splits)
        if t is None:
            print("{} timed out".format(f))
            continue
        _,tree_scores = t
        scores.append(tree_scores)
        ftrees.append(f)
    return scores, ftrees

def get_top_n(scores, trees, n, imp):
    ztrees = zip(trees, scores)
    ztrees = sorted(ztrees, key=lambda x: imp(x[1][1], x[1][0]))
    return [x[0] for x in ztrees[-n:]]

# Scores relative to the previous score
def scores_to_improvements_relative(scores, imp):
    return [[imp(big,small) for small,big in pairwise(s)] for s in scores]

# Scores relative to the score for one cluster
def scores_to_improvements_absolute(scores, imp):
    # The score for one cluster
    improvements = []
    for s in scores:
        s_improvement = []
        small = s[0]
        # Want to ignore 1 since the improvement will always be 1.0
        for big in s[1:]:
            s_improvement.append(imp(big,small))
        improvements.append(s_improvement)
    return improvements

# Just the raw weighted average scores
def scores_to_improvements_raw(scores):
    return [[-1*s for s in series] for series in scores]

def plot_k_improvement(scores, args):
    imp = choose_imp(args)
    if args.absolute:
        improvements = scores_to_improvements_absolute(scores, imp)
        title = "Improvement by Number of Clusters"
        yax = "Improvement"
        xadj = 1
    elif args.raw:
        improvements = scores_to_improvements_raw(scores)
        title = "Weighted Average Score by Number of Clusters"
        yax = "Score"
        xadj = 0
    else:
        improvements = scores_to_improvements_relative(scores, imp)
        title = "Relative Improvement by Number of Clusters"
        yax = "Improvement"
        xadj = 1
    for series in improvements:
        xs = []
        ys = []
        for k, imp in enumerate(series):
            xs.append(k + args.k + xadj)
            ys.append(imp)
        plt.plot(xs,ys, c="blue", marker="o", alpha=0.1)
    plt.ylim((0.95,2.4))
    plt.xlim((1.8,12))
    plt.xlabel("K")
    plt.ylabel(yax)
    plt.title(mk_title(title, args))
    plt.savefig(mk_savename("k_improvement_plot", args), bbox_inches="tight")
    plt.clf()

def get_ki_data(trees, args):
    assert args.k == 1
    mk_score = choose_score(args)
    return get_scores(trees, mk_score, args)

# # # #
# Correlate improvement with both metrics
# # # #

# Can't pickle a lambda, so define this here
def mk_default_list():
    return collections.defaultdict(list)

# 1200 s = 20 minutes
def get_improvements(tree_files, cluster_mk_scores, eval_mk_scores, args, timeout=1200, min_mprs=1000, max_splits=200):
    # Key: clustering method index. Value: list of trees that finished
    ftrees = collections.defaultdict(list)
    # Keys: clustering method index, evaluation method index. Value: list of improvements
    improvements = collections.defaultdict(mk_default_list)
    n = len(tree_files)
    for (i, f) in enumerate(tree_files):
        print("{}: {}/{}".format(f, i+1, n))
        # Get the recon graph + other info
        gene_tree, species_tree, gene_root, recon_g, mpr_count = \
            ClusterUtil.get_tree_info(str(f), args.d,args.t,args.l)
        # Only care about trees with a certain number of MPRs
        if mpr_count < min_mprs:
            continue
        # Get all the scoring functions ready for this tree
        cluster_scores = [mk_cs(species_tree, gene_tree, gene_root) for mk_cs in cluster_mk_scores]
        eval_scores = [mk_es(species_tree, gene_tree, gene_root) for mk_es in eval_mk_scores]
        # Counts MPRs to weight scores
        mpr_counter = ClusterUtil.mk_count_mprs(gene_root)
        # Evaluate the original graph on each eval metric to record improvement
        one_scores = [eval_s(recon_g) for eval_s in eval_scores]
        # Perform the clustering for each cluster score
        for i1, cluster_score in enumerate(cluster_scores):
            t = timeout_cluster(recon_g, gene_root, cluster_score, mpr_count, args, timeout, max_splits)
            if t is None:
                print("{} timed out".format(f))
                continue
            ftrees[i1].append(f)
            graphs,_ = t
            # Evaluate the clustering for each evaluation score
            for i2, eval_score in enumerate(eval_scores):
                one_score = one_scores[i2]
                k_score = ClusterUtil.get_score_nodp(graphs, eval_score, mpr_counter)
                improvement = ClusterUtil.calc_improvement(k_score, one_score)
                improvements[i1][i2].append(improvement)
    return improvements, ftrees

# Make plottable xs and ys for each evaluation metric
# Specific to two dimensional output
def transform_s1_s2(improvements):
    series = []
    for i, evals in enumerate(improvements.values()):
        xs = evals[0]
        ys = evals[1]
        series.append((xs, ys))
    return series

def plot_s1_s2(improvements, args):
    series = transform_s1_s2(improvements)
    colors=["red", "blue"]
    for i,s in enumerate(series):
        xs, ys = s
        new_xs = [1.0/x for x in xs]
        color = colors[i]
        plt.scatter(new_xs, ys, c=color, alpha=0.5)
        slope, intercept, r, p, stderr = scipy.stats.linregress(new_xs, ys)
        print("R: {}, P: {}, E: {}, N: {}".format(r, p, stderr, len(new_xs)))
    #plt.ylim((0, 1))
    #plt.xlim((0, 1))
    plt.xlabel("WAD improvement")
    plt.ylabel("WAS improvement")
    plt.title(mk_title("Relative Improvement for Two Clusters", args))
    plt.savefig(mk_savename("s1_s2_plot", args), bbox_inches="tight")
    plt.clf()

def get_s1_s2_data(trees, args):
    #cluster_names = ["PDV", "Event Support"]
    cluster_mk_scores = [choose_score(args)]
    #cluster_mk_scores = [ClusterUtil.mk_pdv_score, ClusterUtil.mk_support_score]
    eval_mk_scores = [ClusterUtil.mk_pdv_score, ClusterUtil.mk_support_score]
    return get_improvements(trees, cluster_mk_scores, eval_mk_scores, args)

# # # #
# Correlate improvement with the number of splits used
# # # #

def get_n_improvements(tree_files, mk_score, args, timeout=1200, min_mprs=1000, max_splits=200):
    imp = choose_imp(args)
    series = []
    ftrees = []
    n = len(tree_files)
    for (i, f) in enumerate(tree_files):
        print("{}: {}/{}".format(f, i+1, n))
        # Get the recon graph + other info
        gene_tree, species_tree, gene_root, recon_g, mpr_count = \
            ClusterUtil.get_tree_info(str(f), args.d,args.t,args.l)
        #print("MPR count: {}".format(mpr_count))
        # Only care about trees with a certain number of MPRs
        if mpr_count < min_mprs:
            continue
        mpr_counter = ClusterUtil.mk_count_mprs(gene_root)
        score = mk_score(species_tree, gene_tree, gene_root)
        xs = []
        ys = []
        old_ngs = 0
        # Try multiple values of n
        # from 2 to 128
        for n_thres in [2**i for i in range(1,8)]:
            args.nmprs = n_thres
            
            # Timeout
            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(timeout)
            try:
                gs = ClusterUtil.full_split_n(recon_g, gene_root, args.nmprs, mpr_count)
                print("Number of splits: {}".format(len(gs)))
                # Don't bother re-clustering if we split to get the same number of gs (induces the same split)
                if len(gs) == old_ngs or len(gs) > max_splits:
                    continue
                else:
                    old_ngs = len(gs)
                    _, scores = ClusterUtil.combine(gs, score, args.k, mpr_counter)
            except TimeoutError:
                print("{} timed out".format(f))
                continue
            signal.alarm(0)

            true_n = len(scores)
            # Compare two clusters to one cluster
            two_s = scores[1]
            one_s = scores[0]
            improvement = imp(two_s, one_s)
            xs.append(true_n)
            ys.append(improvement)
        ftrees.append(f)
        series.append((xs[:], ys[:]))
    return series, ftrees

def plot_n_improvement(series, args):
    # Different color for each series
    n = len(series)
    norm = matplotlib.colors.Normalize(0, n+1)
    cmap = matplotlib.cm.get_cmap("hsv", n+1)
    for i, (xs, ys) in enumerate(series):
        if args.pdv:
            new_ys = [1.0/y for y in ys]
        else:
            new_ys = ys
        plt.plot(xs, new_ys, c=cmap(norm(i)), marker="o", alpha=0.5)
    plt.ylim((0.95, 1.6))
    plt.xlim((0, 200))
    #plt.ylim(bottom=1)
    #plt.xlim((0,500))
    plt.xlabel("N")
    plt.ylabel("Improvement")
    plt.title(mk_title("Improvement Correlated with Number of Splits", args))
    plt.savefig(mk_savename("n_improvement_plot", args), bbox_inches="tight")
    plt.clf()

def max_change(l):
    max_change = 0
    for l1, l2 in itertools.product(l, repeat=2):
        change = abs(l1 - l2)
        if change > max_change:
            max_change = change
    return max_change

def is_flatline(ys):
    if len(ys) < 2:
        return False
    lastval = ys[-1]
    if ys[-2] == ys[-1]:
        return True

def plot_ni_changes(series, args):
    final_xs = []
    final_ys = []
    max_x = max([max(s[0]) for s in series])
    for threshold in range(max_x):
        # Find the number of changed elements
        flatlined = 0
        for xs, ys in series:
            # All the y-values where the x-value is less than threshold
            good_ys = [ys[i] for i,x in enumerate(xs) if x <= threshold]
            # There are some points after the threshold. Did it reach the flatline before?
            if len(good_ys) > 0 and len(ys) > len(good_ys):
                if is_flatline(ys) and good_ys[-1] == ys[-1]:
                    flatlined += 1
            # All of the points are below the threshold. Did it flatline?
            elif len(good_ys) > 0:
                if is_flatline(ys):
                    flatlined += 1
        # Compute the ratio among families that could have changed.
        final_xs.append(threshold)
        # Ignore families with only one data point
        n = len([s for s in series if len(s) > 1])
        final_ys.append(float(flatlined) / n)
    # Plot it!
    plt.plot(final_xs, final_ys)
    # The ratio is between 0 and 1
    plt.ylim(0,1)
    plt.xlabel("N")
    plt.ylabel("Fraction of Flatlined Families")
    plt.title(mk_title("Flatlined Families as a Function of N", args))
    plt.savefig(mk_savename("ni_change_plot", args), bbox_inches="tight")
    plt.clf()

def plot_ni_variance(series, args):
    final_xs = []
    final_ys = []
    max_x = max([max(s[0]) for s in series])
    for threshold in range(max_x):
        variances = []
        # Find the variance after the threshold for each series
        for xs, ys in series:
            # All the y-values where the x-value is greater than threshold
            good_ys = [ys[i] for i,x in enumerate(xs) if x >= threshold]
            #print(len(good_ys))
            #if len(good_ys) > 0:
            if len(good_ys) > 1:
                variances.append(max_change(good_ys))
                #variances.append(np.var(good_ys))
        # Average them and keep the data
        if len(variances) > 0:
            final_xs.append(threshold)
            final_ys.append(np.mean(variances))
    # Plot it!
    plt.plot(final_xs, final_ys)
    plt.xlabel("N")
    #plt.ylabel("Maximum Variance")
    #plt.title("Maximum Variance of Improvement")
    plt.ylabel("Average of Max Difference of Improvement")
    plt.title("")
    plt.savefig(mk_savename("ni_variance_plot.pdf", args), bbox_inches="tight")
    plt.clf()

def get_ni_data(trees, args):
    assert args.k == 1
    mk_score = choose_score(args)
    return get_n_improvements(trees, mk_score, args)

#MAIN

def main():
    args = process_args()
    p = Path(args.input)
    # If the input is a file, unpickle and use it as data
    if p.is_file():
        with open(str(p), "r") as infile:
            ftrees, data = pickle.load(infile)
    # Otherwise, calculate the data from the trees in the input folder
    else:
        trees = get_tree_files(args.input)
        # Each alg returns some data and a list of tree files that finished
        # We want to keep the list of tree files too,
        # so that we can associate each datum with a single tree.
        if args.s1s2:
            data, ftrees = get_s1_s2_data(trees, args)
        if args.ki:
            data, ftrees = get_ki_data(trees, args)
        if args.ni:
            data, ftrees = get_ni_data(trees, args)
        if args.mpri:
            data, ftrees = get_mpri_data(trees, args)
        if args.time:
            data, ftrees = get_time_data(trees, args)
    # Now generate pretty plots
    if args.s1s2:
        plot_s1_s2(data, args)
    if args.ki:
        plot_k_improvement(data, args)
    if args.ni:
        plot_n_improvement(data, args)
        plot_ni_variance(data, args)
        plot_ni_changes(data, args)
    if args.mpri:
        plot_mprs_improvement(data, args)
    if args.time:
        analyze_times(data, args)
    # Dump the data to the output file
    # NOTE: plotting functions should not alter the data
    if args.output is not None:
        with open(args.output, "w") as outfile:
            pickle.dump((ftrees, data), outfile)

if __name__ == "__main__":
    main()


import ClusterUtil
import ReconciliationVisualization as RV
import HistogramDisplay

import argparse
import signal
from pathlib import Path
import collections

import numpy as np
import matplotlib
# Don't require an X-Server
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
    depth_or_n = parser.add_mutually_exclusive_group(required=True)
    depth_or_n.add_argument("--depth", type=int, metavar="<tree_depth>",
        help="How far down the graph to consider event splits.")
    depth_or_n.add_argument("--nmprs", type=int, metavar="<tree_depth>",
        help="How many MPRs to consider")
    which_plot = parser.add_mutually_exclusive_group(required=True)
    which_plot.add_argument("--d1d2", action="store_true",
        help="Correlate the improvements across distance metrics.")
    which_plot.add_argument("--ki", action="store_true",
        help="Correlate the improvement with the number of clusters.")
    which_plot.add_argument("--ni", action="store_true",
        help="Correlate the improvement with the number of splits used.")
    args = parser.parse_args()
    return args

class TimeoutError(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutError

def get_tree_files(pathstr):
    p = Path(pathstr)
    all_files = [f for f in p.glob("**/*") if f.is_file()]
    tree_files = [f for f in all_files if f.suffix==".newick"]
    return tree_files

# Do the clustering with a timeout
def timeout_cluster(recon_g, gene_root, d, mpr_count, args, timeout):
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(timeout)
    try:
        graphs,scores = cluster(recon_g, gene_root, d, mpr_count, args)
    except TimeoutError:
        return None
    signal.alarm(0)
    return graphs, scores

# Actually perform the clustering
def cluster(recon_g, gene_root, d, mpr_count, args):
    if args.depth is not None:
        graphs, scores = ClusterUtil.cluster_graph(recon_g, gene_root, d, args.depth, args.k)
    elif args.nmprs is not None:
        graphs, scores = ClusterUtil.cluster_graph_n(recon_g, gene_root, d, args.nmprs, mpr_count, args.k)
    else:
        assert False
    return graphs, scores

# # # #
# Correlate improvement with number of clusters
# # # #

def get_scores(tree_files, metric, args, timeout=1200, min_mprs=1000):
    scores = []
    n = len(tree_files)
    for (i, f) in enumerate(tree_files):
        print("{}: {}/{}".format(f, i, n))
        # Get the recon graph + other info
        gene_tree, species_tree, gene_root, recon_g, mpr_count = \
            ClusterUtil.get_tree_info(str(f), args.d,args.t,args.l)
        # Only care about trees with a certain number of MPRs
        if mpr_count < min_mprs:
            continue
        d = metric(species_tree, gene_tree, gene_root)
        t = timeout_cluster(recon_g, gene_root, d, mpr_count, args, timeout)
        if t is None:
            print("{} timed out".format(f))
            continue
        _,tree_scores = t
        scores.append(tree_scores)
    return scores

def scores_to_improvements(scores):
    return [1 - new/float(old) for new,old in scores]

# Make plottable xs and ys from a list of improvements for each tree file
def transform_k_improvement(improvements, initial_k):
    xs = []
    ys = []
    for imps in improvements:
        for k, imp in enumerate(imps):
            # The first index of imps is for k clusters as per the comment
            # in the bottom of combine in ClusterUtil
            xs.append(k + initial_k)
            ys.append(imp)
    return xs, ys

def plot_k_improvement(improvements, initial_k):
    xs, ys = transform_k_improvement(improvements, initial_k)
    plt.scatter(xs, ys, c="blue", alpha=0.5)
    # Improvement is bounded between 0 and 1
    #plt.ylim((0,1))
    plt.xlim((0,10))
    plt.xlabel("Number of Clusters")
    plt.ylabel("Improvement")
    plt.title("Relative improvement")
    plt.savefig("k_improvement_plot.pdf", bbox_inches="tight")
    plt.clf()

def k_improvement_plot(trees, args):
    assert args.k == 1
    metric = ClusterUtil.mk_support_dist
    scores = get_scores(trees, metric, args)
    improvements = [scores_to_improvements(s) for s in scores]
    plot_k_improvement(improvements, args.k)

# # # #
# Correlate improvement with both metrics
# # # #

# 1200 s = 20 minutes
def get_improvements(tree_files, cluster_metrics, eval_metrics, args, timeout=1200, min_mprs=1000):
    # Keys clustering method index, evaluation method index, list of improvements
    improvements = collections.defaultdict(lambda: collections.defaultdict(list))
    n = len(tree_files)
    for (i, f) in enumerate(tree_files):
        print("{}: {}/{}".format(f, i, n))
        # Get the recon graph + other info
        gene_tree, species_tree, gene_root, recon_g, mpr_count = \
            ClusterUtil.get_tree_info(str(f), args.d,args.t,args.l)
        # Only care about trees with a certain number of MPRs
        if mpr_count < min_mprs:
            continue
        # Get all the distances ready
        cluster_ds = [mk_cd(species_tree, gene_tree, gene_root) for mk_cd in cluster_metrics]
        eval_ds = [mk_ed(species_tree, gene_tree, gene_root) for mk_ed in eval_metrics]
        # Counts MPRs to weight scores
        mpr_counter = ClusterUtil.mk_count_mprs(gene_root)
        # Evaluate the original graph on each eval metric to record improvement
        one_scores = [ClusterUtil.get_score_nodp([recon_g], eval_d, mpr_counter) for eval_d in eval_ds]
        # Perform the clustering for each cluster distance
        for i1, cluster_d in enumerate(cluster_ds):
            t = timeout_cluster(recon_g, gene_root, cluster_d, mpr_count, args, timeout)
            if t is None:
                print("{} timed out".format(f))
                continue
            graphs,_ = t
            # Evaluate the clustering for each evaluation distance
            for i2, eval_d in enumerate(eval_ds):
                one_score = one_scores[i2]
                k_score = ClusterUtil.get_score_nodp(graphs, eval_d, mpr_counter)
                improvement = 1 - (k_score / float(one_score))
                improvements[i1][i2].append(improvement)
    return improvements

# Make plottable xs and ys for each evaluation metric
# Specific to two dimensional output
def transform_d1_d2(improvements):
    series = []
    for i, evals in enumerate(improvements.values()):
        xs = evals[0]
        ys = evals[1]
        series.append((xs, ys))
    return series

def plot_d1_d2(improvements):
    series = transform_d1_d2(improvements)
    colors=["red", "blue"]
    for i,s in enumerate(series):
        xs, ys = s
        color = colors[i]
        plt.scatter(xs, ys, c=color, alpha=0.5)
    plt.ylim((0, 1))
    plt.xlim((0, 1))
    plt.xlabel("D1 improvement")
    plt.ylabel("D2 improvement")
    plt.title("Relative improvement")
    plt.savefig("d1_d2_plot.pdf", bbox_inches="tight")
    plt.clf()

def d1_d2_plot(trees, args):
    #cluster_names = ["PDV", "Event Support"]
    #cluster_metrics = [ClusterUtil.mk_pdv_dist, ClusterUtil.mk_support_dist]
    cluster_metrics = [ClusterUtil.mk_support_dist]
    eval_names = ["PDV", "Event Support"]
    eval_metrics = [ClusterUtil.mk_pdv_dist, ClusterUtil.mk_support_dist]
    d1_d2_improvements = get_improvements(trees, cluster_metrics, eval_metrics, args)
    # Now plot them!
    plot_d1_d2(d1_d2_improvements)

# # # #
# Correlate improvement with the number of splits used
# # # #

def get_n_improvements(tree_files, metric, args, timeout=1200, min_mprs=1000):
    series = []
    n = len(tree_files)
    for (i, f) in enumerate(tree_files):
        print("{}: {}/{}".format(f, i, n))
        # Get the recon graph + other info
        gene_tree, species_tree, gene_root, recon_g, mpr_count = \
            ClusterUtil.get_tree_info(str(f), args.d,args.t,args.l)
        #print("MPR count: {}".format(mpr_count))
        # Only care about trees with a certain number of MPRs
        if mpr_count < min_mprs:
            continue
        mpr_counter = ClusterUtil.mk_count_mprs(gene_root)
        d = metric(species_tree, gene_tree, gene_root)
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
                # Don't bother re-clustering if we split to get the same number of gs
                if len(gs) == old_ngs:
                    continue
                else:
                    old_ngs = len(gs)
                    _, scores = ClusterUtil.combine(gs, d, args.k, mpr_counter)
            except TimeoutError:
                print("{} timed out".format(f))
                continue
            signal.alarm(0)

            true_n = args.k + len(scores)
            new, old = scores[0]
            improvement = (1 - new/float(old))
            xs.append(true_n)
            ys.append(improvement)
        series.append((xs[:], ys[:]))
    return series

def plot_n_improvement(series):
    # Different color for each series
    n = len(series)
    norm = matplotlib.colors.Normalize(0, n)
    cmap = matplotlib.cm.get_cmap("hsv", n)
    for i, (xs, ys) in enumerate(series):
        plt.scatter(xs, ys, c=cmap(norm(i)), alpha=0.5)
    # Improvement is bounded between 0 and 1
    #plt.ylim((0,1))
    #plt.xlim((0,500))
    plt.xlabel("Number of Splits")
    plt.ylabel("Improvement")
    plt.title("Relative improvement")
    plt.savefig("n_improvement_plot.pdf", bbox_inches="tight")
    plt.clf()

def n_improvement_plot(trees, args):
    metric = ClusterUtil.mk_support_dist
    series = get_n_improvements(trees, metric, args)
    plot_n_improvement(series)

#MAIN
#TODO: just do all the required calculations ahead of time
# then choose which plot to create?

def main():
    args = process_args()
    trees = get_tree_files(args.input)
    if args.d1d2:
        d1_d2_plot(trees, args)
    if args.ki:
        k_improvement_plot(trees, args)
    if args.ni:
        n_improvement_plot(trees, args)

if __name__ == "__main__":
    main()

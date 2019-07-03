import ClusterUtil
import ReconciliationVisualization as RV
import HistogramDisplay

import argparse
from pathlib import Path

# Ideas
# - Correlate normality with improvement
# - Correlate improvement with depth / no. of graphs before clustering
# - Correlate improvement with the number of clusters
# - Just a histogram of the improvements for each tree
# - Normality distance (i.e. Normality statistic instead of average for the PDV)
#   - How does this correlate with improvements in the other metrics?
# - Generally, for each metric how does it improve the other metrics
#   - Ex. correlate PDV improvement with Support improvement for clusterings created
#       using the PDV (or support) distance.
# - Use "increases average improvement the most" rather than "smallest distance between pairs"
#   for deciding which to merge
#   - is this the same?

# Improvement vs. number of clusters, but improvement is vs. 1 cluster only

def process_args():
    # Required arguments - input file, D T L costs
    parser = argparse.ArgumentParser("")
    parser.add_argument("--input", metavar="<filename>", required=True,
        help="The path to a .newick file with the input trees and tip mapping.")
    parser.add_argument("-d", type=int, metavar="<duplication_cost>", required=True,
        help="The relative cost of a duplication.")
    parser.add_argument("-t", type=int, metavar="<transfer_cost>", required=True,
        help="The relative cost of a transfer.")
    parser.add_argument("-l", type=int, metavar="<loss_cost>", required=True,
        help="The relative cost of a loss.")
    parser.add_argument("-k", type=int, metavar="<number_of_clusters>", required=True,
        help="How many clusters to create.")
    parser.add_argument("--supportdist", action="store_true",
        help="Use the event support distance. Default is PDV distance.")
    depth_or_n = parser.add_mutually_exclusive_group(required=True)
    depth_or_n.add_argument("--depth", type=int, metavar="<tree_depth>",
        help="How far down the graph to consider event splits.")
    depth_or_n.add_argument("--nmprs", type=int, metavar="<tree_depth>",
        help="How many MPRs to consider")
    args = parser.parse_args()
    return args

def main():
    args = process_args()
    # Choose the distance metric
    if args.supportdist:
        mk_d = ClusterUtil.mk_support_dist
    else:
        mk_d = ClusterUtil.mk_pdv_dist
    # Get the recon graph + other info
    gene_tree, species_tree, gene_root, recon_g, mpr_count = \
        ClusterUtil.get_tree_info(args.input, args.d,args.t,args.l)

    # Visualize the graphs
    #RV.visualizeAndSave(recon_g, "original.png")
    #gs = ClusterUtil.full_split(recon_g, gene_root, args.depth)
    #for i, g in enumerate(gs):
    #    RV.visualizeAndSave(g, "{}.png".format(i))
    
    mpr_counter = ClusterUtil.mk_count_mprs(gene_root)
    # Make the distance metric for these specific trees
    d = mk_d(species_tree, gene_tree, gene_root)
    # Actually perform the clustering
    if args.depth is not None:
        graphs,scores = ClusterUtil.cluster_graph(recon_g, gene_root, d, args.depth, args.k)
    elif args.nmprs is not None:
        graphs,scores = ClusterUtil.cluster_graph_n(recon_g, gene_root, d, args.nmprs, mpr_count, args.k)
    else:
        assert False
    #TODO: what do with the graphs
    #TODO: visualize the graph for each of the new gs and the original
    # Visualization
    #get_hist = ClusterUtil.mk_get_hist(species_tree, gene_tree, gene_root)
    #cost_suffix = ".{}-{}-{}".format(args.d, args.t, args.l)
    #p = Path(args.input)
    #orig_p = str(p.with_suffix(cost_suffix + ".pdf"))
    #orig_h = get_hist(recon_g).histogram_dict
    #HistogramDisplay.plot_histogram(orig_p, orig_h, 1, Path(args.input).stem, args.d, args.t, args.l)
    #for i, g in enumerate(graphs):
    #    g_i = "cluster{}".format(i)
    #    g_p = str(p.with_suffix("." + g_i + cost_suffix + ".pdf"))
    #    g_h = get_hist(g).histogram_dict
    #    HistogramDisplay.plot_histogram(g_p, g_h, 1, Path(args.input).stem + g_i, args.d, args.t, args.l)
    # Statistics
    one_score = ClusterUtil.get_score_nodp([recon_g], d, mpr_counter)
    k_score = ClusterUtil.get_score_nodp(graphs, d, mpr_counter)
    improvement_ratio = k_score / float(one_score)
    improvement = 1 - improvement_ratio
    print("Old distance: {}".format(one_score))
    print("New distance: {}".format(k_score))
    print("Improvement:  {}".format(improvement))

def main2():
    args = process_args()
    gene_tree, species_tree, gene_root, recon_g, mpr_count = \
        ClusterUtil.get_tree_info(args.input, args.d,args.t,args.l)
    RV.visualizeAndSave(recon_g, "original.png")
    gs = ClusterUtil.full_split(recon_g, gene_root, args.depth)
    for i, g in enumerate(gs):
        RV.visualizeAndSave(g, "{}.png".format(i))

if __name__ == "__main__":
    main()

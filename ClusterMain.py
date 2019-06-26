import ClusterUtil
import ReconciliationVisualization as RV

import argparse


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
    parser.add_argument("--depth", type=int, metavar="<tree_depth>", required=True,
        help="How far down the graph to consider event splits.")
    parser.add_argument("-k", type=int, metavar="<number_of_clusters>", required=True,
        help="How many clusters to create.")
    parser.add_argument("--supportdist", action="store_true",
        help="Use the event support distance. Default is PDV distance.")
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
    gene_tree, species_tree, gene_tree_root, recon_g = \
        ClusterUtil.get_tree_info(args.input, args.d,args.t,args.l)

    # Visualize the graphs
    #RV.visualizeAndSave(recon_g, "original.png")
    #gs = ClusterUtil.full_split(recon_g, gene_tree_root, args.depth)
    #for i, g in enumerate(gs):
    #    RV.visualizeAndSave(g, "{}.png".format(i))

    # Make the distance metric for these specific trees
    d = mk_d(species_tree, gene_tree, gene_tree_root)
    # Actually perform the clustering
    graphs = ClusterUtil.cluster_graph(recon_g, gene_tree_root, d, args.depth, args.k)
    #TODO: what do with the graphs
    old_d =  d(recon_g, recon_g)
    new_ds = [d(g,g) for g in graphs]
    new_d = sum(new_ds) / float(len(new_ds))
    improvement = new_d / old_d
    print(new_ds)
    print("Old distance: {}".format(old_d))
    print("New distance: {}".format(new_d))
    print("Improvement:  {}".format(improvement))

def main2():
    args = process_args()
    gene_tree, species_tree, gene_tree_root, recon_g = \
        ClusterUtil.get_tree_info(args.input, args.d,args.t,args.l)
    RV.visualizeAndSave(recon_g, "original.png")
    gs = ClusterUtil.full_split(recon_g, gene_tree_root, args.depth)
    for i, g in enumerate(gs):
        RV.visualizeAndSave(g, "{}.png".format(i))

if __name__ == "__main__":
    main()

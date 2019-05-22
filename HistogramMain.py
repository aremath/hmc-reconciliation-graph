
import DTLReconGraph
import HistogramAlg
import Diameter
import HistogramDisplay


def main():
    tree_file = "newickSample/size7/test-size7-no41.newick"
    D = 1
    T = 1
    L = 1

    # From the newick tree create the reconciliation graph
    edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots \
        = DTLReconGraph.reconcile(tree_file, D, T, L)

    # Reformat the host and parasite tree to use it with the histogram algorithm
    gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")
    species_tree, species_tree_root, species_node_count \
        = Diameter.reformat_tree(edge_species_tree, "hTop")

    # Calculate the histogram via histogram algorithm
    diameter_alg_hist = HistogramAlg.diameter_algorithm(
        species_tree, gene_tree, gene_tree_root, dtl_recon_graph, dtl_recon_graph,
        False, False, verify=True)

    HistogramDisplay.plot_histogram("histnonorm.pdf", diameter_alg_hist.histogram_dict, False, False)
    HistogramDisplay.plot_histogram("histxnorm.pdf", diameter_alg_hist.histogram_dict, True, False)
    HistogramDisplay.plot_histogram("histynorm.pdf", diameter_alg_hist.histogram_dict, False, True)
    HistogramDisplay.plot_histogram("histxynorm.pdf", diameter_alg_hist.histogram_dict, True, True)

if __name__ == "__main__":
    main()


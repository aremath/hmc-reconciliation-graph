import HistogramAlg
import DTLMedian
import DTLReconGraph
import Diameter

import itertools
import functools
from collections import deque
import numpy as np

def graph_union(g1, g2):
    """ Return the union of two recon graphs.
    """
    newg = {}
    for k,v in g1.iteritems():
        newg[k] = v[:]
    for k,v in g2.iteritems():
        if k in newg:
            for e in v:
                if e not in newg[k]:
                    newg[k].append(e)
        else:
            newg[k] = v[:]
    return newg

def graph_is_subset(g1, g2):
    """ True iff g1 is a subset of g2
        (contains a subset of mapping / event nodes)
    """
    for k,v in g1.iteritems():
        if k not in g2:
            return False
        for e in v:
            if e not in g2[k]:
                return False
    return True

def assert_pairwise_nonsub(gs):
    """ Assert that the graphs in gs are pairwise not subsets of each other.
    """
    for i1, g1 in enumerate(gs):
        for i2, g2 in enumerate(gs):
            if i2 > i1:
                gis = graph_is_subset(g1, g2) or graph_is_subset(g2, g1)
                assert not gis

def graph_sub(g, root):
    """ Returns a subset of g which is all the nodes reachable from root.
    """
    newg = {}
    assert root in g
    finished = set([root])
    q = deque([root])
    while len(q) > 0:
        p = q.popleft()
        newg[p] = g[p]
        # Add its mapping node children to the queue
        for e, m1, m2 in g[p]:
            # Both children of a tip mapping event are Nones
            if m1 not in finished and m1 != (None, None):
                q.append(m1)
                finished.add(m1)
            # Second child of a loss is (None, None)
            if m2 not in finished and m2 != (None, None):
                q.append(m2)
                finished.add(m2)
    return newg

# Split at higher depth until we have at least n MPRs
# Use full_split_n instead since there may be multiple roots.
#NOTE: UNUSED
def graph_split_n(g, root, n, mpr_count):
    depth=1
    while True:
        gs = graph_split(g, root, depth)
        if len(gs) >= n:
            break
        # We can't get more graphs than there are MPRs
        elif len(gs) == mpr_count:
            break
        depth += 1
    return gs

def graph_split(g, root, depth):
    """ Find a set of sub-graphs by starting at the root and making every
        possible choice of event up to the given depth.
    """
    # Base Case: 0 depth means just descend from the root choosing all events.
    if depth == 0:
        return [graph_sub(g, root)]
    else:
        gs = []
        # Make one choice of event
        for e, m1, m2 in g[root]:
            # Recursively enumerate the possible event choices below that
            # In case e is a tip-mapping
            if m1 == (None, None):
                l = [{}]
            else:
                l = graph_split(g, m1, depth-1)
            # In case e is a loss
            if m2 == (None, None):
                r = [{}]
            else:
                r = graph_split(g, m2, depth-1)
            # Merge each pair of child graphs through the event
            for lg, rg in itertools.product(l, r):
                newg = graph_union(lg, rg)
                newg[root] = [(e, m1, m2)]
                gs.append(newg)
        return gs

# Split at higher depth until n mprs are found
def full_split_n(g, gene_root, n, mpr_count):
    """ full_split, but increases the depth until n splits are found.
        mpr_count is the total number of mprs in g
    """
    depth=1
    while True:
        gs = full_split(g, gene_root, depth)
        if len(gs) >= n:
            break
        # We can't get more graphs than there are MPRs
        elif len(gs) == mpr_count:
            break
        depth += 1
    print("Depth: {}".format(depth))
    return gs

def full_split(g, gene_root, depth):
    """ Use split_gs to split every root to find all of the splits at depth d
        for the entire recon graph.
    """
    # Find the mapping nodes involving the gene root
    roots = [k for k in g.keys() if k[0] == gene_root]
    #TODO: if the top node is lost, then that loss will not be a root
    split_gs = [graph_split(g, r, depth) for r in roots]
    # Flatten
    f = functools.reduce(lambda a, b: a+b, split_gs, [])
    return f

def new_dp(dp, m1, m2):
    """ Helper for combine that updates the DP table for a re-indexed list
    """
    # The new index corresponding to i if the elements at m1 and m2 are removed
    def new_index(i):
        if i < m1:
            return i
        elif i < m2:
            return i - 1
        elif i > m2:
            return i - 2
        assert False
    newdp = {}
    for i, dist in dp.iteritems():
        # update the index for each index that does not involve m1 or m2
        if i not in [m1, m2]:
            newdp[new_index(i)] = dist
    return newdp
    # The new index corresponding to i if the elements at m1 and m2 are removed
    def new_index(i):
        if i < m1:
            return i
        elif i < m2:
            return i - 1
        elif i > m2:
            return i - 2
        assert False
    newdp = {}
    for i, dist in dp.iteritems():
        # update the index for each index that does not involve m1 or m2
        if i not in [m1, m2]:
            newdp[new_index(i)] = dist
    return newdp

def new_u_dp(u_dp, m1, m2):
    """ Re-index a union dp (where the keys are tuples)
    """
    # The new index corresponding to i if the elements at m1 and m2 are removed
    def new_index(i):
        if i < m1:
            return i
        elif i < m2:
            return i - 1
        elif i > m2:
            return i - 2
        assert False
    newdp = {}
    for pair, dist in u_dp.iteritems():
        i1, i2 = pair
        # update the index for each index pair that does not involve m1 or m2
        if i1 not in [m1, m2] and i2 not in [m1, m2]:
            newdp[(new_index(i1), new_index(i2))] = dist
    return newdp

def get_score_vals(graphs, g_score, score_dp, nmprs_dp, mpr_counter):
    """ Compute the relevant score values for a list of graphs
    """
    weighted_score = 0
    total_nmprs = 0
    for i, g in enumerate(graphs):
        if i in score_dp and i in nmprs_dp:
            score = score_dp[i]
            nmprs = nmprs_dp[i]
        else:
            score = g_score(g)
            nmprs = mpr_counter(g)
            # Update DP tables
            score_dp[i] = score
            nmprs_dp[i] = nmprs
        weighted_score += score * nmprs
        total_nmprs += nmprs
    return weighted_score, total_nmprs

def get_score_vals_ignore(graphs, g_score, score_dp, nmprs_dp, mpr_counter, i1, i2):
    """ get_score_vals, but ignore the graphs at indices i1 and i2.
        Cannot just pass graphs without those graphs because the DP
        relies on correct indexing.
    """
    weighted_score = 0
    total_nmprs = 0
    for i, g in enumerate(graphs):
        if i == i1 or i == i2:
            continue
        if i in score_dp and i in nmprs_dp:
            score = score_dp[i]
            nmprs = nmprs_dp[i]
        else:
            score = g_score(g)
            nmprs = mpr_counter(g)
            # Update DP tables
            score_dp[i] = score
            nmprs_dp[i] = nmprs
        weighted_score += score * nmprs
        total_nmprs += nmprs
    return weighted_score, total_nmprs

def get_score(graphs, g_score, score_dp, nmprs_dp, mpr_counter):
    """Compute the weighted average score.
    """
    weighted_score, total_nmprs = get_score_vals(graphs, g_score, score_dp, nmprs_dp, mpr_counter)
    return weighted_score / float(total_nmprs)

def get_score_nodp(graphs, g_score, mpr_counter):
    """ Compute the WAS without using a DP table.
        Used to verify get_score and evaluate a cluster without access to the DP tables.
    """
    weighted_score = 0
    total_nmprs = 0
    for g in graphs:
        score = g_score(g)
        nmprs = mpr_counter(g)
        weighted_score += score * nmprs
        total_nmprs += nmprs
    return weighted_score / float(total_nmprs)

def get_score_merge_indices(graphs, g_score, score_dp, score_u_dp, nmprs_dp, nmprs_u_dp, mpr_counter, i1, i2):
    """ Calculate the WAS for graphs but assuming that i1 and i2 are merged.
        Used to calculate which pair is the best to merge in combine.
    """
    # Calculate the weighted score without the pair to be merged
    weighted_score, total_nmprs = get_score_vals_ignore(graphs, g_score, score_dp, nmprs_dp, \
                                                        mpr_counter, i1, i2)
    # Add in the score for the pair
    pair = (i1, i2)
    if pair in score_u_dp and pair in nmprs_u_dp:
        u_score = score_u_dp[pair]
        u_nmprs = nmprs_u_dp[pair]
    else:
        gu = graph_union(graphs[i1], graphs[i2])
        u_score = g_score(gu)
        u_nmprs = mpr_counter(gu)
        # Store the new values in the tables
        score_u_dp[pair] = u_score
        nmprs_u_dp[pair] = u_nmprs
    weighted_score += u_score * u_nmprs
    total_nmprs += u_nmprs
    return weighted_score / float(total_nmprs)

def combine(split_gs, g_score, k, mpr_counter):
    """ Take a set of graph splits and cluster them by minimizing the WAS at every iteration.
    """
    assert k >= 1
    scores = []
    # Keep track of the score for each graph
    # Key - index into split_gs, value - score for that graph
    score_dp = {}
    # Key - index into split_gs, value - number of MPRS for that graph
    nmprs_dp = {}
    # Key - (i1, i2) indices into split_gs, value - score for g1 u g2
    score_u_dp = {}
    # Key - (i1, i2), value - number of mprs in g1 u g2
    nmprs_u_dp = {}
    # Compute the "distance" between two indices in the graphs table
    # Smaller distance is better since smaller score is better
    def distance(i1, i2):
        return get_score_merge_indices(split_gs, g_score,
            score_dp, score_u_dp, nmprs_dp, nmprs_u_dp, mpr_counter, i1, i2)
    # Merge until there are k (or fewer) splits
    while len(split_gs) > k:
        # Find a pair to merge via smallest distance
        # Start with 0, 1 (guaranteed to exist since k >= 1
        min_dist = distance(0, 1)
        min_i1 = 0
        min_i2 = 1
        for i1, g1 in enumerate(split_gs):
            for i2, g2 in enumerate(split_gs):
                # All combinations like itertools.combinations, but need the indices as well
                if i2 > i1:
                    dist = distance(i1, i2)
                    if dist < min_dist:
                        min_i1 = i1
                        min_i2 = i2
                        min_dist = dist
        # Compute the "old" score (before the merge)
        score = get_score(split_gs, g_score, score_dp, nmprs_dp, mpr_counter)
        scores.append(score)

        # Now merge them
        # Pop in reverse order (i2 > i1 necessarily) to not throw off the previous index
        gm1 = split_gs.pop(min_i2)
        gm2 = split_gs.pop(min_i1)
        gu = graph_union(gm1, gm2)
        split_gs.append(gu)
        
        #Debug:
        #dm1 = score_dp[(min_i1, min_i1)]
        #dm2 = score_dp[(min_i2, min_i2)]
        #du = min_dist * mpr_counter(gu)
        #print("Merge {} and {} to get {}".format(dm1, dm2, du))
        #print("MPR counts {}, {} -> {}".format(mpr_counter(gm1), mpr_counter(gm2), mpr_counter(gu)))
        #print("Is subset: {}".format(graph_is_subset(gm1, gm2) or graph_is_subset(gm2, gm1)))
        #print("Graph sizes {}, {} -> {}".format(len(gm1), len(gm2), len(gu)))

        # Fix up the DP tables now that the list is re-indexed
        # Entries involving gm1 and gm2 are removed and those distances will be
        # re-computed in the next while loop.
        score_dp = new_dp(score_dp, min_i1, min_i2)
        nmprs_dp = new_dp(nmprs_dp, min_i1, min_i2)
        score_u_dp = new_u_dp(score_u_dp, min_i1, min_i2)
        nmprs_u_dp = new_u_dp(nmprs_u_dp, min_i1, min_i2)

    final_score = get_score(split_gs, g_score, score_dp, nmprs_dp, mpr_counter)
    scores.append(final_score)
    # Reverse scores so that the first index of scores is the score for k clusters,
    # the second is for k+1 clusters, etc.
    return split_gs, scores[::-1]

def cluster_graph_n(graph, gene_root, distance, n, mpr_count, k, max_splits):
    """ Find k clusters within MPRs of g by first finding at least n splits
        then merging them by WAS.
    """
    # First split the graph
    gs = full_split_n(graph, gene_root, n, mpr_count)
    print("Number of splits: {}".format(len(gs)))
    if len(gs) > max_splits:
        return None
    mpr_counter = mk_count_mprs(gene_root)
    # Then recombine those splits until we have k graphs
    return combine(gs, distance, k, mpr_counter)

def cluster_graph(graph, gene_root, distance, depth, k, max_splits):
    """ Find k clusters within MPRs of g using a depth-splitting method
        then merging by WAS.
    """
    # First split the graph
    gs = full_split(graph, gene_root, depth)
    print("Number of splits: {}".format(len(gs)))
    if len(gs) > max_splits:
        return None
    mpr_counter = mk_count_mprs(gene_root)
    # Then recombine those splits until we have k graphs
    return combine(gs, distance, k, mpr_counter)

#TODO: this can be improved by keeping the partial DP table around.
# Since the graphs are always the same below a certain level, preserving the table below that level
# would mean less repeated computation.

def mk_pdv_score(species_tree, gene_tree, gene_root):
    """ Makes a score function for a graph by specifying the trees and the root.
        The score is the average pairwise distance.
    """
    def score(g):
        hist = HistogramAlg.diameter_algorithm(species_tree, gene_tree, gene_root, g, g, False, False)
        return hist.mean()
    return score

def avg_event_support(species_tree, gene_tree, g, gene_root):
    """ Compute the average event support for a graph.
    """
    # Compute the event support for each event
    preorder_mapping_nodes = DTLMedian.mapping_node_sort(gene_tree, species_tree, g.keys())
    event_support, count = \
        DTLMedian.generate_scores(list(reversed(preorder_mapping_nodes)), g, gene_root)
    # Take the average over each event
    total_support = 0
    for support in event_support.itervalues():
        total_support += support
    return total_support / len(event_support)

def mk_support_score(species_tree, gene_tree, gene_root):
    """ Make a score function by specifying the trees and the root.
        The score is the negative average event support.
        It is negated because the clustering algorithm minimizes the score,
        but we want to maximize the event support of each cluster.
    """
    def score(g):
        support = avg_event_support(species_tree, gene_tree, g, gene_root)
        # Higher support means closer, so take the reciprocal.
        #return 1.0 / support
        return -1 * support
    return score

def mk_get_pdv_hist(species_tree, gene_tree, gene_root):
    """ Partially apply diameter_algorithm on non-changing arguments
        for convenient use with multiple graphs.
    """
    def get_hist(g):
        h = HistogramAlg.diameter_algorithm(species_tree, gene_tree, gene_root, g, g, False, False)
        return h.histogram_dict
    return get_hist

def event_support_hist(species_tree, gene_tree, gene_root, graph):
    preorder_mapping_nodes = DTLMedian.mapping_node_sort(gene_tree, species_tree, graph.keys())
    event_support, count = \
        DTLMedian.generate_scores(list(reversed(preorder_mapping_nodes)), graph, gene_root)
    supports = event_support.values()
    hist, bins = np.histogram(supports, bins=20, range=(0,1))
    total = np.sum(hist)
    return hist / float(total), bins

def mk_get_support_hist(species_tree, gene_tree, gene_root):
    """ Same as above but for the event support histogram. """
    def get_hist(g):
        return event_support_hist(species_tree, gene_tree, gene_root, g)
    return get_hist

# Create a function that counts the number of mprs
def mk_count_mprs(gene_root):
    """ Partially apply the MPR-counting function on non-changing arguments
        for convenient use with multiple graphs.
    """
    def count_mprs(g):
        # Find the mapping nodes involving the gene root
        roots = [k for k in g.keys() if k[0] == gene_root]
        return DTLReconGraph.count_mprs_wrapper(roots, g)
    return count_mprs

def calc_improvement(big_k, little_k):
    # For Event Support
    return big_k / float(little_k)
    # For PDV
    #return little_k / float(big_k)

# For Event Support
def calc_improvement_support(big_k, little_k):
    return big_k / float(little_k)

# For PDV
def calc_improvement_pdv(big_k, little_k):
    return little_k / float(big_k)

def get_tree_info(newick, d,t,l):
    """ Reconcile the trees and return all the relevant info.
    """
    # From the newick tree create the reconciliation graph
    edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots \
        = DTLReconGraph.reconcile(newick, d, t, l)
    # Reformat the host and parasite tree to use it with the histogram algorithm
    gene_tree, gene_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")
    species_tree, species_tree_root, species_node_count \
        = Diameter.reformat_tree(edge_species_tree, "hTop")
    return gene_tree, species_tree, gene_root, dtl_recon_graph, mpr_count, best_roots


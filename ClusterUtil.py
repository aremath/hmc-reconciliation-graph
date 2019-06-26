import HistogramAlg
import DTLMedian
import DTLReconGraph
import Diameter

import itertools
import functools
from collections import deque

# Union of two reconciliation graphs
def graph_union(g1, g2):
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

# Subset of g of nodes reachable from root
def graph_sub(g, root):
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

# Find a set of sub-graphs by starting at root, and making every possible event
# choice up to a given depth.
def graph_split(g, root, depth):
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

# Split every root, then flatten the list
def full_split(g, gene_root, depth):
    # Find the mapping nodes involving the gene root
    roots = [k for k in g.keys() if k[0] == gene_root]
    #TODO: if the top node is lost, then that loss will not be a root
    split_gs = [graph_split(g, r, depth) for r in roots]
    # Flatten
    f = functools.reduce(lambda a, b: a+b, split_gs, [])
    return f

# Helper for combine that updates the DP table for a re-indexed list
def new_dp(dp, m1, m2):
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
    for pair, dist in dp.iteritems():
        i1, i2 = pair
        # Update the index for each index pair that does not involve m1 or m2
        if i1 not in [m1, m2] and i2 not in [m1, m2]:
            newdp[(new_index(i1), new_index(i2))] = dist
    return newdp

def combine(split_gs, distance, k):
    assert k > 1
    # Use a DP table to keep track of the already-computed distances
    dp = {}
    while len(split_gs) > k:
        # Find a pair to merge via smallest distance
        # Start with 0, 1 (guaranteed to exist since k > 1
        if (0, 1) in dp:
            min_dist = dp[(0, 1)]
        else:
            min_dist = distance(split_gs[0], split_gs[1])
            dp[(0, 1)] = min_dist
        min_i1 = 0
        min_i2 = 1
        for i1, g1 in enumerate(split_gs):
            for i2, g2 in enumerate(split_gs):
                # All combinations like itertools.combinations, but need the indices as well
                if i2 > i1:
                    # Get the distance from the DP table if possible
                    if (i1, i2) in dp:
                        d = dp[(i1, i2)]
                    else:
                        d = distance(g1, g2)
                        dp[(i1, i2)] = d
                    if d < min_dist:
                        min_i1 = i1
                        min_i2 = i2
                        min_dist = d
        # Now merge them
        # Pop in reverse order (i2 > i1 necessarily) to not throw off the previous index
        gm1 = split_gs.pop(min_i2)
        gm2 = split_gs.pop(min_i1)
        gu = graph_union(gm1, gm2)
        split_gs.append(gu)
        # Fix up the DP table now that the list is re-indexed
        # Entries involving gm1 and gm2 are removed and those distances will be
        # re-computed in the next while loop.
        dp = new_dp(dp, min_i1, min_i2)
    return split_gs

# Find k clusters within MPRs of g using a depth-splitting method
def cluster_graph(graph, gene_root, distance, depth, k):
    # First split the graph
    gs = full_split(graph, gene_root, depth)
    print("Number of splits: {}".format(len(gs)))
    # Then recombine those splits until we have k graphs
    return combine(gs, distance, k)

#TODO: this can be improved by keeping the partial DP table around.
# Since the graphs are always the same below a certain level, preserving the table below that level
# would mean less repeated computation.

# NOTE: Neither of these are really true metrics since the distance between g and itself is nonzero.
# Makes a distance on graphs by specifying the trees and the root
# The distance is the average pairwise distance between MPRs in the combined graph
def mk_pdv_dist(species_tree, gene_tree, gene_root):
    def d(g1, g2):
        #hist = HistogramAlg.diameter_algorithm(species_tree, gene_tree, gene_root, g1, g2, False, False)
        gu = graph_union(g1, g2)
        hist = HistogramAlg.diameter_algorithm(species_tree, gene_tree, gene_root, gu, gu, False, False)
        return hist.mean()
    return d

def avg_event_support(species_tree, gene_tree, g, gene_root):
    # Compute the event support for each event
    preorder_mapping_nodes = DTLMedian.mapping_node_sort(gene_tree, species_tree, g.keys())
    event_support, count = \
        DTLMedian.generate_scores(list(reversed(preorder_mapping_nodes)), g, gene_root)
    # Take the average
    total_support = 0
    for support in event_support.itervalues():
        total_support += support
    return total_support / len(event_support)

# The distance here is the reciprocal of the average event support
def mk_support_dist(species_tree, gene_tree, gene_root):
    def d(g1, g2):
        gu = graph_union(g1, g2)
        support_u = avg_event_support(species_tree, gene_tree, gu, gene_root)
        #support_1 = avg_event_support(species_tree, gene_tree, g1, gene_root)
        #support_2 = avg_event_support(species_tree, gene_tree, g2, gene_root)
        #avg_individual = (support_1 + support_2) / 2.0

        # Higher support means closer, so take the reciprocal.
        return 1.0 / support_u
    return d

# Partially apply diameter_algorithm on the non-chaning elements
# This makes use with multiple graphs more convenient
def mk_get_hist(species_tree, gene_tree, gene_root):
    def get_hist(g):
        return HistogramAlg.diameter_algorithm(species_tree, gene_tree, gene_root, g, g, False, False)
    return get_hist

def get_tree_info(newick, d,t,l):
    # From the newick tree create the reconciliation graph
    edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots \
        = DTLReconGraph.reconcile(newick, d, t, l)
    # Reformat the host and parasite tree to use it with the histogram algorithm
    gene_tree, gene_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")
    species_tree, species_tree_root, species_node_count \
        = Diameter.reformat_tree(edge_species_tree, "hTop")
    return gene_tree, species_tree, gene_root, dtl_recon_graph, mpr_count


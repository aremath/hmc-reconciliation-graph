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

def graph_is_subset(g1, g2):
    for k,v in g1.iteritems():
        if k not in g2:
            return False
        for e in v:
            if e not in g2[k]:
                return False
    return True

def assert_pairwise_nonsub(gs):
    for i1, g1 in enumerate(gs):
        for i2, g2 in enumerate(gs):
            if i2 > i1:
                gis = graph_is_subset(g1, g2) or graph_is_subset(g2, g1)
                assert not gis
    assert False

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

# Split at higher depth until we have at least n MPRs
# Use full_split_n instead since there may be multiple roots.
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

# Split at higher depth until n mprs are found
def full_split_n(g, gene_root, n, mpr_count):
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
    for i, dist in dp.iteritems():
        # update the index for each index that does not involve m1 or m2
        if i not in [m1, m2]:
            newdp[new_index(i)] = dist
    return newdp

# Re-index a union dp (where the keys are tuples)
def new_u_dp(u_dp, m1, m2):
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

# Actually compute the relevant score values
def get_score_vals(graphs, g_score, score_dp, nmprs_dp, mpr_counter):
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

# Score is weighted by the number of MPRs
def get_score(graphs, g_score, score_dp, nmprs_dp, mpr_counter):
    weighted_score, total_nmprs = get_score_vals(graphs, g_score, score_dp, nmprs_dp, mpr_counter)
    return weighted_score / float(total_nmprs)

# Used to verify get_score and evaluated a cluster using a different metric
def get_score_nodp(graphs, g_score, mpr_counter):
    weighted_score = 0
    total_nmprs = 0
    for g in graphs:
        score = g_score(g)
        nmprs = mpr_counter(g)
        weighted_score += score * nmprs
        total_nmprs += nmprs
    return weighted_score / float(total_nmprs)

# Score as if the graphs at i1 and i2 are merged
def get_score_merge_indices(graphs, g_score, score_dp, score_u_dp, nmprs_dp, nmprs_u_dp, mpr_counter, i1, i2):
    # New list without the pair to be merged
    newgs = [g for i,g in enumerate(graphs) if i != i1 and i != i2]
    weighted_score, total_nmprs = get_score_vals(newgs, g_score, score_dp, nmprs_dp, mpr_counter)
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

#TODO: decide if we want "improvement over last time"
# or "improvement over one cluster", or something else
def combine(split_gs, g_score, k, mpr_counter):
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
        old_score = get_score(split_gs, g_score, score_dp, nmprs_dp, mpr_counter)

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

        # Compute the "new" score (after the merge)
        new_score = get_score(split_gs, g_score, score_dp, nmprs_dp, mpr_counter)

        # The improvement can actually be negative if we merge a split which is a subset of another split.
        #imp = 1-(old_score/float(new_score))
        #assert imp > 0, imp

        scores.append((old_score, new_score))
    # Reverse scores so that the first index of scores is the score for k clusters,
    # the second is for k+1 clusters, etc.
    return split_gs, scores[::-1]

# Find k clusters within MPRs of g using a depth-splitting method
def cluster_graph_n(graph, gene_root, distance, n, mpr_count, k):
    # First split the graph
    gs = full_split_n(graph, gene_root, n, mpr_count)
    print("Number of splits: {}".format(len(gs)))
    mpr_counter = mk_count_mprs(gene_root)
    # Then recombine those splits until we have k graphs
    return combine(gs, distance, k, mpr_counter)

# Find k clusters within MPRs of g using a depth-splitting method
def cluster_graph(graph, gene_root, distance, depth, k):
    # First split the graph
    gs = full_split(graph, gene_root, depth)
    print("Number of splits: {}".format(len(gs)))
    mpr_counter = mk_count_mprs(gene_root)
    # Then recombine those splits until we have k graphs
    return combine(gs, distance, k, mpr_counter)

#TODO: this can be improved by keeping the partial DP table around.
# Since the graphs are always the same below a certain level, preserving the table below that level
# would mean less repeated computation.

# NOTE: Neither of these are really true metrics since the distance between g and itself is nonzero.
# Makes a distance on graphs by specifying the trees and the root
# The distance is the average pairwise distance between MPRs in the combined graph
def mk_pdv_score(species_tree, gene_tree, gene_root):
    def score(g):
        hist = HistogramAlg.diameter_algorithm(species_tree, gene_tree, gene_root, g, g, False, False)
        return hist.mean()
    return score

def avg_event_support(species_tree, gene_tree, g, gene_root):
    # Compute the event support for each event
    preorder_mapping_nodes = DTLMedian.mapping_node_sort(gene_tree, species_tree, g.keys())
    event_support, count = \
        DTLMedian.generate_scores(list(reversed(preorder_mapping_nodes)), g, gene_root)
    # Take the average over each event
    total_support = 0
    for support in event_support.itervalues():
        total_support += support
    return total_support / len(event_support)

# The score here is inversely related to the average event support
def mk_support_score(species_tree, gene_tree, gene_root):
    def score(g):
        support = avg_event_support(species_tree, gene_tree, g, gene_root)
        # Higher support means closer, so take the reciprocal.
        #return 1.0 / support
        return -1 * support
    return score

# Partially apply diameter_algorithm on the non-chaning elements
# This makes use with multiple graphs more convenient
def mk_get_hist(species_tree, gene_tree, gene_root):
    def get_hist(g):
        return HistogramAlg.diameter_algorithm(species_tree, gene_tree, gene_root, g, g, False, False)
    return get_hist

# Create a function that counts the number of mprs
def mk_count_mprs(gene_root):
    def count_mprs(g):
        # Find the mapping nodes involving the gene root
        roots = [k for k in g.keys() if k[0] == gene_root]
        return DTLReconGraph.count_mprs_wrapper(roots, g)
    return count_mprs

def calc_improvement(big_k, little_k):
    return big_k / float(little_k)

def get_tree_info(newick, d,t,l):
    # From the newick tree create the reconciliation graph
    edge_species_tree, edge_gene_tree, dtl_recon_graph, mpr_count, best_roots \
        = DTLReconGraph.reconcile(newick, d, t, l)
    # Reformat the host and parasite tree to use it with the histogram algorithm
    gene_tree, gene_root, gene_node_count = Diameter.reformat_tree(edge_gene_tree, "pTop")
    species_tree, species_tree_root, species_node_count \
        = Diameter.reformat_tree(edge_species_tree, "hTop")
    return gene_tree, species_tree, gene_root, dtl_recon_graph, mpr_count


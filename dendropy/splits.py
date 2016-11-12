#! /usr/bin/env python

############################################################################
##  splits.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Split calculation and management.
"""

from dendropy import taxa
from dendropy import trees
from dendropy import treegen
from dendropy import utils
from dendropy import get_logger
_LOG = get_logger("dendropy.splits")

def lowest_bit_only(s):
    m = s & (s - 1)
    return m ^ s

__n_bits_set = (0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4)
def count_bits(a):
    '''Returns the number of bits set to one.'''
    global __n_bits_set
    c = long(a)
    if c != a:
        raise ValueError('non-integer argument')
    #if c < 1L:
    if c < 1:
        #if c < 0L:
        if c < 0:
            raise ValueError('negative argument')
        return 0
    n_bits = 0
    while c > 0:
        #i = c & 0x0FL
        i = c & 0x0F
        n_bits += __n_bits_set[i]
        c >>= 4
    return n_bits

def split_to_list(s, mask=-1, one_based=False, ordination_in_mask=False):
    return [i for i in iter_split_indices(s, mask, one_based, ordination_in_mask)]
    
def iter_split_indices(s, mask=-1, one_based=False, ordination_in_mask=False):
    '''returns the index of each bit that is on in `s` and the `mask`
    
        Iy 'one_based` is True then the 0x01 bit is returned as 1 instead of 0.
        If `ordination_in_mask` is True then the indices returned will be the 
            count of the 1's in the mask that are to the right of the bit rather
            than the total number of digits to the right of the bit. Thus, the 
            index will be the index in a taxon block that is the subset of the 
            full set of taxa).
    '''
    currBitIndex = one_based and 1 or 0
    #test_bit = 1L
    test_bit = 1
    maskedSplitRep = s & mask
    standard_ordination = not ordination_in_mask
    while test_bit <= maskedSplitRep:
        if maskedSplitRep & test_bit:
            yield currBitIndex
        if standard_ordination or (mask & test_bit):
            currBitIndex += 1
        test_bit <<=1

def is_trivial_split(split, mask):
    """Returns True if the split occurs in any tree of the taxa `mask` -- if 
    there is only fewer than two 1's or fewer than two 0's in `split` (among 
    all of the that are 1 in mask)."""
    masked = split & mask
    if split == 0 or split == mask:
        return True
    if ((masked - 1) & masked) == 0:
        return True
    cm = (~split) & mask
    if ((cm - 1) & cm) == 0:
        return True
    return False
    
def is_non_singleton_split(split, mask):
    "Returns True if a split is NOT between a leaf and the rest of the taxa."
    # ((split-1) & split) is True (non-zero) only
    # if split is not a power of 2, i.e., if split
    # has more than one bit turned on, i.e., if it
    # is a non-trivial split.    
    return not is_trivial_split(split, mask)

    
def split_as_string(split_mask, width, symbol1=None, symbol2=None):
    "Returns a 'pretty' split representation."
    s = taxa.int_to_bitstring(split_mask).rjust(width, '0')
    if symbol1 is not None:
        s = s.replace('0', symbol1)
    if symbol2 is not None:
        s = s.replace('1', symbol2)
    return s
    
def split_as_string_rev(split_mask, width, symbol1='.', symbol2='*'):
    """
    Returns a 'pretty' split representation, reversed, with first taxon
    on the left (as given by PAUP*)
    """
    return split_as_string(split_mask=split_mask, 
                           width=width,
                           symbol1=symbol1,
                           symbol2=symbol2)[::-1]
    
def split_taxa_list(split_mask, taxa_block, index=0):
    "Returns list of taxa represented by split."
    taxa = []
    while split_mask:
        if split_mask & 1:
            taxa.append(taxa_block[index])
        split_mask = split_mask >> 1
        index += 1
    return taxa

def encode_splits(tree, create_dict=True, delete_degree_two=True):
    """
    Processes splits on a tree, encoding them as bitmask on each edge.    
    Adds the following to each edge:
        - `clade_mask` : a rooted split representation, i.e. a long/bitmask 
            where bits corresponding to indexes of taxa descended from this 
            edge are turned on
    If `create_dict` is True, then the following is added to the tree:
        - `split_edges`:
            [if `tree.is_rooted`]: a dictionary where keys are the
            splits and values are edges.        
            [othersie]: a NormalizedBitmaskDictionary where the keys are the
            normalized (unrooted) split representations and the values
            are edges. A normalized split_mask is where the clade_mask
            is complemented if the right-most bit is not '1' (or just
            the clade_mask otherwise).
    If `delete_degree_two` is True then degree 2 nodes (outdegree one in the
        rooted case) will be deleted as they are encountered (this is required
        if the split_edges dictionary is to refer to all edges in the tree).
    """
    taxa_block = tree.taxa_block
    if taxa_block is None:
        taxa_block = tree.infer_taxa_block()

    if create_dict:
        if tree.is_rooted:
            tree.split_edges = {}
        else:
            atb = taxa_block.all_taxa_bitmask()
            d = utils.NormalizedBitmaskDict(mask=atb)
            tree.split_edges = d
        split_map = tree.split_edges
    if delete_degree_two and (not tree.is_rooted):
        c = tree.seed_node.child_nodes()
        if len(c) == 2:
            tree.deroot()

    for edge in tree.postorder_edge_iter():
        cm = 0
        h = edge.head_node
        child_nodes = h.child_nodes()
        nc = len(child_nodes)
        if nc > 0:
            if nc == 1 and delete_degree_two:
                p = edge.tail_node
                c = child_nodes[0]
                try:
                    c.edge.length += edge.length
                except:
                    pass
                if p is not None:
                    pos = p.child_nodes().index(h)
                    p.add_child(c, pos=pos)
                    p.remove_child(h)
                else:
                    assert h is tree.seed_node
            else:
                for child in child_nodes:
                    cm |= child.edge.clade_mask
        else:
            t = edge.head_node.taxon
            if t:
                cm = taxa_block.taxon_bitmask(t)
        edge.clade_mask = cm
        if create_dict:
            split_map[cm] = edge
                
############################################################################        
## SplitDistribution

class SplitDistribution(object):
    "Collects information regarding splits over multiple trees."
    
    def __init__(self, taxa_block=None):
        "What else?"
        self.total_trees_counted = 0
        self.taxa_block = taxa_block
        self.splits = []
        self.split_counts = {}
        self.split_edge_lengths = {}
        self.split_node_ages = {}
        self.ignore_edge_lengths = False
        self.ignore_node_ages = True
        self.unrooted = True
        self.__split_freqs = None
        self.__trees_counted_for_freqs = 0
        
    def add_split_count(self, split, count):
        if split not in self.splits:
            self.splits.append(split)
            self.split_counts[split] = 0
        self.split_counts[split] += count
        
    def splits_considered(self):
        """
        Returns 4 values:
            total number of splits counted
            total number of unique splits counted
            total number of non-trivial splits counted
            total number of unique non-trivial splits counted
        """
        if not self.split_counts:
            return 0, 0, 0, 0
        num_splits = 0
        num_unique_splits = 0
        num_nt_splits = 0
        num_nt_unique_splits = 0
        taxa_mask = self.taxa_block.all_taxa_bitmask()
        for s in self.split_counts:
            num_unique_splits += 1 
            num_splits += self.split_counts[s]
            if is_non_singleton_split(s, taxa_mask):
                num_nt_unique_splits += 1
                num_nt_splits += self.split_counts[s]
        return num_splits, num_unique_splits, num_nt_splits, num_nt_unique_splits                                                        
        
    def calc_freqs(self):
        "Forces recalculation of frequencies."
        self.__split_freqs = {}
        if self.total_trees_counted == 0:
            total = 1
        else:
            total = self.total_trees_counted
        for split in self.split_counts:
            self.__split_freqs[split] = float(self.split_counts[split]) / total
        self.__trees_counted_for_freqs = self.total_trees_counted            
        return self.__split_freqs
        
    def _get_split_frequencies(self):
        "Returns dictionary of splits : split frequencies."
        if self.__split_freqs is None or self.__trees_counted_for_freqs != self.total_trees_counted:
            self.calc_freqs()
        return self.__split_freqs   
        
    split_frequencies = property(_get_split_frequencies)         

    def count_splits_on_tree(self, tree):
        """
        Counts splits in this tree and add to totals. `tree` must be decorated
        with splits, and no attempt is made to normalize taxa.
        """
        if self.taxa_block is None:
            self.taxa_block = tree.taxa_block
        else:
            assert tree.taxa_block is self.taxa_block
        self.total_trees_counted += 1

        for split, edge in tree.split_edges.iteritems():
            if not self.unrooted:
                split = edge.clade_mask
            if split not in self.split_counts:
                self.splits.append(split)
                self.split_counts[split] = 1
            else:
                self.split_counts[split] += 1     
            if not self.ignore_edge_lengths:
                sel = self.split_edge_lengths.setdefault(split,[])
                if edge.length is not None: 
                    sel.append(tree.split_edges[split].length)
                # for correct behavior when some or all trees have no edge lengths
#                 else:
#                     self.split_edge_lengths[split].append(0.0)
            if not self.ignore_node_ages:
                sna = self.split_node_ages.setdefault(split, [])
                if edge.head_node is not None:
                    sna.append(edge.head_node.distance_from_tip())

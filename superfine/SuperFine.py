'''
    This module contains code for the main SuperFine algorithm.
'''

###########################################################################
##    Copyright 2010 Rahul Suri and Tandy Warnow.
##    This file is part of SuperFine.
##
##    SuperFine is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    SuperFine is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with SuperFine.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################

########################################################################################################################
#   Copyright 2016 Diogo Telmo Neves and Tandy Warnow.                                                                 #
#   This file is part of SuperFine and was adapted from the baseline implementation (see above).                       #
#                                                                                                                      #
#   Besides some refactoring, features to support running maximum likelihood (ML) analyses were added.                 #
#                                                                                                                      #
#   The license is exactly the same of the baseline implementation (see above).                                        #
########################################################################################################################

import os, sys, copy

from newick_modified.tree import *
from spruce.unrooted import *
from spruce.metrics import *
from superfine.adapters import *
from superfine.logger import *


def SuperFine(input, options):
    """Main SuperFine loop."""

    # Read phase/step
    sourceTrees = [parse_tree(sourceTree) for sourceTree in readMultipleTreesFromFile(input)]

    if options.writeData:
        (baseName, _, _) = input.rpartition(".")

    # SCM phase/step
    tree = mergeTrees(sourceTrees, options)

    if options.writeData:
        f = open(baseName + ".scmTree." + options.writeData, 'w')
        f.write(str(tree))
        f.write(';\n')
        f.close()

    # Refinement step/phase
    #
    # key: polytomy node, value: list of bipartitions below the polytomy to
    #   effect; each bipartition is represented as a set of polytomy subtrees,
    #   and the set of all other leaf labels in the tree implicitly represents
    #   the set of leaf labels on the other side of the bipartition
    bipartitionsToAdd = {}
    logger = Logger()

    if options.reconciler == "qmc":
        for polytomy in xfindPolytomies(tree):
            (quartetTrees, delabeling) = encodeSourceTrees(polytomy, sourceTrees, logger, options)

            if quartetTrees:    # not empty list of quartet trees with which to resolve polytomy
                quartetTrees = selectSubset(quartetTrees, options)
                bipartitionsToAdd[polytomy] = reconcileTrees(quartetTrees, delabeling, options)
    elif options.reconciler.endswith("mrp") or options.reconciler.endswith("fml") or options.reconciler.endswith("rml"):
        for polytomy in xfindPolytomies(tree):
            (newSourceTrees, delabeling) = relabelSourceTrees(polytomy, sourceTrees, logger, options)

            if newSourceTrees:  # there are new source trees with which to resolve polytomy
                bipartitionsToAdd[polytomy] = reconcileTrees(newSourceTrees, delabeling, options)

    # add new bipartitions to the original SCM tree
    expandTree(bipartitionsToAdd)

    # print output to stdout, diagnostic info to stderr
    if options.writeData:
        f = open(baseName + ".SuperFineTree." + options.writeData, 'w')
        f.write(str(tree))
        f.write(';\n')
        f.close()

    else:
        print(str(tree) + ';')
    #logger.printInfo()



def mergeTrees(sources, options):
    '''
        Merge source trees using the strict consensus merger.
    '''
    sourceTrees = copy.deepcopy(sources)
    scm = SCMAdapter(sourceTrees, "")
    mergerTree = addDegreeInfo(scm.get_tree())

    return (mergerTree)



def addQuartets(quartetTrees, newQuartetTrees):
    '''
        Add the new quartet trees to the list of quartet trees.
    '''
    for quartet in newQuartetTrees:
        quartetTrees[quartet] = quartetTrees.get(quartet, 0) + 1
    return (quartetTrees)



def relabelSourceTrees(polytomy, sourceTrees, logger, options):
    '''
        Relabel leaves in source trees with labels {1, ..., d} where d is the 
        degree of the polytomy.  Replace subtrees rooted at interior nodes whose
        children are all same-labelled with single so-labelled leaves.  Return 
        a list of such trees and a data structure mapping back to the original 
        labels. 
    '''

    newSourceTrees = []
    (relabeling, delabeling) = buildRelabeling(polytomy)

    for sourceTree in sourceTrees:
        relabeledSourceTree = relabelTree(sourceTree, relabeling)
        relabeledSourceTree = collapseTree(relabeledSourceTree)
        newSourceTrees.append(relabeledSourceTree)

    newSourceTrees = removeUninformativeRelabeledTrees(newSourceTrees)
    return(newSourceTrees, delabeling)



def encodeSourceTrees(polytomy, sourceTrees, logger, options):
    '''
        Encode all source trees into a single list of quartet trees.
    '''
    quartetTrees = {}

    (relabeling, delabeling) = buildRelabeling(polytomy)

    # find quartet trees
    for sourceTree in sourceTrees:
        relabeledSourceTree = relabelTree(sourceTree, relabeling)

        treeToEncode = collapseTree(relabeledSourceTree)

        qtrees = findDisplayedQtrees(treeToEncode)
        quartetTrees = addQuartets(quartetTrees, qtrees)
        del qtrees


    quartetTrees = removeUninformativeQTrees(quartetTrees)

    logger.logInfo(len(quartetTrees))
    return (quartetTrees, delabeling)



def selectSubset(quartetTrees, options):
    '''
        Select a subset of the given set of quartet trees.
    '''
    return (quartetTrees)



def reconcileTrees(trees, delabeling, options):
    '''
        Infer a tree from a set of quartet trees using QMC, MRP, or MRL
        as a black box subroutine.  Map this inferred tree to implied 
        bipartitions in the SCM merger tree using the given delabeling.
    '''
    if options.reconciler == "qmc": # QMC
        quartetTrees = trees
        reconciler = QMCAdapter(quartetTrees)
    elif options.reconciler.endswith("mrp"): # MRP
        sourceTrees = trees
        reconciler = MRPAdapter(sourceTrees, options.numIters, options.reconciler)
    elif options.reconciler.endswith("fml") or options.reconciler.endswith("rml"): # MRL
        sourceTrees = trees
        reconciler = MRLAdapter(sourceTrees, options.reconciler)
    else: # None
        pass

    tree = reconciler.get_tree()
    return findImpliedBipartitions(tree, delabeling)



def buildRelabeling(polytomy):
    '''
        Build a mapping from leaf labels to polytomy group numbers and a mapping
        from polytomy group numbers to polytomy subtrees.
    '''
    relabeling = {}
    delabeling = {}

    # assert: isNonLeaf(polytomy)
    children = [edge[0] for edge in polytomy.get_edges()]

    #for i in xrange(len(children)):
    for i in range(len(children)):
        child = children[i]
        for identifier in child.get_leaves_identifiers():
            relabeling[identifier] = i   # maps leaf labels to integers
        delabeling[i] = child            # maps integers to polytomy subtrees
    return(relabeling, delabeling)



def relabelTree(source, relabeling):
    '''
        Relabel the leaves of a source tree with the given relabeling.
    '''
    sourceTree = copy.deepcopy(source)
    defaultLabel = len(set(relabeling.values()))

    for leaf in sourceTree.get_leaves():
        leaf.set_leaf_identifier(relabeling.get(leaf.identifier, # key
                                                defaultLabel))   # default value
    return(sourceTree)



def collapseTree(sourceTree):
    '''
        Collapse same-labeled sibling leaves into single nodes until a "minimal"
        tree is reached.
    '''

    class SiblingMerger(TreeVisitor):
        '''
            Collapse into leaves any internal nodes whose leaf descendants are 
            all same-labeled.
        '''
        def post_visit_edge(self, src, bootstrap, length, dest):
            dest._leaves_cache = None
            leaves = dest.get_leaves_identifiers()

            if (len(leaves) > 1 and len(set(leaves)) == 1):
                index = src._edges.index((dest, bootstrap, length))
                src._edges[index] = ((Leaf(leaves[0]), bootstrap, length))
                src._leaves_cache = None

    # collapsing
    sourceTree.dfs_traverse(SiblingMerger())
    sourceTree._leaves_cache = None
    sourceTree = pruneOvergrowth(sourceTree)

    # assert(len(sourceTree.get_leaves()) == len(sourceTree.get_leaves_identifiers()))
    return (sourceTree)



def pruneOvergrowth(tree):
    '''
        When a tree is (arbitrarily) rooted at a node within polytomy group X,
        the SiblingMerger class defined in collapseTree() often induces a tree
        with an "overgrowth" of X-labeled leaves. This function defines a
        workaround by pruning such leaves.
    '''

    class TreePruner(TreeVisitor):
        '''Find the subtree whose leaves match the input set.'''
        def __init__(self, nodes):
            self.nodes = set(nodes)
            self.tree = None

        def post_visit_tree(self, tree):
            leaves = tree.get_leaves_identifiers()
            if set(leaves) == self.nodes:
                self.tree = tree

    # calculate frequency counts for leaf identifiers
    leaves = tree.get_leaves_identifiers()
    frequency = {}
    for leaf in leaves:
        frequency[leaf] = frequency.get(leaf, 0) + 1

    # case: source trees to be collapsed to less than 4 leaves are uninformative
    if (len(frequency.keys()) < 4):
        return Tree()

    # find overrepresented leaf labels
    singletons = [leafLabel for leafLabel in frequency.keys() if frequency[leafLabel] == 1]
    duplicates = filter(lambda x : x not in singletons, frequency.keys())
    # assert len(duplicates) <= 1

    # prune away "overgrowth"
    tp = TreePruner(singletons)
    tree.dfs_traverse(tp)
    
    if (tp.tree):
        # add a single root-adjacent leaf for each pruned-away leaf label
        for leafLabel in duplicates:
            tp.tree.add_edge((Leaf(leafLabel), None, None))

        return (tp.tree)
    else:
        return (tree)



def removeUninformativeRelabeledTrees(sourceTrees):
    '''Remove trees with less than four unique taxa.'''
    treesToRemove = []
    for tree in sourceTrees:
        if len(set(tree.get_leaves_identifiers())) < 4:
            treesToRemove.append(tree)

    for tree in treesToRemove:
        sourceTrees.remove(tree)

    return(sourceTrees)



def removeUninformativeQTrees(quartetTrees):
    '''
        Remove quartets which lack 4 unique leaves.
    '''
    treesToRemove = []

    for qTree in quartetTrees.keys():
        taxa = qTree.replace('|', ',').split(',')

        if len(set(taxa)) != 4:
            treesToRemove.append(qTree)

    for qTree in treesToRemove:
        del quartetTrees[qTree]

    return(quartetTrees)



def findImpliedBipartitions(reconciledTree, delabeling):
    '''
        Find the bipartitions implied by the subroutine. Return as a list of
        lists of subtrees.
    '''
    ####################################################################################################################
    # inserted by: diogo telmo neves - dneves@di.uminho.pt
    #
    import re
    reconciledTree = re.sub(r"\)[0-9]+(\.[0-9]+)?", ")", reconciledTree, flags=re.MULTILINE)
    #
    ####################################################################################################################

    bipartitions = []
    tree = parse_tree(reconciledTree)

    for (src, dest) in getInternalEdges(tree):
        # find bipartition induced by the given edge
        setA = [label for label in dest.get_leaves_identifiers()]
        setB = [label for label in tree.get_leaves_identifiers() if label not in setA]

        # Parent of the polytomy (if one exists) is always given by the largest
        #   label number. We want to store the bipartition leaf set which does
        #   *not* include it.
        identifierIntegers = [int(leaf) for leaf in tree.get_leaves_identifiers()]
        maxElement = str(max(identifierIntegers))

        if maxElement not in setA:
            bipartitions.extend([[delabeling[int(label)] for label in setA]])
        else:
            bipartitions.extend([[delabeling[int(label)] for label in setB]])

    return (bipartitions)



def expandTree(bipartitionsToAdd):
    '''Add newly found bipartitions to the original SCM tree.'''

    for polytomy in bipartitionsToAdd.keys():
        bipartitions = bipartitionsToAdd[polytomy]

        # categorize bipartitions by size (i.e. # of subtrees)
        sizes = {}
        for bipartition in bipartitions:
            length = len(bipartition)
            #if (sizes.has_key(length)):
            if (length in sizes.keys()):
                sizes[length].append(bipartition)
            else:
                sizes[length] = [bipartition]

        # iterate over bipartitions, smallest to largest
        for size in sorted(sizes.keys()):
            for bipartition in sizes[size]:

                # find leaves which must be positioned on one side of the new
                #   bipartition
                bipartitionLeaves = []
                for subtree in bipartition:
                    bipartitionLeaves.extend(subtree.get_leaves_identifiers())

                # find branches to relocate below the new edge (i.e. all
                #   polytomy-incident branches s.t. all of their leaf
                #   descendants are in bipartitionLeaves)
                branchesToMove = []
                for edge in polytomy.get_edges():
                    problemLabels = [label for label in edge[0].get_leaves_identifiers()
                                           if label not in bipartitionLeaves]
                    if len(problemLabels) == 0:
                        branchesToMove.append(edge)

                # construct subtree to be positioned below the new edge
                newSubtree = Tree()
                for edge in branchesToMove:
                    newSubtree.add_edge(edge)
                    polytomy._edges.remove(edge)

                # attach new subtree to the polytomy node
                polytomy.add_edge((newSubtree, None, None))

# coding=ISO-8859-1
# -*- coding: ISO-8859-1 -*-

'''This module contains code which interfaces SuperFine with other software.'''

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

import sys
import tempfile
import random
import os
from subprocess import Popen, PIPE
from dendropy.dataio import trees_from_newick
from dendropy.scripts.strict_consensus_merge import strict_consensus_merge
from spruce.mrp import matrixRepresentation, writeRatchetInputFile, getConsensusTreesFromPaupFiles, readTreesFromRatchet
from spruce.unrooted import readNewickFile
from newick_modified.tree import Tree, parse_tree
from matrix_representation.MatrixRepresentation import MatrixRepresentation


def call_command(command, input):
    """Call the command as a subprocess with the given input, and return output and error streams."""
    try:
        pipe = Popen(command, stdin = PIPE, stdout = PIPE, stderr = PIPE)
        (output, err) = pipe.communicate(input)
        return output, err

    except OSError:
        print("Execution of %s failed" % command)
        sys.exit(1)


class SCMAdapter(object):
    """This class is an adapter for the strict consensus merger (SCM) functionality provided by DendroPy."""

    def __init__(self, trees, mergerType):
        self.trees = trees
        # remove uninformative trees
        uninformatives = []
        for t in self.trees:
            if len(set(t.get_leaves_identifiers())) < 4:
                uninformatives.append(t)
        for t in uninformatives:
            self.trees.remove(t)

        self.useGordons = False

        if mergerType == "gordons":
            self.useGordons = True

    def getLeafSets(self):
        leafSets = [set(tree.get_leaves_identifiers()) for tree in self.trees]
        return (leafSets)

    def getOverlap(self, leafSet1, leafSet2):
        intersectionSize = len(leafSet1 & leafSet2)
        return (intersectionSize)

    def getNextPair(self):
        '''
            Return the pair of trees whose leaf sets have the largest             
            intersection; break ties arbitrarily.
        '''
        numTrees = len(self.trees)
        if numTrees == 1:
            raise ValueError("SCM called with only 1 input tree.")

        (index1, index2, currentMax) = (None, None, 3)
        leafSets = self.getLeafSets()

        for i in range(numTrees - 1):
            for j in range(i + 1, numTrees):
                overlap = self.getOverlap(leafSets[i], leafSets[j])

                if (overlap > currentMax):
                    (index1, index2, currentMax) = (i, j, overlap)

        if currentMax == 3: # insufficient overlap for merger
            raise ValueError("Insufficient overlap for SCM step (%d trees left)" % numTrees)
        else:
            return (self.trees[index1], self.trees[index2])

    def pairwiseMerger(self, tree1, tree2):
        data = trees_from_newick((tree1, tree2))
        trees = [i[0] for i in data.trees_blocks]
        output = strict_consensus_merge(trees, gordons_supertree=self.useGordons)
        return (str(output))

    def get_tree(self):
        numMergers = len(self.trees) - 1
        for i in range(numMergers):
            (tree1, tree2) = self.getNextPair()
            self.trees.remove(tree1)
            self.trees.remove(tree2)

            newTree = self.pairwiseMerger(str(tree1), str(tree2))
            self.trees.append(parse_tree(newTree))

        # assert: len(self.trees) == 1
        return (self.trees[0])


class QMCAdapter(object):
    """This class is an adapter for supertree construction functionality provided by Quartets MaxCut (QMC)."""

    def __init__(self, quartetTrees):
        self.trees = quartetTrees

    def get_tree(self):
        inputs = ["%s:%s" % (self.trees[qTree], qTree) for qTree in sorted(self.trees.keys())]
        input = '\n'.join(inputs)
        (output, err) = call_command ("find-cut", input)
        return (output)


class MRPAdapter(object):
    """This class is an adapter for supertree construction functionality using MRP provided by PAUP*."""

    def __init__(self, sourceTrees, numIters = 100, mrpType = 'gmrp'):
        self.trees = sourceTrees
        self.numIters = numIters
        self.mrpType = mrpType

    def get_tree(self):
        matrix = matrixRepresentation(self.trees)
        if (len(matrix[0]) == 0 or len(matrix[1]) == 0):
            return (str(Tree()))

        f = tempfile.NamedTemporaryFile()
        prefix = f.name
        f.close()

        f = open (prefix, 'w')
        writeRatchetInputFile(matrix, f, filePrefix = prefix, numRatchetIterations = self.numIters)
        f.close()

        pipe = Popen("paup -n %s" % prefix, shell = True, stdout = PIPE, stderr = PIPE)
        (out, err) = pipe.communicate()

        if self.mrpType == 'gmrp':
            trees = getConsensusTreesFromPaupFiles(prefix)
            output = trees['gmrp']

        elif self.mrpType == 'rmrp':
            mpTrees = readTreesFromRatchet(prefix + '.tre')
            output = random.choice(mpTrees)

        os.remove(prefix)
        os.remove(prefix + ".log")
        os.remove(prefix + ".gmrp")
        os.remove(prefix + ".smrp")
        os.remove(prefix + ".mmrp")
        os.remove(prefix + ".tre")
        os.remove(prefix + ".tre.nex")

        return output


class MRLAdapter(object):
    """
    This class is an adapter for supertree construction functionality
    using maximum likelihood (ML) methods [1][2][3][4].


    References:
    [1] Nguyen, Nam and Mirarab, Siavash and Warnow, Tandy. MRL and SuperFine+MRL: new supertree methods.
        Algorithms for Molecular Biology (AMB), 2012.
    [2] Price, M.N., Dehal, P.S., and Arkin, A.P. (2009) FastTree: Computing Large Minimum-Evolution Trees
        with Profiles instead of a Distance Matrix. Molecular Biology and Evolution, 2009 26:1641-1650.
    [3] A. Stamatakis. RAxML version 8: A tool for phylogenetic analysis and post-analysis of large phylogenies.
        Bioinformatics, 2014.
    [4] Neves, Diogo Telmo and Sobral, João Luís.
        Parallel SuperFine - A Tool for Fast and Accurate Supertree Estimation: Features and limitations.
        Future Generation Computer Systems, 2016.
    """
    __author__ = "dneves@di.uminho.pt"
    __date__ = "$May 27, 2013 9:37:49 AM$"

    def __init__(self, source_trees, method="fml"):
        """
        Creates an instance of ``MRLAdapter``.

        :param source_trees: The source trees
        :param method: The chosen maximum likelihood (ML) method (default: fml --> run FastTree analysis)
        """
        self.source_trees = source_trees
        self.method = method

    def get_tree(self):
        """
        Returns a ``newick_modified.tree.Tree`` that is computed from the source trees (``self.source_trees``).
        Firstly, a matrix representation (i.e. a supermatrix) is computed from the source trees.
        Then, the chosen maximum likelihood (ML) method is applied over the supermatrix to infer the tree.

        :return: A ``newick_modified.tree.Tree`` that is computed from the source trees (``self.source_trees``).
        """
        supermatrix = MatrixRepresentation(self.source_trees)

        if supermatrix.number_of_sites:
            # just to get the temporary filename
            f = tempfile.NamedTemporaryFile()
            filename = f.name
            f.close()

            # save the supermatrix to a temporary file
            f = open(filename + ".mr", 'w')
            f.write(supermatrix.phylip_format)
            f.flush()
            f.close()

            if self.method == "rml":    # RAxML
                ########################################################################################################
                # TODO: Provide RAxML support                                                                          #
                tree = Tree()   # for now a dumb tree is returned...                                                   #
                ########################################################################################################
            else:   # fml --> FastTree (default ML method)
                pipe = Popen("FastTree -gtr -nosupport -nt < {0}.mr > {0}.tmp".format(filename),
                             shell=True, stdout=PIPE, stderr=PIPE)
                (out, err) = pipe.communicate()
                tree = readNewickFile(filename + ".tmp")

            # house cleaning...
            try:
                os.remove(filename + ".mr")
                filename_reduced = filename + ".mr.reduced"
                # file can be created under certain circumstances
                if os.path.exists(filename_reduced):
                    os.remove(filename_reduced)
                filename_trees = filename + ".trees"
                if os.path.exists(filename_trees):
                    os.remove(filename_trees)
            except OSError:
                sys.stderr.write("MRLAdapter --> get_tree(): Something went wrong while trying to remove files.")

        else:
            tree = Tree()

        tree_unweighted = tree.getUnweighted()
        if not tree_unweighted.endswith(";"):
            tree_unweighted += ";"

        return tree_unweighted

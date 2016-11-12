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

import sys, os, tempfile, random
from subprocess import Popen, PIPE

from dendropy.dataio import trees_from_newick
from dendropy.scripts.strict_consensus_merge import strict_consensus_merge
from newick_modified.tree import parse_tree
from spruce.mrp import *


class SCMAdapter(object):
    '''
        This class is an adapter for the strict consensus merger (SCM) 
        functionality provided by DendroPy.
    '''

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
        else :
            return (self.trees[index1], self.trees[index2])


    def pairwiseMerger(self, tree1, tree2):
        data = trees_from_newick((tree1, tree2))
        trees = [i[0] for i in data.trees_blocks]
        output = strict_consensus_merge(trees, gordons_supertree=self.useGordons)
        return (str(output))


    def getTree(self):
        numMergers = len(self.trees) - 1
        for i in range(numMergers):
            (tree1, tree2) = self.getNextPair()
            self.trees.remove(tree1)
            self.trees.remove(tree2)

            newTree = self.pairwiseMerger(str(tree1), str(tree2))
            self.trees.append(parse_tree(newTree))

        # assert: len(self.trees) == 1
        return (self.trees[0])



def callCommand (command, input):
    '''
        Call the command as a subprocess with the given input, and return output 
        and error streams.
    '''
    try:
        pipe = Popen(command, stdin = PIPE, stdout = PIPE, stderr = PIPE)
        (output, err) = pipe.communicate(input)
        return (output, err)

    except OSError:
        print("Execution of %s failed" % (command))
        sys.exit(1)



class QMCAdapter(object):
    '''
        This class is an adapter for supertree construction functionality 
        provided by Quartets MaxCut (QMC).
    '''

    def __init__(self, quartetTrees):
        self.trees = quartetTrees

    def getTree(self):
        inputs = ["%s:%s" % (self.trees[qTree], qTree) for qTree in sorted(self.trees.keys())]
        input = '\n'.join(inputs)
        (output, err) = callCommand ("find-cut", input)
        return (output)



class MRPAdapter(object):
    '''
        This class is an adapter for supertree construction functionality
        using MRP provided by PAUP*.
    '''

    def __init__(self, sourceTrees, numIters = 100, mrpType = 'gmrp'):
        self.trees = sourceTrees
        self.numIters = numIters
        self.mrpType = mrpType


    def getTree(self):
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

        return(output)



class MRLAdapter(object):
    """
    This class is an adapter for supertree construction functionality using maximum likelihood methods.
    """

    def __init__(self, sourceTrees, method="fml"):
        self.source_trees = sourceTrees
        self.method = method


    def getTree(self):
        if self.method == "rml": # RAxML
            pass
        else: # FastTree
            pass

        tempfile.tempdir = self.tempDirectory
        tempfile.gettempdir()
        f = tempfile.NamedTemporaryFile()
        prefix = f.name
        f.close()

        return None

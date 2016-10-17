'''
    This module defines several error metrics.
'''

###########################################################################
##    Copyright 2010 Rahul Suri and Tandy Warnow.
##    This file is part of spruce.
##
##    spruce is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    spruce is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with spruce.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################

import os
import tempfile
from subprocess import Popen, PIPE

from spruce.unrooted import *
from spruce.mrp import matrixRepresentation
from spruce.mrp import writeMatrix


def getFpFnRfRates (trueTree, estimatedTree, twoWay = False):
    '''
        Return FP, FN, and Robinson-Foulds rates, given true & estimated trees. 
        If the (optional) third argument is set to True, the trees will first 
        be restricted to the intersection of their taxa; otherwise, an exception
        will be raised if their leaf sets are not identical.
    '''

    if set(trueTree.get_leaves_identifiers()) != set(estimatedTree.get_leaves_identifiers()):
        if twoWay:
            intersection = set(trueTree.get_leaves_identifiers()) & set(estimatedTree.get_leaves_identifiers())
            trueTree = restrict(trueTree, intersection)
            estimatedTree = restrict(estimatedTree, intersection)
        else:
            raise Exception('Leaf sets are not identical')

    trueBPs = set(xfindBipartitions(trueTree))
    estimatedBPs = set(xfindBipartitions(estimatedTree))

    falsePositives = estimatedBPs - trueBPs
    falseNegatives = trueBPs - estimatedBPs

    if len(estimatedBPs) == 0:
        fpRate = None
    else:
        fpRate = 1.0*len(falsePositives)/len(estimatedBPs)

    if len(trueBPs) == 0:
        fnRate = None
    else:
        fnRate = 1.0*len(falseNegatives)/len(trueBPs)

    if fpRate == None or fnRate == None:
        rfRate = None
    else:
        rfRate = (fpRate + fnRate)/2.0

    return(fpRate, fnRate, rfRate)



def getRawFpFn (trueTree, estimatedTree, twoWay = False):
    '''
        Return actual numbers of FP,FN, and TP, *not* rates, given true & 
        estimated trees.  If the (optional) third argument is set to True, the 
        trees will first be restricted to the intersection of their taxa; 
        otherwise, an exception will be raised if their leaf sets are not 
        identical.
    '''

    if set(trueTree.get_leaves_identifiers()) != set(estimatedTree.get_leaves_identifiers()):
        if twoWay:
            intersection = set(trueTree.get_leaves_identifiers()) & set(estimatedTree.get_leaves_identifiers())
            trueTree = restrict(trueTree, intersection)
            estimatedTree = restrict(estimatedTree, intersection)
        else:
            raise Exception('Leaf sets are not identical.')

    trueBPs = set(xfindBipartitions(trueTree))
    estimatedBPs = set(xfindBipartitions(estimatedTree))

    falsePositives = estimatedBPs - trueBPs
    falseNegatives = trueBPs - estimatedBPs

    return(len(falsePositives), len(falseNegatives), len(trueBPs), len(estimatedBPs))



def sumErrorsAcrossSources(sourceTrees, tree, normalized = False, twoWay = False):
    '''
        Return the sums of FP, FN, and Robinson-Foulds scores of the tree 
        restricted to the leaf set of each source tree and scored against that 
        source tree.  If the (optional) third argument is set to True, these 
        scores are normalized.  If the (optional) fourth argument is set to 
        True, each source tree is also restricted to the estimated tree's leaf 
        set prior to scoring; otherwise, an exception is raised if the source 
        tree contains taxa that the estimated tree doesn't.
    '''

    sum = {}
    sum["fp"] = 0
    sum["fn"] = 0
    sum["rf"] = 0
    sum["trueBPs"] = 0
    sum["estimatedBPs"] = 0

    for source in sourceTrees:
        restrictedTree = restrict(tree, set(source.get_leaves_identifiers()))

        (fp, fn, trueBPs, estimatedBPs) = getRawFpFn (source, restrictedTree, twoWay)
        rf = fp + fn

        sum["fp"] += fp
        sum["fn"] += fn
        sum["rf"] += rf
        sum["trueBPs"] += trueBPs
        sum["estimatedBPs"] += estimatedBPs

    if normalized:
        if sum["estimatedBPs"] == 0:
            normalizedSumFP = None
        else: 
            normalizedSumFP = 1.0*sum["fp"]/sum["estimatedBPs"]

        if sum["trueBPs"] == 0:
            normalizedSumFN = None
            normalizedSumRF = None 
        else:
            normalizedSumFN = 1.0*sum["fn"]/sum["trueBPs"]
            normalizedSumRF = 1.0*sum["rf"]/(2*sum["trueBPs"])

        return (normalizedSumFP, normalizedSumFN, normalizedSumRF)

    else:
        return(sum["fp"], sum["fn"], sum["rf"])



def getResolution (tree):
    '''Return the percent resolution of the tree.'''

    numBipartitions = len(set(xfindBipartitions(tree)))
    numTaxa = len(set(tree.get_leaves_identifiers()))
    resolution = (1.0 * numBipartitions)/(numTaxa - 3)

    return (resolution)



def getParsimonyScores(sourceTrees, supertrees, rooted = False):
    '''
        Return the parsimony scores of the supertrees, given source trees.  
        Note: PAUP* must be installed and executable with a call to "paup -n" 
        for this function.
    '''

    rooting = 'U'
    if (rooted):
        rooting = 'R'

    f = tempfile.NamedTemporaryFile()
    tempName = f.name
    f.close()

    matrix = matrixRepresentation(sourceTrees)
    f = open (tempName, 'w')
    writeMatrix(matrix, f)
    del matrix

    f.write("set maxtrees = %d;\n" % len(supertrees))
    f.write("Begin trees;\n")
    for i in xrange(len(supertrees)):
        f.write("\ttree %d = [&%s] %s;\n" % (i, rooting, supertrees[i]))
    f.write("end;\n\n")

    instructions = ['begin paup;',
                    'set criterion = parsimony;',
                    'pscores all/ scorefile = %s.out.nex replace = yes;' % tempName,
                    'quit;',
                    'end;', 
                    'quit warntsave = no;']
    f.write('\n'.join(instructions))
    f.close()

    pipe = Popen("paup -n %s" % tempName, shell = True, stdout = PIPE, stderr = PIPE)
    (out, err) = pipe.communicate()

    f = open(tempName + ".out.nex", 'r')
    lines = f.readlines()
    f.close()

    scores = []

    for line in lines[1:]:
        scores.append(int(line.split()[1]))

    os.remove(tempName)
    os.remove(tempName + ".out.nex")

    return (scores)

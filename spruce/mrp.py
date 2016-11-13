'''
    This module defines functionality for working with matrices of binary 
    character encodings of trees.  Some of the functionality relies on having 
    PAUP* installed and accessible via a call of the form "paup -n".
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

from newick_modified.tree import *
from spruce.unrooted import *

import random
import os
import tempfile
from subprocess import Popen, PIPE

def matrixRepresentation (trees):
    '''
        Return a matrix of binary characters representing the set of trees.  The
        matrix returned is a sparse representation: a 2-tuple containing a list 
        of taxa and a list of mappings (one mapping for each character) from 
        taxon name to value. Unknown values (typically '?' in Nexus files) are 
        represented implicitly.
    '''

    taxa = tuple(sorted(list(set([identifier for tree in trees
                                             for identifier in tree.get_leaves_identifiers()]))))
    columns = []

    for tree in trees:
        for bp in xfindBipartitions(tree):
            columns.append({})
            currentColumn = columns[-1]

            for taxon in taxa:
                if taxon in bp[0] and taxon in bp[1]:
                    continue
                elif taxon in bp[0]:
                    currentColumn[taxon] = 0
                elif taxon in bp[1]:
                    currentColumn[taxon] = 1

    return (taxa, columns)



def readMatrixFromFile (file):
    '''
        Read a matrix (into the data structure described above) from the 
        specified Nexus file.
    '''

    f = open (file, 'r')
    lines = f.readlines()
    f.close()

    (startLine, endLine) = (None, None)

    # find sequence lines
    #for i in xrange(len(lines)):
    for i in range(len(lines)):
        if lines[i].strip().lower() == 'matrix':
            startLine = i + 1
        elif lines[i].strip() == ';':
            endLine = i - 1

    # find number of characters
    for line in lines[startLine: endLine + 1]:
        if line.strip() == '':
            continue
        (taxon, sequence) = line.split()
        sequenceLength = len(sequence)
        #columns = [{} for i in xrange(sequenceLength)]
        columns = [{} for i in range(sequenceLength)]
        taxa = []
        break

    # populate data structures
    for line in lines[startLine: endLine + 1]:
        if line.strip() == '':
            continue
        (taxon, sequence) = line.split()
        taxa.append(taxon)

        #zeroes = [i for i in xrange(len(sequence)) if sequence[i] == '0']
        zeroes = [i for i in range(len(sequence)) if sequence[i] == '0']
        for index in zeroes:
            columns[index][taxon] = 0

        #ones = [i for i in xrange(len(sequence)) if sequence[i] == '1']
        ones = [i for i in range(len(sequence)) if sequence[i] == '1']
        for index in ones:
            columns[index][taxon] = 1

    taxa = tuple(sorted(taxa))
    return(taxa, columns)



def writeMatrix (matrix, fileStream, format = "nexus"):
    '''
        Write the given matrix to the given stream.  Default format is Nexus, 
        though Phylip is also supported via the (optional) third parameter.
    '''

    (taxa, columns) = (matrix[0], matrix[1])
    numTaxa = len(taxa)
    numChars = len(columns)
    missingLabel = '?'

    if format == "nexus":
        fileStream.write("#NEXUS\nbegin characters;\n")
        fileStream.write("\tdimensions newtaxa ntax = %d nchar = %d;\n" % (numTaxa, numChars))
        fileStream.write("\tformat missing = %s;\n" % missingLabel)
        fileStream.write("\tmatrix\n")

    elif format == "phylip":
        fileStream.write("\t%d\t%d\n" % (numTaxa, numChars))

    for taxon in taxa:
        values = [str(column.get(taxon, missingLabel)) for column in columns]
        values = ''.join(values)

        fileStream.write("\t'%s'\t" % str(taxon))
        fileStream.write(values)
        fileStream.write("\n")

    if format == "nexus":
        fileStream.write(";\nend;\n\n")



def getMatrixString (matrix, format = "nexus"):
    '''
        Return the given matrix as a string.  Default format is Nexus, though 
        Phylip is also supported via the (optional) second parameter.
    '''

    (taxa, columns) = (matrix[0], matrix[1])
    numTaxa = len(taxa)
    numChars = len(columns)
    missingLabel = '?'

    lines = []

    if format == "nexus":
        ################################################################################################################
        # changed by: diogo telmo neves - dneves@di.uminho.pt
        #
        lines.append("#NEXUS")
        lines.append("begin taxa;")
        lines.append("\tdimensions ntax={};".format(numTaxa))
        lines.append("\ttaxlabels {};".format(' '.join([str(taxon) for taxon in taxa])))
        lines.append("end;")
        lines.append("begin characters;")
        lines.append("\tdimensions newtaxa ntax = %d nchar = %d;" % (numTaxa, numChars))
        lines.append("\tformat missing = %s;" % missingLabel)
        lines.append("\tmatrix")
        #
        ################################################################################################################

    elif format == "phylip":
        lines.append("\t%d\t%d" % (numTaxa, numChars))

    for taxon in taxa:
        values = [str(column.get(taxon, missingLabel)) for column in columns]
        values = ''.join(values)

        lines.append("\t'%s'\t%s" % (str(taxon), values))

    if format == "nexus":
        lines.append(";\nend;\n")

    return '\n'.join(lines)



def writeRatchetInputFile (matrix, 
                           fileStream,
                           filePrefix = "ratchet", 
                           numRatchetIterations = 100, 
                           percentToUpweight = .25, 
                           weight = 2, 
                           startingTree = None):
    '''
        Write the given matrix and PAUP* commands for performing a ratcheted MRP
        analysis to the given stream.  The third parameter sets the prefix for 
        names of files generated by PAUP* upon running the file.  Other 
        parameters determine the details of the ratchet analysis.
    '''

    (taxa, columns) = (matrix[0], matrix[1])
    numTaxa = len(taxa)
    numChars = len(columns)

    numCharactersToSelect = int(numChars * percentToUpweight)
    randomSeed = random.randint(0, 10000)

    logFile = filePrefix + ".log"
    treeFile = filePrefix + ".tre"

    strictConsensusTreeFile = filePrefix + ".smrp"
    majorityConsensusTreeFile = filePrefix + ".mmrp"
    greedyConsensusTreeFile = filePrefix + ".gmrp"

    lines = getMatrixString(matrix)

    treeBlock = []
    if (startingTree):
        treeBlock = ['begin trees;',
                     '\ttree 1 = [&U] %s;' % startingTree,
                     'end;']

    ratchetInfo = ['[Ratchet parameters: NumChar = %d' % numChars,
                   ' Inclusion Percentage = %f (%d characters)' % (percentToUpweight, numCharactersToSelect), 
                   ' Replicates = %d' % numRatchetIterations,
                   ' Final search = no]\n']

    paupBlock = ['begin paup;',
                 '\tset autoclose = yes warntree = no warnreset = no notifybeep = no monitor = yes taxlabels = full;',
                 '\tlog file = %s replace;' % logFile,
                 '\tset criterion = parsimony;',
                 '\tpset collapse = no;',
                 '\n\t[!][!*** Replicate 0 (initial tree) ***]']

    if not startingTree:
        paupBlock += ['\thsearch addseq = random nreps = 1 rseed = %d swap = TBR multrees = no dstatus = 60;' % randomSeed]

    paupBlock += ['\tsavetrees file = %s format = altnex replace;' % treeFile,
                  '\tsavetrees file = %s.nex format = nexus replace;\n' % treeFile]

    fileStream.write("%s\n%s\n%s\n%s" % (lines, 
                                         '\n'.join(treeBlock), 
                                         ' '.join(ratchetInfo), 
                                         '\n'.join(paupBlock)))

    listOfIndices = range(1, numChars + 1)

    #for replicate in xrange(numRatchetIterations):
    for replicate in range(numRatchetIterations):
        selectedIndices = map(str, sorted(random.sample(listOfIndices, numCharactersToSelect)))

        replicateBlock = ['\n\t[!][!*** Replicate #%d ***]' % (replicate + 1),
                          '\tweights %d: %s;' % (weight, ' '.join(selectedIndices)),
                          '\thsearch start = current swap = TBR multrees = no dstatus = 60;',
                          '\tweights 1: all;',
                          '\thsearch start = current swap = TBR multrees = no dstatus = 60;',
                          '\tsavetrees file = %s format = altnex append;' % treeFile,
                          '\tsavetrees file = %s.nex format = nexus append;\n' % treeFile]
        fileStream.write('\n'.join(replicateBlock))

    consensusBlock = ['\n\t[!][!*** Determining consensus trees ***]',
                      '\tset MaxTrees = %d;' % (2*numRatchetIterations + 1),
                      '\tgettrees file = %s allblocks = yes warntree = no;' % treeFile,
                      '\tset criterion = parsimony;',
                      '\tcondense collapse = no deldupes = yes;',
                      '\tfilter best = yes;',
                      '\tcontree all / strict = yes treefile = %s replace;' % strictConsensusTreeFile,
                      '\tcontree all / majrule = yes strict = no treefile = %s replace;' % majorityConsensusTreeFile,
                      '\tcontree all / majrule = yes strict = no le50 = yes treefile = %s replace;' % greedyConsensusTreeFile,
                      '\tsavetrees file = %s replace = yes format = altnex;' % treeFile,
                      '\n\tlog stop;',
                      'end;\n',
                      'quit warntsave = no;']
    fileStream.write('\n'.join(consensusBlock))



def readTreesFromRatchet (file):
    '''Read trees from a Nexus file.'''

    f = open (file, 'r')
    lines = f.readlines()
    f.close()

    trees = []
    for line in lines:
        line = line.strip()
        if line.startswith('tree') or line.startswith('Tree') or line.startswith('TREE'):
            (_, _, tree) = line.partition(']')
            tree = tree.strip()
            trees.append(tree)
    return (trees)



def getConsensusTreesFromPaupFiles (filePrefix):
    '''Read consensus trees generated by PAUP* as per commands above.'''

    consensuses = ['gmrp', 'mmrp', 'smrp']
    trees = {}
    for consensus in consensuses:
        f = open (filePrefix + "." + consensus, 'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            if line.startswith('\ttree'):
                (a, b, c) = line.partition('(')
                tree = b + c.strip()
                trees[consensus] = tree

    return (trees)



def getGreedyConsensus (sourceTrees, supertrees, rooted = False):
    '''
        Return greedy consensus of supertrees (actually only the plenary ones), 
        given source trees.  Note: this function requires PAUP* to be installed 
        and callable via a "paup -n" command.
    '''

    taxa = set([identifier for tree in sourceTrees for identifier in tree.get_leaves_identifiers()])
    plenarySupertrees = [tree for tree in supertrees if len(set(tree.get_leaves_identifiers())) == len(taxa)]
    supertrees = plenarySupertrees

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
    #for i in xrange(len(supertrees)):
    for i in range(len(supertrees)):
        f.write("\ttree %d = [&%s] %s;\n" % (i, rooting, supertrees[i]))
    f.write("end;\n\n")

    instructions = ['begin paup;',
                    'contree all / strict=no showtree=no majrule=yes le50 = yes treefile = %s.out.nex replace = yes;' % tempName,
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

    os.remove(tempName)
    os.remove(tempName + ".out.nex")

    return (lines)

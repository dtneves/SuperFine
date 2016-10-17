#!/usr/bin/env python

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

from optparse import OptionParser
from spruce.unrooted import *
from spruce.metrics import *

import sys

parser = OptionParser(usage = "usage: %prog [options]")
parser.add_option("-t", "--trees", dest = "estimate_tree", help = "read estimate trees from FILE", metavar = "FILE")
parser.add_option("-s", "--sources", dest = "source_trees", help = "read source trees from FILE", metavar = "FILE")

(options, args) = parser.parse_args()
if (options.estimate_tree == None or options.source_trees == None):
    parser.error("values for both options must be specified\ntry running with the --help flag")

estimates = [parse_tree (tree) for tree in readMultipleTreesFromFile(options.estimate_tree)]
sources = [parse_tree (tree) for tree in readMultipleTreesFromFile(options.source_trees)]

sizes = [len(estimate.get_leaves()) for estimate in estimates]
taxa = set([identifier for tree in sources for identifier in tree.get_leaves_identifiers()])

numUnique = len(set(sizes))
fullCount = len(taxa)

plenaryTrees = [tree for tree in estimates if len(tree.get_leaves()) == fullCount]
maxCount = len(plenaryTrees)

if numUnique > 1 or sizes[0] != fullCount:
    sys.stderr.write('WARNING: %d differnet sizes among %d trees, %d maximal (%d total taxa), %s\n' % (numUnique, len(sizes), maxCount, fullCount, options.estimate_tree))
    sys.stdout.write(';\n'.join(map(str, plenaryTrees)))
    sys.stdout.write(';\n\n') 

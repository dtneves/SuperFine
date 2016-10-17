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

desc = '''
           This script reads a tree and a set of source trees and prints the 
           Sum-FP, Sum-FN, and Sum-RF scores of the tree.
       '''

parser = OptionParser(usage = "usage: %prog [options]", description = desc)
parser.add_option("-t", "--tree", dest = "estimate_tree", help = "read estimate tree from FILE", metavar = "FILE")
parser.add_option("-s", "--sources", dest = "sources_file", help = "read source trees from FILE", metavar = "FILE")

(options, args) = parser.parse_args()
if (options.estimate_tree == None or options.sources_file == None):
    parser.error("values for both options must be specified\ntry running with the --help flag")

sources = [parse_tree(source) for source in readMultipleTreesFromFile(options.sources_file)]
estimate = readNewickFile(options.estimate_tree)

print sumErrorsAcrossSources(sources, estimate, normalized = True)

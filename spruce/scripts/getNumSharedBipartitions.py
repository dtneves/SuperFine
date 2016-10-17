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

import sys
from optparse import OptionParser
from spruce.unrooted import *
from spruce.metrics import *

desc = '''This script print the number of bipartitions common to both input trees.'''

parser = OptionParser(usage = "usage: %prog tree1 tree2", description = desc)

(options, args) = parser.parse_args()

if len(args) != 2:
    parser.error("Incorrect number of arguments. Try the -h flag for help.")

tree1 = readNewickFile(args[0])
tree2 = readNewickFile(args[1])

tree1_bipartitions = set(xfindBipartitions(tree1))
tree2_bipartitions = set(xfindBipartitions(tree2))

intersection = tree1_bipartitions & tree2_bipartitions

print len(intersection)

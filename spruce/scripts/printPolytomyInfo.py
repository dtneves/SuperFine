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

desc = '''
           This script reads an input tree and prints the degrees of polytomies 
           and the number of such occurrences in the tree.
       '''

parser = OptionParser(usage = "usage: %prog [options] < input > output", 
                      description = desc)
parser.add_option("-i", "--in", dest = "in_tree", help = "read tree from FILE", metavar = "FILE")
parser.add_option("-o", "--out", dest = "outfile", help = "write output to FILE", metavar = "FILE")

(options, args) = parser.parse_args()

if (options.in_tree == None):
    inputStream = sys.stdin
else:
    inputStream = open(options.in_tree, 'r')

if (options.outfile == None):
    outputStream = sys.stdout
else:
    outputStream = open(options.outfile, 'w')

data = inputStream.read()
if data.strip().endswith(':0.0;'):
    (data, _, _) = data.strip().rpartition(':0.0;')

tree = parse_tree(data)
addDegreeInfo(tree)

occurrences = {}
for polytomy in xfindPolytomies(tree):
    occurrences[polytomy.degree] = occurrences.get(polytomy.degree, 0) + 1

lines = ['%d\t%d' % (key, occurrences[key]) for key in sorted(occurrences.keys())]

outputStream.write('\n'.join(lines))
outputStream.write('\n')

inputStream.close()
outputStream.close()

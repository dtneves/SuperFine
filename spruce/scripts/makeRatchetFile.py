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
import sys
from newick_modified.tree import *
from spruce.unrooted import *
from spruce.mrp import *

desc = '''
           This script writes a Nexus file with a matrix of binary character 
           encodings of the trees read from input followed by PAUP* commands 
           for performing a ratcheted MRP analysis.
       '''

parser = OptionParser(usage = "usage: %prog [options] < input > output", 
                      description = desc)
parser.add_option("-i", "--input", dest = "input", metavar = "FILE",
                  help = "read Newick input from FILE")
parser.add_option("-o", "--output", dest = "output", metavar = "FILE",
                  help = "write Nexus output to FILE")
parser.add_option("-p", "--prefix", dest = "prefix", metavar = "STRING", 
                  default = "ratchet",
                  help = "use STRING as prefix for PAUP* output files, "
                         "[default: '%default']")

(options, args) = parser.parse_args()

if options.input == None:
    inputStream = sys.stdin
else:
    inputStream = open(options.input, 'r')

if options.output == None:
    outputStream = sys.stdout
else:
    outputStream = open(options.output, 'w')

sourceTrees = [parse_tree(sourceTree) for sourceTree in readMultipleTreesFromStream(inputStream)]
matrix = matrixRepresentation(sourceTrees)
writeRatchetInputFile(matrix, outputStream, filePrefix = options.prefix)

inputStream.close()
outputStream.close()

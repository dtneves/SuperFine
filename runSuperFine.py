#!/usr/bin/env python

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

from optparse import OptionParser, OptionGroup
from superfine.SuperFine import SuperFine


def parse_options(command_line=None):
    """Parse command line for options"""

    desc = "This script runs the superfine algorithm on a set of input trees given in a file, in Newick format."

    parser = OptionParser(usage="usage: %prog [options] input_trees_file > output",
                          version="%prog 1.0", description=desc)

    parser.set_defaults(reconciler="qmc", numIters=100, writeData=None)

    group4InfoString = "These options enable selection of the supertree algorithm " \
                       "to be used as a subroutine within superfine for resolving polytomies.  " \
                       "The 'qmc' option requires that Quartets MaxCut be installed; " \
                       "the 'gmrp' (greedy consensus of MP trees) and 'rmrp' (random MP tree) options " \
                       "require that PAUP* be installed; " \
                       "the 'fml' option requires that FastTree be installed; and " \
                       "the 'rml' option requires that RAxML be installed.  " \
                       "Note that the selected subroutine's binary must be in the system's executable search path."
    group4 = OptionGroup(parser, "Quartet Tree Reconciliation Options".upper(), group4InfoString)
    group4.add_option("-r", "--reconcile", choices=("qmc", "gmrp", "rmrp", "fml", "rml"),
                      dest="reconciler", metavar="ALG",
                      help="use ALG to reconcile relabeled trees, "
                           "where ALG is one of {qmc, gmrp, rmrp} [default: %default]")
    group4.add_option("-n", "--numIters", type="int", dest="numIters", metavar="N",
                      help="use N ratchet iterations when resolving with MRP [default: %default]")
    parser.add_option_group(group4)

    group5InfoString = ' '.join(["This option causes output of both the final",
                                 "superfine tree and the corresponding SCM tree",
                                 "to be written to disk.  Both are written to", 
                                 "the directory containing the source trees", 
                                 "file, and the specified suffix is appended to", 
                                 "each one's name."])
    group5 = OptionGroup(parser, "Data Output Options".upper(), group5InfoString)
    group5.add_option("-w", "--write", dest="writeData", metavar="SUFFIX",
                      help="write merger tree and final tree to disk in same "
                           "directory as source trees file, append .SUFFIX "
                           "to written file names [default: %default]")
    parser.add_option_group(group5)

    if command_line:
         (options, args) = parser.parse_args(command_line)
    else:
        (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error("Incorrect number of arguments. Try the -h flag for help.")

    input = args[0]

    return (input, options)


# MAIN
if __name__ == '__main__':
    (input, options) = parse_options()

    SuperFine(input, options)

'''This module contains logging utility functions for SuperFine.'''

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

import sys

class Logger(object):
    '''This class logs useful information.'''

    def __init__(self):
        self.unresolvablePolytomies = 0
        self.resolvablePolytomies = 0

    def logInfo(self, quartetTrees):
        '''Increment resolvable or unresolvable count.'''
        if not quartetTrees:
            self.unresolvablePolytomies += 1
        else:
            self.resolvablePolytomies += 1

    def printInfo(self):
        '''Print diagnostic info to stderr.'''
        # print info on resolvables
        if self.resolvablePolytomies == 1:
            sys.stderr.write("1 polytomy successfully resolved.\n")
        else:
            sys.stderr.write(self.resolvablePolytomies.__repr__() + " polytomies successfully resolved.\n")

        # print info on unresolvables
        if self.unresolvablePolytomies == 0:
            pass
        elif self.unresolvablePolytomies == 1:
            sys.stderr.write("1 polytomy could *not* be resolved.\n")
        else:
            sys.stderr.write(self.unresolvablePolytomies.__repr__() + " polytomies could *not* be resolved.\n")

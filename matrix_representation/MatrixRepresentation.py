#! /usr/bin/python

########################################################################################################################
#   Copyright 2013 Diogo Telmo Neves and Tandy Warnow.                                                                 #
#   This file is part of Matrix Representation.                                                                        #
#                                                                                                                      #
#   Matrix Representation is free software: you can redistribute it and/or                                             #
#   modify it under the terms of the GNU General Public License as published by                                        #
#   the Free Software Foundation, either version 3 of the License, or (at your option) any later version.              #
#                                                                                                                      #
#   Matrix Representation is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;                 #
#   without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                          #
#   See the GNU General Public License for more details.                                                               #
#                                                                                                                      #
#   You should have received a copy of the GNU General Public License along with Matrix Representation.                #
#   If not, see <http://www.gnu.org/licenses/>.                                                                        #
########################################################################################################################

__author__ = "dneves@di.uminho.pt"
__date__ = "$May 26, 2013 8:19:51 AM$"

import os
from spruce.unrooted import xfindBipartitions
from newick_modified.tree import Tree, parse_tree
from platform import system


class MatrixRepresentation(object):
    """
    This class allows one to get a matrix representation (aka supermatrix) from a set of overlapping source trees.
    See [1][2][3][4] to learn more about::
        (i) how to represent a set of overlapping source trees as a supermatrix; and
        (ii) the utility of creating such supermatrix.

    An example on how to use some of the functionality provided by this class is as follows::

    ``
    filename = "../datasets/biological/seabirds/kennedy.source_trees_manual"
    source_trees = source_trees_from_file(filename)
    supermatrix = MatrixRepresentation(source_trees)

    print(supermatrix)
    ``


    References:
    [1] Baum B. Combining trees as a way of combining data sets for phylogenetic inference, and
        the desirability of combining gene trees. Taxon. 1992;41:3-10.
    [2] Ragan MA. Phylogenetic inference based on matrix representation of trees. Mol Phylog Evol. 1992;1:53-58.
    [3] Swenson M, Suri R, Linder C, Warnow T. An experimental study of Quartets MaxCut and other supertree methods.
        Alg Mol Bio. 2011;6(7) Special issue for selected papers from WABI 2010.
    [4] Kupczok A, Schmidt H, von Haeseler A. Accuracy of phylogeny reconstruction methods
        combining overlapping gene data sets. Alg Mol Bio. 2010;5:37-53.
    """

    def __init__(self, source_trees, unknown_site_as_indel=False, nucleotide=False, validate=False):
        """
        Creates a new instance of ``MatrixRepresentation``.

        :param source_trees: A list of overlapping source trees from which
        a Matrix Representation (aka supermatrix) [1][2] will be computed,
        where each tree is an instance of ``newick.modified.tree.Tree``
        :param unknown_site_as_indel: Useful for tools (e.g. FastTree [3]) that give an awkward warning
        saying that the question mark symbol - ? - is an unknown symbol.
        Thus, if the given ``unknown_site_as_indel`` is ``False`` then unknown sites will be treated as expected
        (i.e. any unknown site is represented by a question mark symbol - ``?``),
        otherwise unknown sites will be treated as gaps (i.e. any unknown site is represented by an indel - ``-``)
        :param nucleotide: If ``False`` then a matrix over ``{0, 1, ?}`` or ``{0, 1, -}``
        will be computed accordingly to how the taxon appears on a given bipartition;
        otherwise a matrix over ``{A, C, ?}`` or over ``{A, C, -}``
        :param validate: Allows to bypass the validation of parameters ``source_trees``, ``unknown_site_as_indel``, and
        ``nucleotide``.  In other words, it allows to say:  - In user we trust!  Thus, performance is improved.


        References:
        [1] Baum B. Combining trees as a way of combining data sets for phylogenetic inference, and
            the desirability of combining gene trees. Taxon. 1992;41:3-10.
        [2] Ragan MA. Phylogenetic inference based on matrix representation of trees. Mol Phylog Evol. 1992;1:53-58.
        [3] Price, M.N., Dehal, P.S., and Arkin, A.P. (2009) FastTree: Computing Large Minimum-Evolution Trees
            with Profiles instead of a Distance Matrix. Molecular Biology and Evolution. 2009 26:1641-1650.
        """
        # bypass, or NOT, parameters validation
        if validate:
            # taking care...
            if not source_trees or not isinstance(source_trees, list):
                raise TypeError("The given \"source_trees\" has to be a non-empty list.")
            for source_tree in source_trees:
                if not isinstance(source_tree, Tree):
                    raise TypeError("Each element of the given \"source_trees\" list has to be "
                                    "an instance of (newick.modified.tree.)Tree.")
            if not isinstance(unknown_site_as_indel, bool):
                raise TypeError("The given \"unknown_site_as_indel\" is NOT an instance of Boolean.")
            if not isinstance(nucleotide, bool):
                raise TypeError("The given \"nucleotide\" is NOT an instance of Boolean.")

        self.__SOURCE_TREES = source_trees
        self.__UNKNOWN_SITE = '?' if (not unknown_site_as_indel) else '-'
        if not nucleotide:
            self.__ZERO_OR_A = '0'
            self.__ONE_OR_C = '1'
        else:
            # A - adenine nucleobase
            self.__ZERO_OR_A = 'A'
            # C - cytosine nucleobase
            self.__ONE_OR_C = 'C'
        self.__number_of_species = 0
        self.__number_of_sites = 0
        # the matrix representation of the given source trees
        self.__SUPERMATRIX = self.__compute_supermatrix()

    @staticmethod
    def source_trees_from_file(filename, validate=False):
        """
        Returns a list of phylogenetic trees, where each tree is an instance of ``newick.modified.tree.Tree``,
        after parsing each one contained in the file that has the given filename.
        It is expected that each phylogenetic tree, contained in the file, is represented in the Newick format [1][2].

        :param filename: The filename of a file containing phylogenetic trees,
        where each tree is represented in the Newick format [1][2].
        :param validate: Allows to bypass the validation of parameter ``filename``.
        In other words, it allows to say:  - In user we trust!  Thus, performance is improved.
        :return: A list of trees, where each tree is an instance of ``newick.modified.tree.Tree``


        References:
        [1] http://evolution.genetics.washington.edu/phylip/newicktree.html
        [2] https://en.wikipedia.org/wiki/Newick_format
        """
        if validate:
            if not filename or not isinstance(filename, str):
                raise TypeError("The given \"filename\" has to be a string representing a full path.")
            if not os.path.exists(filename):
                raise TypeError("The file \"{}\" does not exist.".format(filename))
            if os.path.isdir(filename):
                raise TypeError("The file \"{}\" is a directory.".format(filename))

        f = open(filename)
        source_trees_strs = f.readlines()
        f.close()

        return [parse_tree(source_tree_str) for source_tree_str in source_trees_strs]

    @staticmethod
    def supported_formats():
        """
        Returns a list with the supported file formats, which are:: FASTA, NEXUS, PHYLIP, and RAW.

        :return: A list with the supported file formats.
        """
        return ["FASTA", "NEXUS", "PHYLIP", "RAW"]

    def __compute_supermatrix(self):
        """
        Computes a matrix representation [1][2] from the source trees (aka supermatrix).

        The supermatrix is represented as a dictionary where::
            - each dictionary's key is a taxon, and
            - and the dictionary's value of a given dictionary's key - taxon - is a list of all sites of that taxon.


        References:
        [1] Baum B. Combining trees as a way of combining data sets for phylogenetic inference, and
            the desirability of combining gene trees. Taxon. 1992;41:3-10.
        [2] Ragan MA. Phylogenetic inference based on matrix representation of trees. Mol Phylog Evol. 1992;1:53-58.
        """
        supermatrix = {}
        taxa = set()
        bipartitions_lists = [] # used to store lists of bipartitions, one list per each source tree
        bipartitions_lists_append = bipartitions_lists.append
        _len = len

        # warm-up...
        for source_tree in self.__SOURCE_TREES:
            # full set of taxa
            #taxa |= set([taxon for taxon in source_tree.leaves_identifiers])
            taxa |= set(source_tree.leaves_identifiers)
            bipartitions_list = []
            bipartitions_list_append = bipartitions_list.append

            # computes
            for bipartitions in xfindBipartitions(source_tree):
                bipartitions_list_append(bipartitions)

            bipartitions_lists_append(bipartitions_list)
            self.__number_of_sites += _len(bipartitions_list)

        self.__number_of_species = _len(taxa)

        # computes the supermatrix
        for taxon in taxa:
            sites = []
            sites_append = sites.append

            for bipartitions_list in bipartitions_lists:
                for bipartitions in bipartitions_list:
                    # depending on the source trees, some taxon, of the full set of taxa,
                    # may not be found in any bipartition of one or more source trees
                    if (taxon in bipartitions[0] and taxon in bipartitions[1]) \
                            or ((taxon not in bipartitions[0]) and (taxon not in bipartitions[1])):
                        sites_append(self.__UNKNOWN_SITE)
                    elif (taxon in bipartitions[0]):
                        sites_append(self.__ZERO_OR_A)
                    elif (taxon in bipartitions[1]):
                        sites_append(self.__ONE_OR_C)

            supermatrix[taxon] = sites

        return supermatrix

    def __get_number_of_species(self):
        """
        Returns the number of species that the supermatrix has.

        :return: The number of species that the supermatrix has.
        """
        return self.__number_of_species

    def __get_number_of_sites(self):
        """
        Returns the number of sites that the supermatrix has.

        :return: The number of sites that the supermatrix has.
        """
        return self.__number_of_sites

    def __get_supermatrix(self):
        """
        Returns the supermatrix.

        :return: The supermatrix.
        """
        return self.__SUPERMATRIX

    def __fasta_format(self):
        """
        Returns the supermatrix represented in FASTA format [1].


        References:
        [1] http://en.wikipedia.org/wiki/FASTA_format
        """
        return "\n".join([">{0}\n{1}".format(taxon, "".join(sites))
                          for (taxon, sites) in self.__SUPERMATRIX.items()])

    def __nexus_format(self):
        """
        Returns the supermatrix represented in NEXUS format [1][2].


        References:
        [1] Maddison DR, Swofford DL, Maddison WP. NEXUS: An extensible file format for systematic information.
            Systematic Biology. 1997;46,4:590-621.
        [2] http://en.wikipedia.org/wiki/Nexus_file
        """
        NEXUS = "#NEXUS\nBegin data;\nDimensions ntax={0} nchar={1};\nFormat missing={2};\nMatrix\n{3}\n;\nEnd;\n"
        MATRIX = "\n".join(["{0} {1}".format(taxon, "".join(sites))
                            for (taxon, sites) in self.__SUPERMATRIX.items()])

        return NEXUS.format(self.__number_of_species, self.__number_of_sites, self.__UNKNOWN_SITE, MATRIX)

    def __phylip_format(self):
        """
        Returns the supermatrix represented in PHYLIP (sequential) format [1][2].
        The PHYLIP format can be seen as a simplified version of the NEXUS format [3].


        References:
        [1] Joseph Felsenstein (August 2003). Inferring Phylogenies. Sinauer Associates. ISBN 0-87893-177-5.
        [2] http://en.wikipedia.org/wiki/PHYLIP
        [3] Maddison DR, Swofford DL, Maddison WP. NEXUS: An extensible file format for systematic information.
            Systematic Biology. 1997;46,4:590-621.
        """
        MATRIX = "\n".join(["{0} {1}".format(taxon, "".join(sites))
                            for (taxon, sites) in self.__SUPERMATRIX.items()])

        return "{0} {1}\n{2}".format(self.__number_of_species, self.__number_of_sites, MATRIX)

    def __raw_format(self):
        """
        Returns the supermatrix represented in RAW format.
        The returned string has as many lines as the number of species, where each line has the followin format::
         - starts with the species name;
         - then follows a semicolon (``;``);
         - after follows the counterpart supermatrix row of that species; and
         - ends with a newline (``\n``).
        An example is as follows::

        ``
        A;00001
        B;00011
        C;10000
        D;11000
        ``

        :return: The supermatrix represented in RAW format.
        """
        return "".join(["{};{}\n".format(species, "".join(self.matrix[species]))
                        for species in sorted(self.matrix.keys())])

    def to_string(self, supported_format="RAW"):
        """
        Returns a string representation of the supermatrix.

        :param supported_format: If not "FASTA", "NEXUS", or "PHYLIP" then a RAW representation will be returned
        (see ``__raw_format`` method).
        :return:
        """
        if supported_format == "FASTA":
            return self.fasta_format
        elif supported_format == "NEXUS":
            return self.nexus_format
        elif supported_format == "PHYLIP":
            return self.phylip_format
        else:
            return self.raw_format

    def __str__(self):
        """
        A string representation of the calling instance.
        In other words, a string representation of the supermatrix.

        :return: A string representation of the calling instance.
        In other words, a string representation of the supermatrix.
        """
        #return super.__str__(self)
        return self.to_string()

    def write_supermatrix_to_file(self, filename, supported_format="RAW", validate=False):
        """
        Writes the ``supermatrix`` to a file that has the given ``filename``

        :param filename: The path to the file plus the name of the file (i.e. a complete path)
        :param supported_format: If ``None`` the given ``supermatrix`` will be written to the given ``filename``
        in its pure form; otherwise the given ``supported_format`` will be used.
        :param validate: Allows to bypass the validation of parameter ``filename``.
        In other words, it allows to say:  - In user we trust!  Thus, performance is improved.
        """
        # bypass, or NOT, parameters validation
        if validate:
            if not filename or not isinstance(filename, str):
                raise TypeError("The given \"filename\" has to be a string representing a full path.")

            if system() == "Windows":
                SEP = '\\'
            else:
                SEP = '/'

            (left, sep, right) = filename.rpartition(SEP)

            if left:
                directory = "{0}{1}".format(left, SEP)
            else:
                directory = "./"
            if not right:
                raise ValueError("\"{}\" is not a valid filename.".format(filename))
            # taking care...
            if not os.path.exists(directory):
                try:
                    os.makedirs(directory)
                except:
                    raise ValueError(
                        "\"{}\" does not allow to create the necessary set of directories.".format(filename))

        # write it...
        f = open(filename, 'w')
        f.write(self.to_string(supported_format))
        f.close()

    number_of_species = property(__get_number_of_species, None, None,
                                "Returns the number of species that the supermatrix has.")
    number_of_sites = property(__get_number_of_sites, None, None,
                                "Returns the number of sites that the supermatrix has.")
    matrix = property(__get_supermatrix, None, None,
                                "Returns the supermatrix.")
    fasta_format = property(__fasta_format, None, None,
                                "Returns the supermatrix represented in FASTA format.")
    nexus_format = property(__nexus_format, None, None,
                                "Returns the supermatrix represented in NEXUS format.")
    phylip_format = property(__phylip_format, None, None,
                                "Returns the supermatrix represented in PHYLIP (sequential) format.")
    raw_format = property(__raw_format, None, None,
                                "Returns the supermatrix represented in RAW format.")


if __name__ == "__main__":
    filename = "../datasets/biological/seabirds/kennedy.source_trees_manual"
    source_trees = MatrixRepresentation.source_trees_from_file(filename)
    supported_formats = MatrixRepresentation.supported_formats()    # 0: FASTA, 1: NEXUS, 2: PHYLIP, and 3: RAW
    supermatrix = MatrixRepresentation(source_trees)

    #print(supermatrix.matrix)
    #print(supermatrix.fasta_format)
    #print(supermatrix.phylip_format)
    #print(supermatrix.nexus_format)
    #print(supermatrix.raw_format)
    #print(supermatrix)
    supermatrix.write_supermatrix_to_file("out.mr", supported_format=supported_formats[2])

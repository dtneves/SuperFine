#! /usr/bin/env python

############################################################################
##  test_tree_io.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Tests input/output of trees from files.
"""

import unittest
import datetime
import logging
import tempfile
import os
from optparse import OptionGroup
from optparse import OptionParser
import sys
if sys.version_info > (3, ):
    from io import StringIO
else:
    from cStringIO import StringIO

from dendropy import get_logger
from dendropy.datasets import Dataset
from dendropy import get_logging_level

import dendropy.tests
_LOG = get_logger("TreeParsingAndWriting")

from dendropy import taxa
from dendropy import trees
from dendropy import utils
from dendropy.splits import encode_splits
from dendropy import treedists
from dendropy import datasets

### MODULES THAT WE ARE TESTING ###
from dendropy import nexus
# from dendropy import nexml
### MODULES THAT WE ARE TESTING ###

def iterate_on_trees(tree_files, tf_iterator=nexus.iterate_over_trees):
    "Test (supposedly) memory-economical iteration on trees."
    logging_level = get_logging_level()
    total_tree_files = len(tree_files)
    total_trees = 0
    start_time = datetime.datetime.now()
    if logging_level > logging.INFO:
        minimal_logging = True
    else:
        minimal_logging = False
    _LOG.info("\n*** ITERATOR: <%s> ***" % tf_iterator.__name__)
    for tree_file_idx, tree_file in enumerate(tree_files):
        _LOG.info("   - %s" % os.path.basename(tree_file))
        for tree_idx, tree in enumerate(tf_iterator(file_obj=open(tree_file,'r'))):
            if not minimal_logging:
                _LOG.debug("\n%s" % str(tree))
        total_trees += (tree_idx + 1)
    if not minimal_logging:
        _LOG.debug("\n")
    end_time = datetime.datetime.now()
    _LOG.info("Trees Read: %s" % total_trees)
    _LOG.info("Start time: %s" % start_time)
    _LOG.info("  End time: %s" % end_time)
    run_time = end_time-start_time
    _LOG.info("  Run time: %s" % utils.pretty_print_timedelta(run_time))
    return run_time

def compare_parse_performance(tree_files, methods):
    _LOG.info("\nRunning iterators for (speed) performance comparison ...")
    results = {}
    for method in methods:
        results[method] = iterate_on_trees(tree_files=tree_files, tf_iterator=method)
    _LOG.info("\n---")
    for m1 in methods:
        for m2 in methods[methods.index(m1)+1:]:
            t1 = results[m1]
            t2 = results[m2]
            if t1 >= t2:
                diff = t1 - t2
                diff_sign = "+"
            else:
                diff = t2 - t1
                diff_sign = "-"
            diff_seconds = diff.seconds + float(diff.microseconds)/1000000
            _LOG.info("<%s> vs. <%s> = %s%s seconds " % (m1.__name__, m2.__name__, diff_sign, diff_seconds))

def test_tree_iter_performance(format,
                               heavy=False,
                               wait_to_start=False):
    "Test speed of (supposedly) memory-economical iteration on trees."
    sources = dendropy.tests.data_source_path(format=format, heavy=heavy)
    if wait_to_start:
        raw_input("Hit [ENTER] to begin iterating over trees: ")
    iterate_on_trees(sources)

# def get_anolis_consensus_tree():
#     leaf_nodes = {
#             "Anolis_ahli": 0.2642,
#             "Anolis_aliniger": 0.16,
#             "Anolis_alutaceus": 0.1619,
#             "Anolis_angusticeps": 0.0857,
#             "Anolis_bahorucoensis": 0.2267,
#             "Anolis_barahonae": 0.2115,
#             "Anolis_brevirostris": 0.1801,
#             "Anolis_coelestinus": 0.1932,
#             "Anolis_cristatellus": 0.2144,
#             "Anolis_cuvieri": 0.1687,
#             "Anolis_distichus": 0.1151,
#             "Anolis_equestris": 0.0227,
#             "Anolis_garmani": 0.1068,
#             "Anolis_grahami": 0.0864,
#             "Anolis_insolitus": 0.2439,
#             "Anolis_krugi": 0.1573,
#             "Anolis_lineatopus": 0.1957,
#             "Anolis_loysiana": 0.1836,
#             "Anolis_luteogularis": 0.0306,
#             "Anolis_marcanoi": 0.2359,
#             "Anolis_occultus": 0.4231,
#             "Anolis_olssoni": 0.2569,
#             "Anolis_ophiolepis": 0.0945,
#             "Anolis_paternus": 0.0595,
#             "Anolis_sagrei": 0.0968,
#             "Anolis_strahmi": 0.1978,
#             "Anolis_stratulus": 0.1973,
#             "Anolis_valencienni": 0.1643,
#             "Anolis_vanidicus": 0.206,
#             "Diplolaemus_darwinii": 0.3182,
#     }
#
#     taxa_block = taxa.TaxaBlock()
#     leaf_nodes = []
#     for tax_label in leaf_nodes:
#         taxon = taxa_block.add_taxon(oid="TAXON_"+tax_label, label=tax_label)
#         node = trees.Node(oid="tax_label_"+tax_label, taxon=taxon)
#         node.edge.length = leaf_nodes[tax_label]

def read_newick_tree(tree_filepath):
    "Wrapper to read and return a tree from a single-tree NEWICK file."
    f = open(tree_filepath, 'r')
    tstr = f.read()
    _LOG.info('Reading "%s"' % os.path.basename(tree_filepath))
    _LOG.debug(tstr)
    tree = nexus.parse_newick_string(tstr)
    leaf_nodes = tree.leaf_nodes()
    _LOG.info("%d leaf_nodes on tree: %s" % (len(leaf_nodes), (", ".join([str(n.taxon) for n in leaf_nodes]))))
    return tree

def write_newick_tree(tree, tree_filepath):
    "Wrapper to write a single tree to a NEWICK file."
    nw = nexus.NewickWriter()
    f = open(tree_filepath, 'w')
    tstr = nw.compose_tree(tree)
    _LOG.info('\nWriting "%s"' % os.path.basename(tree_filepath))
    f.write(tstr)

def read_nexus_tree(tree_filepath):
    "Wrapper to read and return a tree from a single-tree NEWICK file."
    f = open(tree_filepath, 'r')
    tstr = f.read()
    _LOG.info('Reading "%s"' % os.path.basename(tree_filepath))
    _LOG.debug(tstr)
    reader = nexus.NexusReader()
    dataset = reader.read_dataset(StringIO(tstr))
    tree = dataset.trees_blocks[0][0]
    leaf_nodes = tree.leaf_nodes()
    _LOG.info("%d leaf_nodes on tree: %s" % (len(leaf_nodes), (", ".join([str(n.taxon) for n in leaf_nodes]))))
    return tree

def write_nexus_tree(tree, tree_filepath):
    "Wrapper to write a single tree to a NEWICK file."
    d = datasets.Dataset()
    taxa_block = tree.infer_taxa_block()
    tree_block = d.add_trees_block(taxa_block=taxa_block)
    tree_block.append(tree)
    nw = nexus.NexusWriter()
    _LOG.info('\nWriting "%s"' % os.path.basename(tree_filepath))
    nw.write_dataset(d, open(tree_filepath, 'w'))

class TreeIOTest(unittest.TestCase):

    def testChangeTranslate(self):
        f = """#NEXUS
Begin taxa ;
    dimensions ntax = 4;
    taxlabels a b c d ;
end;
begin trees;
    translate 
        1 a,
        2 b,
        3 c,
        4 d;
    tree t = (1,2,(3,4));
end;
begin trees;
    translate 
        1 d,
        2 b,
        3 c,
        4 a;
    tree t = (4,2,(3,1));
end;
"""
        d = Dataset()
        d.read(StringIO(f), format="NEXUS")
        t = d.trees_blocks[0][0]
        s = d.trees_blocks[1][0]
        self.assertEqual(t.taxa_block, s.taxa_block)
        encode_splits(s)
        encode_splits(t)
        self.assertEqual(treedists.symmetric_difference(t, s), 0)

    def testNoTranslate(self):
        f = """#NEXUS
Begin taxa ;
    dimensions ntax = 4;
    taxlabels a b c d ;
end;
begin trees;
    tree t = (1,2,(3,4));
    tree s =  (a,b,(d,c));
end;
"""
        d = Dataset()
        d.read(StringIO(f), format="NEXUS")
        t = d.trees_blocks[0][0]
        s = d.trees_blocks[0][1]
        self.assertEqual(t.taxa_block, s.taxa_block)
        encode_splits(s)
        encode_splits(t)
        self.assertEqual(treedists.symmetric_difference(t, s), 0)

    def testParseSpacy(self):
        f = """#NEXUS


BEGIN TAXA;
    DIMENSIONS NTAX=30;
    TAXLABELS
        'Anolis ahli'
        'Anolis garmani'
        'Anolis grahami'
        'Anolis valencienni'
        'Anolis lineatopus'
        'Anolis aliniger'
        'Anolis coelestinus'
        'Anolis bahorucoensis'
        'Anolis equestris'
        'Anolis luteogularis'
        'Anolis occultus'
        'Anolis barahonae'
        'Anolis cuvieri'
        'Anolis insolitus'
        'Anolis olssoni'
        'Anolis brevirostris'
        'Anolis distichus'
        'Anolis cristatellus'
        'Anolis krugi'
        'Anolis stratulus'
        'Anolis alutaceus'
        'Anolis vanidicus'
        'Anolis angusticeps'
        'Anolis paternus'
        'Anolis loysiana'
        'Anolis marcanoi'
        'Anolis strahmi'
        'Diplolaemus darwinii'
        'Anolis ophiolepis'
        'Anolis sagrei'
  ;
END;

BEGIN TREES;
    tree 'con 50 majrule' = [&U] ('Anolis ahli':0.2642130000,((('Anolis garmani':0.1068380000,'Anolis grahami':0.0863670000)1.00:0.069511,'Anolis valencienni':0.1642630000)0.87:0.020752,'Anolis lineatopus':0.1957260000)1.00:0.077682,((((((('Anolis aliniger':0.1600010000,'Anolis coelestinus':0.1932310000)1.00:0.071920,'Anolis bahorucoensis':0.2266880000)0.68:0.023043,('Anolis equestris':0.0227020000,'Anolis luteogularis':0.0306410000)1.00:0.198165,'Anolis occultus':0.4231200000)0.89:0.056277,('Anolis barahonae':0.2114890000,'Anolis cuvieri':0.1686700000)1.00:0.084190,('Anolis insolitus':0.2438820000,'Anolis olssoni':0.2568770000)1.00:0.050618)0.86:0.031679,(('Anolis brevirostris':0.1801300000,'Anolis distichus':0.1151360000)1.00:0.123136,(('Anolis cristatellus':0.2144360000,'Anolis krugi':0.1573300000)0.93:0.036788,'Anolis stratulus':0.1973470000)1.00:0.081037)1.00:0.056582)0.77:0.021826,(('Anolis alutaceus':0.1619060000,'Anolis vanidicus':0.2059960000)1.00:0.118216,(('Anolis angusticeps':0.0857100000,'Anolis paternus':0.0595110000)1.00:0.153413,'Anolis loysiana':0.1836280000)1.00:0.042858)1.00:0.057139,('Anolis marcanoi':0.2359120000,'Anolis strahmi':0.1977660000)1.00:0.141032,'Diplolaemus darwinii':0.6364930000)1.00:0.067869,('Anolis ophiolepis':0.0945010000,'Anolis sagrei':0.0967580000)1.00:0.179398)0.96:0.044895);
END;
"""
        r2 = StringIO(f)
        temp_dataset2 = nexus.NexusReader().read_dataset(file_obj=r2)                                
                


    def setUp(self):
        self.formats = ["newick",
                        "nexus",
#                         "nexml",
                       ]
        self.readers = {"newick": nexus.NewickReader,        
                        "nexus": nexus.NexusReader,
#                         "nexml": nexml.NexmlReader,
                       }        
        self.writers = {"newick": nexus.NewickWriter,        
                        "nexus": nexus.NexusWriter,
#                         "nexml": nexml.NexmlWriter,
                       }
        self.tree_data = ["anolis.mbcon",
                          #"anolis.mcmct",
                          #"primates.mcmct",
                         ]
        self.tree_files = {}
        for format in self.readers:
            self.tree_files[format] = []
            for td in self.tree_data:
                fp = td + ".trees." + format
                _LOG.debug("About to parse %s" % fp)
                self.tree_files[format].append(dendropy.tests.data_source_path(fp))

    def round_trip_tree_file(self,
                             tree_filepath,
                             reader_class,
                             writer_class):
        "Round-trips a treefile."
        reader = reader_class()
        _LOG.info("\nDATA FILE: \"%s\"" % os.path.basename(tree_filepath))
        dataset = reader.read_dataset(file_obj=open(tree_filepath, "r"))      
        for tb_idx, trees_block in enumerate(dataset.trees_blocks):
            for t_idx, tree in enumerate(trees_block):
            
                _LOG.info("*** Tree %d of %d from tree block %d of %d in \"%s\""
                            % (t_idx+1,
                            len(trees_block),
                            tb_idx+1,
                            len(dataset.trees_blocks),
                            os.path.basename(tree_filepath)
                            ))
                
                _LOG.debug("\nORIGINAL TREE >>>\n%s\n<<< ORIGINAL TREE" 
                            % tree.compose_newick()
                              )                
                # write ...
                _LOG.info("(writing out)")
                temp_dataset = datasets.Dataset()
                temp_trees_block = trees.TreesBlock(taxa_block=trees_block.taxa_block)
                temp_trees_block.append(tree)
                temp_dataset.add_trees_block(trees_block=temp_trees_block)
                writer = writer_class()
                result1 = StringIO()
                writer.write_dataset(temp_dataset, result1)
                result1 = result1.getvalue()                               
                _LOG.debug("\nWRITE OUT >>>\n%s\n<<< WRITE OUT" % result1)

                # read back ...
                _LOG.info("(reading back)")           
                r2 = StringIO(result1)
                #r2 = open("/Users/jeet/Documents/Projects/Phyloinformatics/DendroPy/dendropy/dendropy/tests/data/anolis.mbcon.trees.nexus", "r")
                temp_dataset2 = reader.read_dataset(file_obj=r2)                                
                tree2 = temp_dataset2.trees_blocks[0][0]
                result2 = StringIO()
                writer.write_dataset(temp_dataset, result2)
                result2 = result2.getvalue()                
                _LOG.debug("\nREAD IN >>>\n%s\n<<< READ IN" % result2)      
                
                # compare ...
                _LOG.debug("\nREPARSED TREE >>>\n%s\n<<< REPARSED TREE\n" 
                            % tree.compose_newick()
                              )                    
                assert result1 == result2, \
                       "Reparsed tree strings do not match:\n\n" \
                                                       +"FIRST >>>\n%s\n<<< FIRST\n\nSECOND >>>\n%s\n<<< SECOND" % (result1, result2)
                _LOG.info("(reparsed tree string match)")                            

    def testTreeFileParse(self):
        for format in self.formats:
            _LOG.info('\n[Testing %s format parsing: <%s>, <%s>]'  % (format.upper(),
                                                                          self.readers[format].__name__,
                                                                          self.writers[format].__name__))
            for tfile in self.tree_files[format]:
                self.round_trip_tree_file(tfile, self.readers[format], self.writers[format])
    def testStoreEdgeLens(self):
        n = '((((((t4:0.06759,t32:0.06759):0.198252,t39:0.265842):0.135924,((t9:0.244134,(((t23:0.014408,t49:0.014408):0.040121,t16:0.05453):0.156614,t2:0.211144):0.03299):0.013224,t34:0.257358):0.144408):0.112116,(((t45:0.110713,(t47:0.019022,t8:0.019022):0.09169):0.163042,((t1:0.168924,(((((t15:0.012356,t30:0.012356):0.000247,t18:0.012603):0.037913,t22:0.050516):0.076193,(t44:0.071301,t46:0.071301):0.055407):0.037072,(((((t50:0.00000,t29:0.00000):0.01744,t35:0.017441):0.066422,t10:0.083863):0.047231,((t6:0.012709,(t26:0.00805,t40:0.00805):0.004659):0.043941,t11:0.05665):0.074443):0.008316,t31:0.13941):0.024371):0.005144):0.025169,t33:0.194093):0.079662):0.183823,(t48:0.343218,((t41:0.032738,t27:0.032738):0.229887,((t5:0.030394,t43:0.030394):0.204863,((((t14:0.028794,t24:0.028794):0.002007,t3:0.030801):0.181488,t38:0.212289):0.017427,(t17:0.01869,t25:0.01869):0.211027):0.005541):0.027368):0.080592):0.11436):0.056304):0.078832,(t21:0.107754,t13:0.107754):0.48496):0.114273,((t36:0.352531,(((t12:0.042324,t7:0.042324):0.155519,t19:0.197843):0.016322,(t37:0.12614,t28:0.12614):0.088025):0.138366):0.147704,(t42:0.088633,t20:0.088633):0.411601):0.206753)'
        d = Dataset()
        tree = d.trees_from_string(string=n, format="NEWICK")[0]
        for nd in tree.postorder_node_iter():
            if nd is not tree.seed_node:
                if nd.edge.length is None:
                    _LOG.info("%s has edge length of None" % trees.format_node(nd))
                    self.assertTrue(nd.edge.length is not None)

        ts = str(tree)
        d = Dataset()
        tree = d.trees_from_string(string=ts, format="NEWICK")[0]
        for nd in tree.postorder_node_iter():
            if nd is not tree.seed_node:
                if nd.edge.length is None:
                    _LOG.warn("%s has edge length of None" % trees.format_node(nd))
                    self.assertTrue(nd.edge.length is not None)

    def testTaxaWithUnderscoreRead(self):
        rd = dendropy.tests.data_source_path("rana.nex")
        rt = dendropy.tests.data_source_path("rana.tre")
        d = Dataset()
        d.read(open(rd, "rU"), format="NEXUS")
        self.assertEqual(len(d.taxa_blocks[0]), 64)
        d.read_trees(open(rt, "rU"), format="NEXUS")
        self.assertEqual(len(d.taxa_blocks[0]), 64)

def main_local():
    "Main CLI handler."

    parser = OptionParser(add_help_option=True)

    parser.add_option('-p', '--performance',
                      action="store_const",
                      dest="format_type",
                      const="all",
                      help="evaluate NEXUS and NEWICK parsing performance")

    parser.add_option('--NEXUS',
                      action="store_const",
                      dest="format_type",
                      const="NEXUS",
                      help="evaluate NEXUS format parsing performance")

    parser.add_option('--NEWICK',
                      action="store_const",
                      dest="format_type",
                      const="NEWICK",
                      help="evaluate NEWICK format parsing performance")

    parser.add_option('-H', '--heavy',
                        action='store_true',
                        dest='heavy',
                        default=False,
                        help='run heavy (large file) version of performance tests')

    parser.add_option('-w', '--wait',
                        action='store_true',
                        dest='wait',
                        default=False,
                        help='wait for user confirmation before starting runs')

    (opts, args) = parser.parse_args()

    if opts.format_type == "NEXUS":
        filename_filter = "*.trees.nexus"
    elif opts.format_type == "NEWICK":
        filename_filter = "*.trees.newick"
    else:
        filename_filter = "*.trees.*"
    test_tree_iter_performance(filename_filter=filename_filter,
                               heavy=opts.heavy,
                               wait_to_start=opts.wait)


import sys
if __name__ == "__main__":
    if len(sys.argv) == 1:
        unittest.main()
    else:
        main_local()


    #compare_heavy(nexus.iterate_over_trees, "*.newick.tre")
    #compare_heavy(nexus.tree_iter, "*.newick.tre")

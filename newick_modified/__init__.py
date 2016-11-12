'''
A Python module for parsing Newick files.

Copyright (C) 2003-2008, Thomas Mailund <mailund@birc.au.dk>
'''

# convenience inclusion of namespace...
from newick_modified.lexer import LexerError
from newick_modified.parser import *
from newick_modified.tree import parse_tree

if __name__ == '__main__':
    import unittest
    import lexertest
    import parsertest

    test_suite = unittest.TestSuite()
    test_suite.addTest(lexertest.test_suite)
    test_suite.addTest(parsertest.test_suite)

    unittest.TextTestRunner(verbosity=1).run(test_suite)



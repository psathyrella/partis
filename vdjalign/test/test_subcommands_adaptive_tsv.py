import unittest

import pysam

from ..subcommands import adaptive_tsv


def suite():
    s = unittest.TestSuite()
    classes = []
    for cls in classes:
        s.addTests(unittest.makeSuite(cls))
    return s

import unittest

import pysam

from ..subcommands import adaptive_csv


class AnnotateVTestCase(unittest.TestCase):

    def setUp(self):
        self.count_map = {'1': 121, '2': 40}
        self.tid_cysteine_map = {2: 285}
        self.read = pysam.AlignedRead()
        read = self.read

        read.qname = '1'
        read.pos = 206
        read.tid = 2
        read.seq = 'CATCTCCAGAGACAACGCCGAGAACTCACTGTTTCTGCAGATGTACACCCTGAGAGTCGAGGACACGGCTATGTATTACTGTGTGAGACCTCCTGGGGGGTACAACTGGTTCGACCCCTGGGGCCAGGGA'
        read.qual = 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
        read.cigar = [(0, 83), (4, 47)]
        read.tags = [('NM', 7), ('AS', 48)]

        self.assertFalse(self.read.is_secondary)
        self.assertFalse(self.read.is_unmapped)
        self.assertGreater(self.read.aend, 285 + 2)

    def test_annotate(self):
        expected_count = self.count_map[self.read.qname]
        annotated = adaptive_csv.annotate_v_aligned_read(self.tid_cysteine_map, self.count_map, self.read)[0]
        self.assertEqual(expected_count, annotated.opt(adaptive_csv.TAG_COUNT))
        frame = annotated.opt(adaptive_csv.TAG_FRAME)
        self.assertEqual(2, frame)
        cdr3_start = annotated.opt(adaptive_csv.TAG_CDR3_START)
        self.assertEqual(79, cdr3_start)
        self.assertIn(annotated.seq[cdr3_start:cdr3_start+3], ('TGT', 'TGG'))


def suite():
    s = unittest.TestSuite()
    classes = [AnnotateVTestCase]
    for cls in classes:
        s.addTests(unittest.makeSuite(cls))
    return s

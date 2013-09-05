import unittest

from ..ssw import align

class AlignTestCase(unittest.TestCase):

    def test_basic(self):
        res = align('ACCGT', 'TACAGTTTA', False)
        expected = {'mismatches': 1, 'query_begin': 0, 'query_end': 4,
                'cigar_string': '5M', 'sw_score_next_best': 0, 'ref_begin': 1,
                'ref_end_next_best': -1,
                'ref_end': 5, 'sw_score': 6}
        self.assertEqual(expected, res)

    def test_to_pairwise(self):
        res = align('TACATTT', 'ATACAGTTTA')
        self.assertEqual('4M1D3M', res['cigar_string'])
        self.assertEqual('TACA-TTT', res['query_align'])
        self.assertEqual('TACAGTTT', res['ref_align'])

def suite():
    s = unittest.TestSuite()
    s.addTests(unittest.makeSuite(AlignTestCase))
    return s

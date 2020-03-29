###
#   This class contains unit and integration tests for lifesci.bed_utils.
#   The structure of this file is taken from an old blog post here:
#       http://blog.jameskyle.org/2010/10/nose-unit-testing-quick-start/
###

import collections
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import lifesci.bed_utils as bed_utils
import pyllars.logging_utils as logging_utils
import pyllars.math_utils as math_utils

from lifesci.bed_utils import interval_overlap
from lifesci.bed_utils import transcript_overlap

from nose.tools import assert_equal
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from nose.tools import raises

class TestBedUtils(object):

    @classmethod
    def build_a_b_bed_dfs(self):
        # these are values for a set of intervals which are either independent
        # (use ab_info) or grouped into transcripts (use ab_blocks)
        self.a_starts = [10, 30, 40, 70, 90,  110, 130, 155]
        self.a_ends   = [20, 50, 60, 80, 100, 120, 140, 160]
        self.a_info   = [1,  2,  3,  4,  5,   6,   7,   8]
        self.a_blocks = [1,  1,  2,  2,  3,   3,   3,   4]

        self.b_starts = [15, 30, 40, 70, 85,  90,  95,  110, 110, 130, 130, 145]
        self.b_ends   = [20, 50, 60, 80, 100, 100, 115, 120, 120, 140, 140, 155]
        self.b_info   = [1,  2,  3,  4,  5,   6,   7,   8,   9,   10,  11,  12]
        self.b_blocks = [1,  1,  2,  3,  4,   5,   6,   5,   4,   5,   4,   7]

        # now, build a bed6 data frame for a and b
        a_df = pd.DataFrame()
        a_df['start'] = self.a_starts
        a_df['end'] = self.a_ends
        a_df['id'] = self.a_blocks
        a_df['seqname'] = '1'
        a_df['score'] = 0
        a_df['strand'] = '+'
        self.a_df = a_df

        # now, build a bed6 data frame for a and b
        b_df = pd.DataFrame()
        b_df['start'] = self.b_starts
        b_df['end'] = self.b_ends
        b_df['id'] = self.b_blocks
        b_df['seqname'] = '1'
        b_df['score'] = 0
        b_df['strand'] = '+'
        self.b_df = b_df

        a_df['length'] = a_df['end'] -  a_df['start']
        b_df['length'] = b_df['end'] -  b_df['start']

        a_groups = self.a_df.groupby('id')
        a_lengths = a_groups['length'].sum()
        self.a_lengths_map = a_lengths.to_dict()

        b_groups = self.b_df.groupby('id')
        b_lengths = b_groups['length'].sum()
        self.b_lengths_map = b_lengths.to_dict()

        # now, build the bed6 data frames for testing subtraction
        a_starts = np.array([10, 30, 55, 70, 90,  120])
        a_ends   = np.array([20, 40, 60, 80, 101, 130])
        a_ids    = np.array([1,  1,  2,  3,  4,   4])

        b_starts = np.array([40, 65, 100])
        b_ends   = np.array([50, 75, 110])
        b_ids    = np.array([1,  1,  2])

        a_bed = pd.DataFrame()
        a_bed['start'] = a_starts
        a_bed['end'] = a_ends
        a_bed['id'] = a_ids
        a_bed['seqname'] = '1'
        a_bed['score'] = 0
        a_bed['strand'] = '+'


        b_bed = pd.DataFrame()
        b_bed['start'] = b_starts
        b_bed['end'] = b_ends
        b_bed['id'] = b_ids
        b_bed['seqname'] = '1'
        b_bed['score'] = 0
        b_bed['strand'] = '+'

        self.a_diff_df = a_bed
        self.b_diff_df = b_bed

    @classmethod
    def build_complex_bed_entry(self):
        exon_lengths = np.array([20, 20, 70, 70, 70, 70, 30])
        exon_genomic_relative_starts = np.array([0, 40, 90, 190, 290, 390, 490])

        exon_lengths_str = ','.join(str(l) for l in exon_lengths)
        exon_genomic_relative_starts_str = ','.join(str(l) for l in exon_genomic_relative_starts)

        self.bed_entry = {
            'seqname': '1',
            'start': 10,
            'end': 530,
            'id': 'test_me',
            'score': 0,
            'strand': '+',
            'thick_start': 120,
            'thick_end': 450,
            'color': 0,
            'num_exons': 7,
            'exon_lengths': exon_lengths_str,
            'exon_genomic_relative_starts': exon_genomic_relative_starts_str    
        }
        


    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        TestBedUtils.build_a_b_bed_dfs()
        TestBedUtils.build_complex_bed_entry()
        

    @classmethod
    def teardown_class(klass):
        """This method is run once for each class _after_ all tests are run"""

    def setUp(self):
        """This method is run once before _each_ test method is executed"""

    def teardown(self):
        """This method is run once after _each_ test method is executed"""

    def test_get_interval_overlaps(self):
        overlaps = bed_utils.get_interval_overlaps(
            self.a_starts, 
            self.a_ends, 
            self.a_info, 
            self.b_starts, 
            self.b_ends, 
            self.b_info)

        expected_overlaps = [
            interval_overlap(a_info=1, b_info=1, overlap=5),
            interval_overlap(a_info=2, b_info=2, overlap=20),
            interval_overlap(a_info=2, b_info=3, overlap=10),
            interval_overlap(a_info=3, b_info=2, overlap=10),
            interval_overlap(a_info=3, b_info=3, overlap=20),
            interval_overlap(a_info=4, b_info=4, overlap=10),
            interval_overlap(a_info=5, b_info=5, overlap=10),
            interval_overlap(a_info=5, b_info=6, overlap=10),
            interval_overlap(a_info=5, b_info=7, overlap=5),
            interval_overlap(a_info=6, b_info=7, overlap=5),
            interval_overlap(a_info=6, b_info=8, overlap=10),
            interval_overlap(a_info=6, b_info=9, overlap=10),
            interval_overlap(a_info=7, b_info=10, overlap=10),
            interval_overlap(a_info=7, b_info=11, overlap=10)
        ]

        assert_equal(overlaps, expected_overlaps)

    def test_get_transcript_overlaps(self):
        overlaps = bed_utils.get_interval_overlaps(
            self.a_starts, 
            self.a_ends, 
            self.a_blocks, 
            self.b_starts, 
            self.b_ends, 
            self.b_blocks)

        transcript_overlaps = bed_utils.get_transcript_overlaps(overlaps)

        expected_transcript_overlaps = {
            (1, 1): 25,
            (1, 2): 10,
            (2, 1): 10,
            (2, 2): 20,
            (2, 3): 10,
            (3, 4): 30,
            (3, 5): 30,
            (3, 6): 10
        }

        assert_equal(transcript_overlaps, expected_transcript_overlaps)

    def test_get_transcript_overlap_fractions(self):
        overlaps = bed_utils.get_interval_overlaps(
            self.a_starts, 
            self.a_ends, 
            self.a_blocks, 
            self.b_starts, 
            self.b_ends, 
            self.b_blocks)

        transcript_overlaps = bed_utils.get_transcript_overlaps(overlaps)

        transcript_fractions = bed_utils.get_transcript_overlap_fractions(
            transcript_overlaps, 
            self.a_lengths_map, 
            self.b_lengths_map)

        expected_transcript_fractions = [
            transcript_overlap(a_info=1, b_info=2, overlap=10, 
                a_fraction=0.33333333333333331, b_fraction=0.5),
            transcript_overlap(a_info=2, b_info=3, overlap=10, 
                a_fraction=0.33333333333333331, b_fraction=1.0),
            transcript_overlap(a_info=3, b_info=6, overlap=10, 
                a_fraction=0.33333333333333331, b_fraction=0.5),
            transcript_overlap(a_info=2, b_info=2, overlap=20, 
                a_fraction=0.66666666666666663, b_fraction=1.0),
            transcript_overlap(a_info=3, b_info=4, overlap=30, 
                a_fraction=1.0, b_fraction=0.8571428571428571),
            transcript_overlap(a_info=1, b_info=1, overlap=25, 
                a_fraction=0.83333333333333337, b_fraction=1.0),
            transcript_overlap(a_info=2, b_info=1, overlap=10, 
                a_fraction=0.33333333333333331, b_fraction=0.40000000000000002),
            transcript_overlap(a_info=3, b_info=5, overlap=30, 
                a_fraction=1.0, b_fraction=1.0)
        ]

        assert_equal(transcript_fractions, expected_transcript_fractions)

    def test_get_bed_overlaps(self):
        bed_overlaps = bed_utils.get_bed_overlaps(self.a_df, self.b_df, min_a_overlap=0.33, min_b_overlap=0.5)

        expected_bed_overlaps = [
            transcript_overlap(a_info=1, b_info=2, overlap=10, 
                a_fraction=0.33333333333333331, b_fraction=0.5),
            transcript_overlap(a_info=2, b_info=3, overlap=10, 
                a_fraction=0.33333333333333331, b_fraction=1.0),
            transcript_overlap(a_info=3, b_info=6, overlap=10, 
                a_fraction=0.33333333333333331, b_fraction=0.5),
            transcript_overlap(a_info=2, b_info=2, overlap=20, 
                a_fraction=0.66666666666666663, b_fraction=1.0),
            transcript_overlap(a_info=3, b_info=4, overlap=30, 
                a_fraction=1.0, b_fraction=0.8571428571428571),
            transcript_overlap(a_info=1, b_info=1, overlap=25, 
                a_fraction=0.83333333333333337, b_fraction=1.0),
            transcript_overlap(a_info=3, b_info=5, overlap=30, 
                a_fraction=1.0, b_fraction=1.0)
        ]

        assert_equal(bed_overlaps, expected_bed_overlaps)

    def test_retain_thick_only(self):
        thick_only_bed_entry = bed_utils.retain_thick_only(self.bed_entry)

        expected_exon_lengths = np.array([50, 70, 70, 50])
        expected_exon_genomic_relative_starts = np.array([0, 80, 180, 280])

        expected_exon_lengths_str = ','.join(str(l) for l in expected_exon_lengths)
        expected_exon_genomic_relative_starts_str = ','.join(str(l) 
                                    for l in expected_exon_genomic_relative_starts)

        expected_bed_entry = {
            'seqname': '1',
            'start': 120,
            'end': 450,
            'id': 'test_me',
            'score': 0,
            'strand': '+',
            'thick_start': 120,
            'thick_end': 450,
            'color': 0,
            'num_exons': 4,
            'exon_lengths': expected_exon_lengths_str,
            'exon_genomic_relative_starts': expected_exon_genomic_relative_starts_str    
        }

        assert_equal(thick_only_bed_entry, expected_bed_entry)

    def test_retain_before_thick_only(self):
        before_thick_only_bed_entry = bed_utils.retain_before_thick_only(self.bed_entry)

        expected_exon_lengths = np.array([20, 20, 20])
        expected_exon_genomic_relative_starts = np.array([0, 40, 90])

        expected_exon_lengths_str = ','.join(str(l) for l in expected_exon_lengths)
        expected_exon_genomic_relative_starts_str = ','.join(str(l) 
                                    for l in expected_exon_genomic_relative_starts)

        expected_bed_entry = {
            'seqname': '1',
            'start': 10,
            'end': 120,
            'id': 'test_me',
            'score': 0,
            'strand': '+',
            'thick_start': -1,
            'thick_end': -1,
            'color': 0,
            'num_exons': 3,
            'exon_lengths': expected_exon_lengths_str,
            'exon_genomic_relative_starts': expected_exon_genomic_relative_starts_str    
        }

        assert_equal(before_thick_only_bed_entry, expected_bed_entry)

    def test_retain_after_thick_only(self):
        after_thick_only_bed_entry = bed_utils.retain_after_thick_only(self.bed_entry)

        expected_exon_lengths = np.array([20, 30])
        expected_exon_genomic_relative_starts = np.array([0, 50])

        expected_exon_lengths_str = ','.join(str(l) for l in expected_exon_lengths)
        expected_exon_genomic_relative_starts_str = ','.join(str(l) 
                                    for l in expected_exon_genomic_relative_starts)

        expected_bed_entry = {
            'seqname': '1',
            'start': 450,
            'end': 530,
            'id': 'test_me',
            'score': 0,
            'strand': '+',
            'thick_start': -1,
            'thick_end': -1,
            'color': 0,
            'num_exons': 2,
            'exon_lengths': expected_exon_lengths_str,
            'exon_genomic_relative_starts': expected_exon_genomic_relative_starts_str    
        }

        assert_equal(after_thick_only_bed_entry, expected_bed_entry)

    def test_subtract_bed(self):
        # first, check for exact overlaps
        expected_result = {1, 2}
        result = bed_utils.subtract_bed(self.a_diff_df, self.b_diff_df)
        assert_equal(result, expected_result)

        # now, test requiring some more overlap of A
        expected_result = {1, 2, 4}
        min_a_overlap = 0.2
        result = bed_utils.subtract_bed(self.a_diff_df, self.b_diff_df, min_a_overlap=min_a_overlap)
        assert_equal(result, expected_result)

        # and requiring more substantial overlap of B
        expected_result = {1, 2, 3, 4}
        min_b_overlap = 0.3
        result = bed_utils.subtract_bed(self.a_diff_df, self.b_diff_df, min_b_overlap=min_b_overlap)
        assert_equal(result, expected_result)

    def test_get_entries_with_upstream_overlaps(self):
        # construct the data
        a_starts  = np.array([20,  50,  60,  65,  58,  60])
        a_ends    = np.array([30,  60,  70,  75,  63,  78])
        a_strands = np.array(['+', '+', '-', '-', '+', '-'])
        a_ids     = np.array([1,   2,   3,   4,   5,   6])

        b_starts  = np.array([9,   9,   40,  75])
        b_ends    = np.array([20,  20,  55,  80])
        b_strands = np.array(['+', '-', '+', '-'])
        b_ids     = np.array([1,   2,   3,   4])

        bed_a = pd.DataFrame()
        bed_a['start']   = a_starts
        bed_a['end']     = a_ends
        bed_a['strand']  = a_strands
        bed_a['id']      = a_ids
        bed_a['score']   = 0
        bed_a['seqname'] = '1'

        bed_b = pd.DataFrame()
        bed_b['start']   = b_starts
        bed_b['end']     = b_ends
        bed_b['strand']  = b_strands
        bed_b['id']      = b_ids
        bed_b['score']   = 0
        bed_b['seqname'] = '1'

        upstream_window = 5
        allow_overlaps = True

        expected_result = [
            bed_utils.transcript_overlap(a_info=1, b_info=1, overlap=5, 
                a_fraction=1.0, b_fraction=0.45454545454545453),
            bed_utils.transcript_overlap(a_info=2, b_info=3, overlap=5, 
                a_fraction=1.0, b_fraction=0.33333333333333331),
            bed_utils.transcript_overlap(a_info=4, b_info=4, overlap=5, 
                a_fraction=1.0, b_fraction=1.0),
            bed_utils.transcript_overlap(a_info=5, b_info=3, overlap=2, 
                a_fraction=0.40000000000000002, b_fraction=0.13333333333333333),
            bed_utils.transcript_overlap(a_info=6, b_info=4, overlap=2, 
                a_fraction=0.40000000000000002, b_fraction=0.40000000000000002)
        ]

        result = bed_utils.get_entries_with_upstream_overlaps(bed_a, bed_b, 
            upstream_window, allow_overlaps=allow_overlaps)

        assert_equal(expected_result, result)

        upstream_window = 5
        allow_overlaps = False

        expected_result = [
            bed_utils.transcript_overlap(a_info=1, b_info=1, overlap=5, 
                a_fraction=1.0, b_fraction=0.45454545454545453),
            bed_utils.transcript_overlap(a_info=4, b_info=4, overlap=5, 
                a_fraction=1.0, b_fraction=1.0),
            bed_utils.transcript_overlap(a_info=5, b_info=3, overlap=2, 
                a_fraction=0.40000000000000002, b_fraction=0.13333333333333333)
        ]

        result = bed_utils.get_entries_with_upstream_overlaps(bed_a, bed_b, 
            upstream_window, allow_overlaps=allow_overlaps)

        assert_equal(expected_result, result)

        a_starts  = np.array([20,  40,  45,  40])
        a_ends    = np.array([30,  50,  60,  50])
        a_strands = np.array(['+', '+', '-', '+'])
        a_ids     = np.array(['a1','a1','a2','a3'])

        b_starts  = np.array([5,   5,   15,  15,  35,  50,  60,  40])
        b_ends    = np.array([10,  10,  20,  20,  40,  55,  65,  45])
        b_strands = np.array(['+', '-', '+', '-', '+', '-', '-', '+'])
        b_ids     = np.array(['b1','b2','b1','b2','b3','b4','b4','b5'])

        bed_a = pd.DataFrame()
        bed_a['start']   = a_starts
        bed_a['end']     = a_ends
        bed_a['strand']  = a_strands
        bed_a['id']      = a_ids
        bed_a['score']   = 0
        bed_a['seqname'] = '1'

        bed_b = pd.DataFrame()
        bed_b['start']   = b_starts
        bed_b['end']     = b_ends
        bed_b['strand']  = b_strands
        bed_b['id']      = b_ids
        bed_b['score']   = 0
        bed_b['seqname'] = '1'

        exons = pd.concat([bed_a, bed_b])

        # now, create the bed12 data frame
        a1 = {
            'start': 20,
            'end': 50,
            'id': 'a1',
            'strand': '+',
            'num_exons': 2,
            'exon_lengths': "10,10",
            'exon_genomic_relative_starts': "0, 20"
        }

        a2 = {
            'start': 45,
            'end': 60,
            'id': 'a2',
            'strand': '-',
            'num_exons': 1,
            'exon_lengths': "15",
            'exon_genomic_relative_starts': "0"
        }

        a3 = {
            'start': 40,
            'end': 50,
            'id': 'a3',
            'strand': '+',
            'num_exons': 1,
            'exon_lengths': 10,
            'exon_genomic_relative_starts': '0'
        }

        b1 = {
            'start': 5,
            'end': 20,
            'id': 'b1',
            'strand': '+',
            'num_exons': 2,
            'exon_lengths': "5,5",
            'exon_genomic_relative_starts': "0, 10"
        }

        b2 = {
            'start': 5,
            'end': 20,
            'id': 'b2',
            'strand': '-',
            'num_exons': 2,
            'exon_lengths': "5,5",
            'exon_genomic_relative_starts': "0, 10"
        }

        b3 = {
            'start': 35,
            'end': 40,
            'id': 'b3',
            'strand': '+',
            'num_exons': 1,
            'exon_lengths': "5",
            'exon_genomic_relative_starts': "0"
        }

        b4 = {
            'start': 50,
            'end': 65,
            'id': 'b4',
            'strand': '-',
            'num_exons': 2,
            'exon_lengths': "5,5",
            'exon_genomic_relative_starts': "0, 10"
        }

        b5 = {
            'start': 40,
            'end': 45,
            'id': 'b5',
            'strand': '+',
            'num_exons': 1,
            'exon_lengths': "5",
            'exon_genomic_relative_starts': "0"
        }

        bed_a = pd.DataFrame([a1, a2, a3])
        bed_b = pd.DataFrame([b1, b2, b3, b4, b5])

        bed_a['seqname'] = '1'
        bed_a['score'] = 0
        bed_a['thick_start'] = bed_a['start']
        bed_a['thick_end'] = bed_a['end']
        bed_a['color'] = 0


        bed_b['seqname'] = '1'
        bed_b['score'] = 0
        bed_b['thick_start'] = bed_b['start']
        bed_b['thick_end'] = bed_b['end']
        bed_b['color'] = 0

        upstream_window = 5
        allow_overlaps = True

        expected_results = [
            bed_utils.transcript_overlap(a_info='a1', b_info='b1', overlap=5, 
                a_fraction=1.0, b_fraction=0.5),
            bed_utils.transcript_overlap(a_info='a2', b_info='b4', overlap=5, 
                a_fraction=1.0, b_fraction=0.5),
            bed_utils.transcript_overlap(a_info='a3', b_info='b3', overlap=5, 
                a_fraction=1.0, b_fraction=1.0),
        ]

        results = bed_utils.get_entries_with_upstream_overlaps(bed_a, bed_b, 
            upstream_window, allow_overlaps=allow_overlaps, exons=exons)

        assert_equal(results, expected_results)

        upstream_window = 5
        allow_overlaps = False

        expected_results = [
            bed_utils.transcript_overlap(a_info='a1', b_info='b1', overlap=5, 
                a_fraction=1.0, b_fraction=0.5),
            bed_utils.transcript_overlap(a_info='a3', b_info='b3', overlap=5, 
                a_fraction=1.0, b_fraction=1.0)
        ]

        results = bed_utils.get_entries_with_upstream_overlaps(bed_a, bed_b, 
            upstream_window, allow_overlaps=allow_overlaps, exons=exons)

        assert_equal(results, expected_results)





###
#   This class contains unit and integration tests for lifesci.gtf_utils.
#   The structure of this file is taken from an old blog post here:
#       http://blog.jameskyle.org/2010/10/nose-unit-testing-quick-start/
###

import numpy as np
import pandas as pd

import lifesci.bed_utils as bed_utils
import lifesci.gtf_utils as gtf_utils

import pandas.util.testing

class TestGtfUtils(object):

    @classmethod
    def build_exected_gtf_entries(self):
        def get_expected_gtf_feature(start:int, end:int, frame, feature:str, 
                attributes:str=("gene_id my_gene_id; my_numeric_attribute 10; "
                "my_other_attribute \'my other value\'; transcript_id "
                "my_transcript;")):
            exon = {
                "seqname": "chr1",
                "source": "my_source",
                "feature": feature,
                "start": start,
                "end": end,
                "score": 0,
                "frame": frame,
                "strand": "+",
                "attributes": attributes
            }
            
            return pd.Series(exon)

        expected_gtf_exon_starts = np.array([11, 51, 81])
        expected_gtf_exon_ends = np.array([40, 60, 100])
        it = zip(expected_gtf_exon_starts, expected_gtf_exon_ends)

        expected_gtf_exons = [
            get_expected_gtf_feature(start, end, ".", "exon") for start, end in it
        ]

        expected_gtf_cds_starts  = np.array([21, 51, 81])
        expected_gtf_cds_ends = np.array([40, 60, 90])
        #expected_gtf_cds_frame = np.array([0,2,0])
        #it = zip(expected_gtf_cds_starts, expected_gtf_cds_ends, expected_gtf_cds_frame)
        it = zip(expected_gtf_cds_starts, expected_gtf_cds_ends)

        expected_gtf_cds = [
            get_expected_gtf_feature(start, end, ".", "CDS") for start, end in it
        ]

        expected_gtf = expected_gtf_exons + expected_gtf_cds
        expected_gtf = pd.DataFrame(expected_gtf)
        expected_gtf = expected_gtf.sort_values('start')
        expected_gtf = expected_gtf.reset_index(drop=True)
        self.expected_gtf = expected_gtf[gtf_utils.gtf_field_names]

    @classmethod
    def build_bed_entry(self):
        bed_entry = {
            "seqname": "chr1",
            "start": 10,
            "end": 100,
            "id": "my_transcript",
            "score": 0,
            "strand": "+",
            "thick_start": 20,
            "thick_end": 90,
            "color": 0,
            "num_exons": 3,
            "exon_lengths": "30,10,20",
            "exon_genomic_relative_starts": "0,40,70",
            "gene_id": "my_gene_id",
            "my_other_attribute": "my other value",
            "my_numeric_attribute": 10
        }

        self.bed_entry=pd.Series(bed_entry)


    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        TestGtfUtils.build_exected_gtf_entries()
        TestGtfUtils.build_bed_entry()
        

    @classmethod
    def teardown_class(klass):
        """This method is run once for each class _after_ all tests are run"""

    def setUp(self):
        """This method is run once before _each_ test method is executed"""

    def teardown(self):
        """This method is run once after _each_ test method is executed"""

    def test_get_gtf_entries(self):
        
        source = "my_source"
        gtf_entries = gtf_utils.get_gtf_entries(self.bed_entry, source)
        pandas.util.testing.assert_frame_equal(self.expected_gtf, gtf_entries)

import unittest

import numpy as np
from numpy import testing as npt

from newmap.unique_counts import update_upper_search_bound


class TestUpperBoundChange(unittest.TestCase):

    def setUp(self):
        pass

    def test_center_gap(self):
        max_length = 5
        ambiguity_mask = \
            np.array([True, True, False, False, False, True, True])
        num_kmers = ambiguity_mask.size
        sequence_length = num_kmers

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)

        npt.assert_array_equal(upper_bound, [5, 5, 3, 2, 1, 5, 5])

    def test_large_center_gap(self):
        max_length = 4
        ambiguity_mask = \
            np.array([True, False, False, False, False, False, True])
        num_kmers = ambiguity_mask.size
        sequence_length = num_kmers

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)

        npt.assert_array_equal(upper_bound, [4, 4, 4, 3, 2, 1, 4])

    def test_start_gap(self):
        max_length = 4
        ambiguity_mask = \
            np.array([False, False, False, False, False,
                      True, True,
                      False,
                      True])
        num_kmers = ambiguity_mask.size
        sequence_length = num_kmers

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)

        npt.assert_array_equal(upper_bound, [4, 4, 3, 2, 1,
                                             4, 4,
                                             1,
                                             4])

    def test_end_gap(self):
        max_length = 4
        ambiguity_mask = \
            np.array([True, True,
                      False, False,
                      True,
                      False, False, False, False, False])
        num_kmers = ambiguity_mask.size
        sequence_length = num_kmers

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)

        npt.assert_array_equal(upper_bound, [4, 4,
                                             2, 1,
                                             4,
                                             4, 4, 3, 2, 1])

    def test_only_start_gap(self):
        max_length = 4
        ambiguity_mask = \
            np.array([False, False,
                      True, True,])
        num_kmers = ambiguity_mask.size
        sequence_length = num_kmers

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)

        npt.assert_array_equal(upper_bound, [2, 1,
                                             4, 4,])

    def test_only_end_gap(self):
        max_length = 4
        ambiguity_mask = \
            np.array([True, True,
                      False, False,])
        num_kmers = ambiguity_mask.size
        sequence_length = num_kmers

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)

        npt.assert_array_equal(upper_bound, [4, 4,
                                             2, 1,])

    def test_all_gaps(self):
        max_length = 500
        ambiguity_mask = \
            np.array([False, False,
                      True, True,
                      False, False, False,
                      True, True, True,
                      False, False])
        num_kmers = ambiguity_mask.size
        sequence_length = num_kmers

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)

        npt.assert_array_equal(upper_bound, [2, 1,
                                             500, 500,
                                             3, 2, 1,
                                             500, 500, 500,
                                             2, 1])

    def test_end_sequence(self):
        max_length = 50
        ambiguity_mask = \
            np.array([False, False, False, False])
        num_kmers = ambiguity_mask.size
        sequence_length = num_kmers

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)

        npt.assert_array_equal(upper_bound, [4, 3, 2, 1])

    def test_end_sequence_small_lookahead(self):
        max_length = 4
        ambiguity_mask = \
            np.array([False, False, False, False, False, False])
        num_kmers = ambiguity_mask.size
        sequence_length = num_kmers + 1

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)

        npt.assert_array_equal(upper_bound, [4, 4, 4, 4, 3, 2])

    def test_all_ambiguous(self):
        max_length = 50
        ambiguity_mask = \
            np.array([True, True, True, True])
        num_kmers = ambiguity_mask.size
        sequence_length = num_kmers

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)

        npt.assert_array_equal(upper_bound, [50, 50, 50, 50])

    def test_no_ambiguous_max_lookahead(self):
        max_length = 50
        ambiguity_mask = \
            np.array([True, True, True, True])
        num_kmers = ambiguity_mask.size
        lookahead = 50
        sequence_length = num_kmers + (lookahead - 1)

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)

        npt.assert_array_equal(upper_bound, [50, 50, 50, 50])

    def test_no_ambiguous_too_much_lookahead(self):
        max_length = 50
        ambiguity_mask = \
            np.array([True, True, True, True])
        num_kmers = ambiguity_mask.size
        lookahead = 51
        sequence_length = num_kmers + (lookahead - 1)

        upper_bound = np.full(num_kmers, max_length)
        self.assertRaises(AssertionError, update_upper_search_bound,
                          upper_bound, ambiguity_mask, max_length,
                          sequence_length)

    def test_end_gap_with_lookahead_long_length(self):
        max_length = 50
        ambiguity_mask = \
            np.array([True, True, False, False, False, False])
        num_kmers = ambiguity_mask.size
        lookahead = 50
        sequence_length = num_kmers + (lookahead - 1)

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)
        npt.assert_array_equal(upper_bound, [50, 50, 50, 50, 50, 50])

    def test_end_gap_with_short_lookahead_length(self):
        max_length = 50
        ambiguity_mask = \
            np.array([True, True, False, False, False, False])
        num_kmers = ambiguity_mask.size
        lookahead = 3
        sequence_length = num_kmers + lookahead

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)
        npt.assert_array_equal(upper_bound, [50, 50, 7, 6, 5, 4])

    def test_end_gap_with_lookahead_short_max_length(self):
        max_length = 5
        ambiguity_mask = \
            np.array([True, True, False, False, False, False])
        num_kmers = ambiguity_mask.size
        lookahead = 2
        sequence_length = num_kmers + lookahead

        upper_bound = np.full(num_kmers, max_length)
        update_upper_search_bound(upper_bound, ambiguity_mask, max_length,
                                  sequence_length)
        npt.assert_array_equal(upper_bound, [5, 5, 5, 5, 4, 3])


if __name__ == '__main__':
    unittest.main()

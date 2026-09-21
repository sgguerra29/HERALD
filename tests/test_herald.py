# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""Unit tests for HERALD's two decision functions.

Run from the repository root:
    python3 -m unittest discover -s tests -v

Uses only the standard library, so no test dependency beyond tqdm.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from FragAnalyzer import FragAnalyzer  # noqa: E402


def sam_line(name, ref, *tags):
    """Build a split SAM record with the given optional tags."""
    return [name, "0", ref, "1", "60", "50M", "*", "0", "0", "ACGT", "*"] + list(tags)


def block(*pairs):
    """Build the tab-joined rows check_consecutive_references consumes.

    Each pair is (fragment_number, genome). Rows are emitted in fragment
    order, matching how post_process sorts before calling.
    """
    return [f"read_frag_{n}\t{genome}\tAS:i:100" for n, genome in pairs]


class TestGetAlignmentScore(unittest.TestCase):
    """Alignment score lookup by tag prefix rather than column position."""

    def test_reads_tag_at_column_14(self):
        # The layout minimap2 produced in our benchmark runs.
        rec = sam_line("r", "chr1", "NM:i:3", "ms:i:100", "AS:i:250")
        self.assertEqual(FragAnalyzer.get_alignment_score(rec), 250)

    def test_reads_tag_at_first_optional_column(self):
        rec = sam_line("r", "chr1", "AS:i:250", "NM:i:3")
        self.assertEqual(FragAnalyzer.get_alignment_score(rec), 250)

    def test_reads_tag_at_a_late_column(self):
        rec = sam_line("r", "chr1", "NM:i:3", "ms:i:1", "nn:i:0",
                       "tp:A:P", "cm:i:9", "s1:i:80", "AS:i:42")
        self.assertEqual(FragAnalyzer.get_alignment_score(rec), 42)

    def test_missing_tag_returns_sentinel(self):
        rec = sam_line("r", "chr1", "NM:i:3", "MD:Z:50")
        self.assertEqual(FragAnalyzer.get_alignment_score(rec), -1)

    def test_no_optional_tags_returns_sentinel(self):
        self.assertEqual(FragAnalyzer.get_alignment_score(sam_line("r", "chr1")), -1)

    def test_malformed_value_returns_sentinel(self):
        rec = sam_line("r", "chr1", "AS:i:notanumber")
        self.assertEqual(FragAnalyzer.get_alignment_score(rec), -1)

    def test_negative_score_is_preserved(self):
        # A real negative score must stay distinct from the -1 sentinel's meaning, 
        # and must still compare correctly during best-hit selection.
        rec = sam_line("r", "chr1", "AS:i:-8")
        self.assertEqual(FragAnalyzer.get_alignment_score(rec), -8)

    def test_trailing_newline_is_tolerated(self):
        # The last field of an unstripped SAM line carries a newline.
        rec = sam_line("r", "chr1", "NM:i:3", "AS:i:120\n")
        self.assertEqual(FragAnalyzer.get_alignment_score(rec), 120)

    def test_ignores_similar_tags(self):
        # ZS:i: and a tag whose value contains AS must not be mistaken for it.
        rec = sam_line("r", "chr1", "ZS:i:999", "MD:Z:AS:i:5", "AS:i:7")
        self.assertEqual(FragAnalyzer.get_alignment_score(rec), 7)

    def test_mandatory_columns_are_not_scanned(self):
        # A read named "AS:i:5" must not be read as a score.
        rec = sam_line("AS:i:5", "chr1", "NM:i:0")
        self.assertEqual(FragAnalyzer.get_alignment_score(rec), -1)


class TestCheckConsecutiveReferences(unittest.TestCase):
    """Contiguity filter: one background genome, every other genome in one run."""

    def test_aba_flanked_insert_accepted(self):
        # host host donor donor host host — the canonical flanked layout.
        self.assertTrue(FragAnalyzer.check_consecutive_references(
            block((1, "host"), (2, "host"), (3, "donor"),
                  (4, "donor"), (5, "host"), (6, "host"))))

    def test_ab_terminal_insert_accepted(self):
        # donor at the left end, host filling the rest.
        self.assertTrue(FragAnalyzer.check_consecutive_references(
            block((1, "donor"), (2, "donor"), (3, "host"), (4, "host"))))

    def test_ba_terminal_insert_accepted(self):
        self.assertTrue(FragAnalyzer.check_consecutive_references(
            block((1, "host"), (2, "host"), (3, "donor"), (4, "donor"))))

    def test_single_fragment_insert_accepted(self):
        self.assertTrue(FragAnalyzer.check_consecutive_references(
            block((1, "host"), (2, "donor"), (3, "host"))))

    def test_alternating_genomes_rejected(self):
        # No choice of background leaves the other genome in one run.
        self.assertFalse(FragAnalyzer.check_consecutive_references(
            block((1, "host"), (2, "donor"), (3, "host"), (4, "donor"))))

    def test_split_donor_rejected(self):
        # Donor appears in two separate runs; host is already the background.
        self.assertFalse(FragAnalyzer.check_consecutive_references(
            block((1, "host"), (2, "donor"), (3, "host"),
                  (4, "donor"), (5, "host"))))

    def test_two_distinct_inserts_accepted(self):
        # Three genomes: host as background, two individually contiguous inserts from different donors.
        self.assertTrue(FragAnalyzer.check_consecutive_references(
            block((1, "host"), (2, "donorA"), (3, "donorA"),
                  (4, "host"), (5, "donorB"), (6, "host"))))

    def test_two_split_donors_rejected(self):
        self.assertFalse(FragAnalyzer.check_consecutive_references(
            block((1, "host"), (2, "donorA"), (3, "donorB"),
                  (4, "donorA"), (5, "host"))))

    def test_single_genome_accepted(self):
        # Vacuously true. Such reads are filtered earlier by the genome-count check in recollect_fragments, 
        # so this documents the boundary.
        self.assertTrue(FragAnalyzer.check_consecutive_references(
            block((1, "host"), (2, "host"), (3, "host"))))

    def test_gap_in_fragment_numbers_within_a_run(self):
        # Fragment 3 failed to realign, so the donor run is 2 then 4. The filter treats that as non-contiguous and rejects it.
        self.assertFalse(FragAnalyzer.check_consecutive_references(
            block((1, "host"), (2, "donor"), (4, "donor"), (5, "host"))))

    def test_gap_in_the_background_genome_is_tolerated(self):
        # The background genome is exempt from the contiguity requirement, so a missing host fragment does not reject the read.
        self.assertTrue(FragAnalyzer.check_consecutive_references(
            block((1, "host"), (3, "donor"), (4, "donor"), (6, "host"))))


class TestValidateArgs(unittest.TestCase):
    """Default resolution for --min-match-frac and -m."""

    @staticmethod
    def make_args(**overrides):
        import argparse
        defaults = dict(size=10, min_frags=None, min_match_frac=None,
                        max_genomes=2, no_post_process=False)
        defaults.update(overrides)
        return argparse.Namespace(**defaults)

    def setUp(self):
        import launcher
        self.launcher = launcher

    def test_default_is_one_minus_one_over_n(self):
        args = self.launcher.validate_args(self.make_args(size=10))
        self.assertAlmostEqual(args.min_match_frac, 0.9)

    def test_default_tracks_fragment_count(self):
        for size, expected in [(2, 0.5), (4, 0.75), (15, 1 - 1 / 15), (20, 0.95)]:
            with self.subTest(size=size):
                args = self.launcher.validate_args(self.make_args(size=size))
                self.assertAlmostEqual(args.min_match_frac, expected)

    def test_default_always_in_range(self):
        # -s > 1 is enforced separately, so the computed default can never fall outside the range the validator accepts.
        for size in range(2, 200):
            args = self.launcher.validate_args(self.make_args(size=size))
            self.assertTrue(0.0 <= args.min_match_frac <= 1.0)

    def test_explicit_value_is_not_overridden(self):
        args = self.launcher.validate_args(self.make_args(size=10, min_match_frac=0.55))
        self.assertAlmostEqual(args.min_match_frac, 0.55)

    def test_min_frags_defaults_to_size(self):
        args = self.launcher.validate_args(self.make_args(size=12))
        self.assertEqual(args.min_frags, 12)


class TestGetFragmentNumber(unittest.TestCase):

    def test_parses_trailing_number(self):
        self.assertEqual(FragAnalyzer.get_fragment_number(["read_frag_7", "chr1"]), 7)

    def test_tolerates_underscores_in_read_name(self):
        self.assertEqual(
            FragAnalyzer.get_fragment_number(["NORMAL_CP012153.2_read3_frag_12", "chr1"]), 12)


if __name__ == "__main__":
    unittest.main()


class TestEndToEnd(unittest.TestCase):
    """Both stages over the fixture SAMs in tests/data.

    """

    DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")

    def setUp(self):
        import shutil
        import tempfile
        self.tmp = tempfile.mkdtemp()
        shutil.copytree(self.SRC, os.path.join(self.tmp, "src"))
        self.launcher = os.path.join(self.tmp, "src", "launcher.py")

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def run_herald(self, *args):
        import subprocess
        result = subprocess.run(
            [sys.executable, self.launcher] + list(args),
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout

    def test_stage_one_fragments_only_unaligned_reads(self):
        self.run_herald("-i", os.path.join(self.DATA, "stage1_reads.sam"),
                        "-s", "4", "--mmf", "0.9")
        with open(os.path.join(self.tmp, "tmp", "fragments.fasta")) as fh:
            names = [ln[1:].strip() for ln in fh if ln.startswith(">")]

        # The two clean reads align above the threshold and are not fragmented.
        self.assertNotIn("read1_clean_host_frag_1", names)
        self.assertNotIn("read4_clean_host_frag_1", names)
        # The two insert-bearing reads fall below it and are cut into 4 each.
        for read in ("read2_flanked_insert", "read3_terminal_insert"):
            for n in range(1, 5):
                self.assertIn(f"{read}_frag_{n}", names)
        self.assertEqual(len(names), 8)

    def test_stage_one_threshold_controls_admission(self):
        # At a low threshold the partially-aligned reads count as aligned and nothing is fragmented.
        self.run_herald("-i", os.path.join(self.DATA, "stage1_reads.sam"),
                        "-s", "4", "--mmf", "0.3")
        with open(os.path.join(self.tmp, "tmp", "fragments.fasta")) as fh:
            self.assertEqual(fh.read().strip(), "")

    def test_stage_two_post_processed_output(self):
        self.run_herald("-i", os.path.join(self.DATA, "stage2_realigned_fragments.sam"),
                        "-s", "4", "--results", "-o", "fixture")
        out = os.path.join(self.tmp, "Output", "Post_Proc_Fragment_Results_fixture.txt")
        with open(out) as fh:
            text = fh.read()

        # Contiguous inserts survive post-processing.
        self.assertIn("read2_flanked_insert_frag_1", text)
        self.assertIn("read3_terminal_insert_frag_1", text)
        # Alternating genomes do not.
        self.assertNotIn("read5_alternating", text)
        # A read whose fragments all hit one genome is not a candidate.
        self.assertNotIn("read6_single_genome", text)
        # The secondary alignment must not have won on score.
        self.assertIn("read2_flanked_insert_frag_2\tdonor_genome\tAS:i:140", text)
        # The raw intermediate is removed once post-processing succeeds.
        self.assertFalse(os.path.exists(
            os.path.join(self.tmp, "Output", "Fragment_Results_fixture.txt")))

    def test_stage_two_without_post_processing(self):
        self.run_herald("-i", os.path.join(self.DATA, "stage2_realigned_fragments.sam"),
                        "-s", "4", "--results", "-o", "raw", "--no-post-process")
        out = os.path.join(self.tmp, "Output", "Fragment_Results_raw.txt")
        with open(out) as fh:
            text = fh.read()
        # Without the filter, the alternating read is retained as a candidate.
        self.assertIn("read5_alternating", text)

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""HERALD command-line entry point.

Two passes over the same tool:
  1. Slicing: reads an aligned .sam, writes tmp/fragments.fasta.
  2. Result processing (--results): reads the realigned .sam, writes Output/.
"""

import argparse
import sys
from FragAnalyzer import FragAnalyzer
from Fragmenter import Fragmenter


def setup_parser():
    """Define and parse command-line arguments."""
    parser = argparse.ArgumentParser(description='HERALD: horizontal gene transfer detection')

    parser.add_argument("--results", help="Result processing mode", action='store_true')
    parser.add_argument("-o", "--output", help="Output file name", action='store')
    parser.add_argument("-i", "--input",
                        help="Input file. Either initial sam file or results to process (with --results)",
                        required=True)

    group = parser.add_argument_group('Fragmenting Options')
    group.add_argument("-s", "--size", help="Fragments per sequence [default: 10]", type=int,
                       default=10, action="store")
    group.add_argument("-m", "--min-frags", help="Minimum matched fragments for candidacy [default: Equal to -s]",
                       type=int, default=None, action="store")
    group.add_argument("--min-match-frac", "--mmf", dest="min_match_frac",
                       help="Minimum fraction of the read that must match (M bases / read length) to count as aligned [default: 1 - 1/s, i.e. 0.9 at the default 10 fragments]",
                       type=float, default=None, action="store")
    group.add_argument("--max-genomes", help="Maximum number of different genomes fragments can align to for a read to be considered an HGT candidate [default: 2]", type=int, default=2)
    group.add_argument("--no-post-process", help="Skip the consecutive-genome post-processing filter. Post-processing assumes every candidate read has exactly -s matched fragments, so this is useful if you set -m below -s.", action="store_true")

    args = parser.parse_args()
    return args


def validate_args(args):
    """Range-check arguments and apply the -m and --mmf defaults.

    Exits with status 1 on invalid input.
    """
    if args.size <= 1:
        print("Error: -s/--size: Value must be greater than 1.")
        sys.exit(1)

    if args.min_frags is None:
        args.min_frags = args.size

    if args.min_match_frac is None:
        # Match the admission gate to the fragmenter: a read is admitted once the unmatched portion 
        # approaches one fragment, which is the smallest amount of donor sequence fragmentation can resolve.
        args.min_match_frac = 1 - 1 / args.size

    if args.min_frags > args.size:
        print("Error: -m/--min-frags: Value must be less than or equal to -s/--size.")
        sys.exit(1)

    if args.max_genomes < 2:
        print("Error: --max-genomes: Value must be at least 2.")
        sys.exit(1)

    if not 0.0 <= args.min_match_frac <= 1.0:
        print("Error: --min-match-frac: Value must be between 0.0 and 1.0.")
        sys.exit(1)

    if args.min_frags < args.size and not args.no_post_process:
        print("Warning: post-processing assumes every candidate read has exactly -s matched fragments. "
              "With -m set below -s, post-processing results may be unreliable. "
              "Consider passing --no-post-process.")

    return args


def fragment_launcher(args):
    """Dispatch to the slicing or result-processing stage."""
    frag_size = args.size
    min_frags = args.min_frags
    result_mode = args.results is True
    input_file = args.input
    output_file = args.output

    if result_mode:
        print("Launching result processor")
        frag_analyzer = FragAnalyzer(input_file, frag_size, min_frags, output_file, args.max_genomes,
                                     not args.no_post_process)
        frag_analyzer.recollect_fragments()
    else:
        print("Launching fragment slicer")
        print(f"Minimum matched fraction: {args.min_match_frac:g}")
        frag = Fragmenter(input_file, frag_size, min_frags, args.min_match_frac)
        frag.find_unaligned_seqs()
        frag.fragment_seq()


def main():
    args = validate_args(setup_parser())
    try:
        fragment_launcher(args)
    except (ValueError, FileNotFoundError) as err:
        # Expected user-input problems: report plainly rather than as a traceback.
        print(f"Error: {err}")
        sys.exit(1)


if __name__ == "__main__":
    main()

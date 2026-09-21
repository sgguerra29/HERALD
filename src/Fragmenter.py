# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""pt 1: find unaligned reads in a SAM file and slice them into fragments.

Output is a FASTA of fragments written to tmp/fragments.fasta. Intended to be
realigned with an external aligner before being passed to FragAnalyzer.
"""

import os
import re
import tqdm

# Repository root, resolved from this file rather than the working directory.
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class Fragmenter:
    """Splits poorly-aligned reads into equal-length fragments for realignment.

    Args:
        input_file: Path to an aligned .sam file.
        fragments_per_seq: Number of fragments to cut each read into.
        min_matched_frags: Accepted for interface symmetry with FragAnalyzer;
            fragment generation does not use it.
        min_match_frac: Fraction of a read's bases that must be CIGAR M
            operations for the read to count as aligned. Reads below this are
            treated as unaligned and fragmented.
    """

    def __init__(self, input_file, fragments_per_seq, min_matched_frags, min_match_frac=0.55):
        self.min_match_frac = min_match_frac
        self.sequence_dict = {}
        self.unaligned_seq_list = []
        self.parent_dir = PROJECT_ROOT
        self.temp_dir = os.path.join(PROJECT_ROOT, "tmp")

        self.sam_file = input_file
        self.fragments_per_seq = fragments_per_seq

    @staticmethod
    def get_matched_bases(cigar):
        """Sum all M operations in a CIGAR string."""
        return sum(int(n) for n in re.findall(r'(\d+)M', cigar))

    def find_unaligned_seqs(self):
        """Populate self.unaligned_seq_list with (name, sequence) for unaligned reads.

        A read may appear on several SAM lines. Alignment status is taken only
        from the primary line; the longest observed sequence is retained.
        Consumes self.sequence_dict.

        Raises:
            ValueError: If the input file is not a .sam file.
        """
        print("Processing SAM input:")
        if not self.sam_file.endswith(".sam"):
            raise ValueError(
                f"Input must be a .sam file, got: {self.sam_file}. "
                "Convert BAM/CRAM input with 'samtools view -h' first."
            )

        with tqdm.tqdm(total=os.path.getsize(self.sam_file)) as pbar:
            with open(self.sam_file, "r") as infile:
                for line in infile:
                    pbar.update(len(line))
                    if line[0] != '@':
                        line = line.split('\t')
                        sequence_name = line[0]
                        flag = int(line[1])
                        is_secondary = flag & 256
                        is_supplementary = flag & 2048
                        is_primary = not is_secondary and not is_supplementary
                        cigar = line[5]
                        read_len = len(line[9]) if line[9] != '*' else 0
                        matched = Fragmenter.get_matched_bases(cigar) if cigar != '*' else 0
                        frac_matched = matched / read_len if read_len > 0 else 0
                        is_aligned = (line[2] != '*') and (flag & 4 == 0) and (frac_matched >= self.min_match_frac)
                        sequence = line[9]

                        if sequence_name in self.sequence_dict:
                            existing_seq = self.sequence_dict[sequence_name][1]

                            # Always keep the longest sequence seen
                            if sequence and sequence != '*':
                                new_seq = sequence if (existing_seq is None or len(sequence) > len(existing_seq)) else existing_seq
                            else:
                                new_seq = existing_seq

                            # Only update alignment status from primary alignments
                            if is_primary:
                                updated_aligned = is_aligned
                            else:
                                updated_aligned = self.sequence_dict[sequence_name][0]

                            self.sequence_dict[sequence_name] = [updated_aligned, new_seq, sequence_name]

                            # If sequence is aligned we don't need to store its sequence.
                            if self.sequence_dict[sequence_name][0]:
                                self.sequence_dict[sequence_name][1] = None

                        else:
                            seq_to_store = sequence if (sequence and sequence != '*') else None
                            aligned_to_store = is_aligned if is_primary else False
                            # An aligned read never gets fragmented, so don't hold its sequence.
                            if aligned_to_store:
                                seq_to_store = None
                            self.sequence_dict[sequence_name] = [aligned_to_store, seq_to_store, sequence_name]

        for sequence in self.sequence_dict.values():
            is_aligned = sequence[0]
            seq_name = sequence[2]
            sequence = sequence[1]
            if not is_aligned:
                self.unaligned_seq_list.append((seq_name, sequence))

        del self.sequence_dict
        print("Completed SAM Processing")

    def fragment_seq(self):
        """Write each unaligned read as fragments_per_seq FASTA records.

        Fragments are named "<read>_frag_<n>", 1-indexed; FragAnalyzer parses that suffix to 
        recover fragment order. Fragment length is the floor of read length over fragment count,
        with the final fragment absorbing the remainder. Consumes self.unaligned_seq_list.
        """
        print("Slicing unaligned sequences:")
        os.makedirs(self.temp_dir, exist_ok=True)
        skipped = 0
        out_path = os.path.join(self.temp_dir, "fragments.fasta")
        with open(out_path, "w") as outfile:
            for sequence_name, sequence in tqdm.tqdm(self.unaligned_seq_list):
                if not sequence:
                    skipped += 1
                    continue
                frag_len = len(sequence) // self.fragments_per_seq

                for i in range(self.fragments_per_seq):
                    start_idx = i * frag_len
                    end_idx = start_idx + frag_len if i != self.fragments_per_seq - 1 else len(sequence)

                    label = i + 1
                    outfile.write(f">{sequence_name}_frag_{label}\n")
                    outfile.write(f"{sequence[start_idx:end_idx]}\n")

        if skipped > 0:
            print(f"Warning: skipped {skipped} reads with missing sequence data.")
        del self.unaligned_seq_list
        print("Completed Fragment Sequence Processing")
        print(f"Fasta for alignment generated at: {out_path}")

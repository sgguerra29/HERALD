# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""pt 2: score realigned fragments and report horizontal-transfer candidates.

Consumes the SAM produced by realigning tmp/fragments.fasta and writes results to Output/.
"""

import os
import tqdm
import datetime
from SequenceData import SequenceData

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class FragAnalyzer:
    """group realigned fragments by read and filters for HGT candidates.

    Args:
        input_file: Path to the .sam file of realigned fragments.
        frag_size: Fragments per read used at slicing time. Recorded for
            reference; candidate filtering uses min_frags.
        min_frags: Minimum aligned fragments a read needs to be a candidate.
        output_file: Suffix for the result filename. if falsy, a timestamp. 
        max_genomes: Upper bound on distinct genomes a candidate read may
            align to.
        post_process_enabled: Apply the contiguity filter and keep only the
            filtered output.
    """

    def __init__(self, input_file, frag_size, min_frags, output_file, max_genomes=2, post_process_enabled=True):
        identifier = "".join(str(datetime.datetime.now())[:-5].split(":"))
        self.parent_dir = PROJECT_ROOT
        self.temp_dir = os.path.join(PROJECT_ROOT, "tmp")
        self.input_file = input_file
        self.max_genomes = max_genomes
        self.post_process_enabled = post_process_enabled

        if output_file:
            self.output_file = f"Fragment_Results_{output_file}"
        else:
            self.output_file = f"Fragment_Results_{identifier}"
        self.min_matched_frags = min_frags
        self.frag_size = frag_size
        self.realigned_seqs = {}

    @staticmethod
    def get_fragment_number(record):
        """Recover the 1-indexed fragment number from a split SAM record."""
        return int(record[0].split('_')[-1])

    @staticmethod
    def get_alignment_score(fields):
        """Return the AS:i: alignment score from a split SAM record.

        Returns -1 when the tag is absent or unparseable, which ranks the record below any scored
        alignment during best-hit selection.
        """
        for field in fields[11:]:
            if field.startswith("AS:i:"):
                try:
                    return int(field[5:])
                except ValueError:
                    return -1
        return -1

    @staticmethod
    def check_consecutive_references(dataset):
        """
        looks for way to designate one  genome as the "background." can be split across separate
        flanking runs (e.g. a host genome flanking a foreign insert) — such that every OTHER 
        genome's fragments form a single, uninterrupted run of fragment numbers (a clean insert).

        with 2 genomes it reduces to "host flanks insert"; with 3+ genomes it allows multiple
        distinct, individually-contiguous inserts (from different donor genomes) within a single read.
        """
        references = [line.split("\t")[1] for line in dataset]
        frag_numbers = [FragAnalyzer.get_fragment_number(line.split("\t")) for line in dataset]

        # Group fragment numbers by genome, preserving fragment-number order
        # (dataset is pre-sorted by fragment number before this is called).
        genome_frag_numbers = {}
        for ref, frag_num in zip(references, frag_numbers):
            genome_frag_numbers.setdefault(ref, []).append(frag_num)

        def is_consecutive(numbers):
            return numbers == list(range(numbers[0], numbers[-1] + 1))

        genomes = list(genome_frag_numbers.keys())

        # Try each genome in turn as the "background" genome.
        for background in genomes:
            if all(is_consecutive(genome_frag_numbers[g]) for g in genomes if g != background):
                return True

        return False

    def recollect_fragments(self):
        """Read the realigned SAM, select candidates, and write results.

        Each fragment is kept once, at its best alignment score. A read is reported when it has
        at least min_matched_frags aligned fragments spanning between 2 and max_genomes distinct genomes.
        Reads are written as blank-line-separated blocks of "fragment<TAB>genome<TAB>AS:i:score" rows.
        """
        print("Analyzing realigned fragments")
        with tqdm.tqdm(total=os.path.getsize(f"{self.input_file}")) as pbar:
            with open(f"{self.input_file}", "r") as infile:
                for line in infile:
                    pbar.update(len(line))

                    if line[0] != '@':
                        line = line.rstrip("\n").split('\t')
                        flag = int(line[1])
                        is_secondary = flag & 256
                        is_supplementary = flag & 2048
                        if is_secondary or is_supplementary:
                            continue
                        seq_name = line[0]
                        genome = line[2]
                        if genome == '*':
                            continue
                        alignment_score = FragAnalyzer.get_alignment_score(line)

                        if seq_name in self.realigned_seqs.keys():
                            old_score = self.realigned_seqs[seq_name][1]
                            if alignment_score > old_score:
                                self.realigned_seqs[seq_name] = [genome, alignment_score]
                                continue
                            else:
                                continue
                        else:
                            self.realigned_seqs[seq_name] = [genome, alignment_score]

        seq_dict = {}
        for frag_name, entry in self.realigned_seqs.items():
            seq_name = frag_name.split("_frag_")[0]
            genome = self.realigned_seqs[frag_name][0]
            alignment_score = self.realigned_seqs[frag_name][1]
            if seq_name not in seq_dict.keys():
                seq_dict[seq_name] = SequenceData(seq_name)

            seq_dict[seq_name].genomes.add(genome)
            seq_dict[seq_name].fragments.add(frag_name)
            seq_dict[seq_name].pairs.append([frag_name, genome, alignment_score])

        output = ""
        for entry in seq_dict.values():
            frag_count = len(entry.fragments)
            genome_count = len(entry.genomes)

            # Not enough aligned fragments or not enough species -> no results
            if genome_count < 2 or genome_count > self.max_genomes or frag_count < self.min_matched_frags:
                continue

            for pair in sorted(entry.pairs, key=FragAnalyzer.get_fragment_number):
                output_str = f"{pair[0]}\t{pair[1]}\tAS:i:{str(pair[2])}\n"
                output += output_str
            output += "\n"

        output_dir = os.path.join(self.parent_dir, "Output")
        os.makedirs(output_dir, exist_ok=True)
        raw_path = os.path.join(output_dir, f"{self.output_file}.txt")
        with open(raw_path, "w") as f:
            f.write(output)

        if self.post_process_enabled:
            self.post_process()
            # The raw file is just an intermediate once post-processing has
            # succeeded — only the filtered file is meant to be the deliverable.
            os.remove(raw_path)
            print(f"Results generated at: {os.path.join(output_dir, f'Post_Proc_{self.output_file}.txt')}")
        else:
            print(f"Results generated at: {raw_path}")

    def post_process(self):
        """Filter the raw result file down to reads with contiguous inserts.

        Re-reads the blocks written by recollect_fragments and keeps only those
        passing check_consecutive_references, writing them to Output/Post_Proc_<output_file>.txt.
        """

        def read_file():
            valid_res = []
            entry = []

            def finalize_entry():
                if not entry:
                    return
                sorted_entry = sorted(entry, key=FragAnalyzer.get_fragment_number)
                str_entry = ["\t".join(a) for a in sorted_entry]
                if FragAnalyzer.check_consecutive_references(str_entry):
                    valid_res.append(sorted_entry)

            with open(os.path.join(self.parent_dir, "Output", f"{self.output_file}.txt"), "r") as f:
                for line in f:
                    if line.strip() == "":
                        # Blank line marks the end of the current entry, regardless of how many fragments it contains.
                        finalize_entry()
                        entry = []
                        continue
                    entry.append(line.split("\t"))

            # Handle a final entry that isn't followed by a trailing blank line.
            finalize_entry()
            return valid_res

        results = read_file()
        prepped_results = []
        for read in results:
            rows = []
            for row in read:
                row = "\t".join(row)
                rows.append(row)
            prepped_results.append("".join(rows))

        out_path = os.path.join(self.parent_dir, "Output", f"Post_Proc_{self.output_file}.txt")
        with open(out_path, "w") as out:
            out.write("\n".join(prepped_results))

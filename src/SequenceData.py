# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


class SequenceData:
    """Per-read accumulator for realigned fragments.

    Attributes:
        sequence: Read name (the portion of a fragment name before "_frag_").
        genomes: Distinct reference names the read's fragments aligned to.
        fragments: Distinct fragment names that aligned.
        pairs: [fragment_name, genome, alignment_score] rows, unordered.
    """

    def __init__(self, sequence):
        self.sequence = sequence
        self.genomes = set()
        self.fragments = set()
        self.pairs = []

    def __repr__(self):
        return f"({self.sequence}, {self.genomes}, {self.fragments}, {self.pairs})"

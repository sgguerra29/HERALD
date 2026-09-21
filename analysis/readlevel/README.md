# Read-level module analysis (modules 1 and 2)

`results_single_frag.py` counts, per read and broken down by read length, how
many fragments aligned to a genome other than the one the read was simulated
from, bucketed as 0, 1, or more than 1 foreign fragment.

Single-fragment correction described in Methods 2.7: candidates supported by exactly one foreign fragment were treated as spurious. Both raw nd corrected counts are reported. HERALD does not perform this correction
itself, automatically or parametrically.

The counts were first done by hand. This script was written once the files grew too large to count manually, and the manual results were then checked against it and agreed.
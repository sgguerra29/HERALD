<a id="readme-top"></a>

[![Python][python-shield]][Python-url]
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![License: MPL 2.0][license-shield]][license-url]

<br />
<div align="center">
<h3 align="center">HERALD</h3>

  <p align="center">
    Fragment-based detection of candidate horizontal gene transfer events in aligned sequencing reads
    <br />
    <a href="https://github.com/sgguerra29/HERALD/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/sgguerra29/HERALD/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>

<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a></li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#parameters">Parameters</a></li>
    <li><a href="#output-format">Output Format</a></li>
    <li><a href="#requirements-on-the-aligner">Requirements on the Aligner</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

## About The Project

HERALD identifies reads that may span a horizontal gene transfer junction.

Reads that fail to align cleanly to the reference are cut into equal-length
fragments and realigned. A read whose fragments land on more than one genome is
a candidate: the pattern expected from a junction is a run of fragments matching
one genome flanked by fragments matching another. An optional post-processing
filter keeps only reads whose fragment-to-genome assignments form contiguous
runs, which removes candidates whose genome assignments alternate in ways no
single insertion could produce.

Detection runs in two passes, with an alignment step of your choice in between.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Repository Layout

```
src/            HERALD itself. Depends only on tqdm.
run.sh          Wrapper for src/launcher.py.
tests/          Unit and end-to-end tests, with small SAM fixtures.
benchmarks/     Simulation scripts used to generate the manuscript test data.
analysis/       Scored results, scoring code, and the R figure scripts.
  scoring/      HERALD output + ground truth -> the scored TSVs (module 4)
  pacbio/       Candidate classification and cross-checks (module 3)
  readlevel/    Single-fragment correction (modules 1 and 2)
```

`benchmarks/` and `analysis/` are provided for transparency and are not needed
to run HERALD. Neither is imported by anything in `src/`.

## Getting Started

### Prerequisites

* [Python >= 3.7](https://www.python.org/downloads/)
* A short-read or long-read aligner. Results reported with this tool were
  produced using [minimap2](https://github.com/lh3/minimap2) 2.30-r1287.

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/sgguerra29/HERALD.git
   ```
2. Install dependencies
   ```sh
   pip3 install tqdm
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Usage

1. From the HERALD directory, slice unaligned reads into 10 fragments each:
   ```sh
   ./run.sh --size 10 --min-frags 10 --input <input_file.sam>
   ```
2. Fragments for realignment are written to `tmp/fragments.fasta`.
3. Realign `tmp/fragments.fasta` with your aligner of choice. Use settings that
   report a single primary alignment per fragment and emit an `AS:i:` tag; see
   [Requirements on the Aligner](#requirements-on-the-aligner).
4. Process the realigned output:
   ```sh
   ./run.sh --size 10 --min-frags 10 --input <re-aligned.sam> --results
   ```
5. Results are written to `Output/Post_Proc_Fragment_Results_<id>.txt`, or to
   `Output/Fragment_Results_<id>.txt` when `--no-post-process` is passed.

`run.sh` can be invoked from any working directory; `tmp/` and `Output/` are
always resolved relative to the HERALD source tree and created if missing.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Parameters

| Argument | Function | Default |
|---|---|---|
| `-h, --help` | Displays arguments and usage | -- |
| `-i, --input` | Input aligned `.sam` file. Applies to both passes. **Mandatory.** | -- |
| `--results` | Run the result processor rather than the fragment slicer. Pass alongside the same parameters used in pass 1. | off |
| `-o, --output` | Output file identifier. Applies to `--results`. | timestamp |
| `-s, --size` | Fragments to create per read. Must be greater than 1. | 10 |
| `-m, --min-frags` | Fragments that must realign for a read to qualify as a candidate. Must be less than or equal to `-s`. Can be changed at result-processing time. | equal to `-s` |
| `--min-match-frac`, `--mmf` | Fraction of a read's bases that must be CIGAR `M` operations for the read to count as aligned. Reads below this are fragmented. Must be between 0.0 and 1.0. | `1 - 1/s` |
| `--max-genomes` | Maximum distinct genomes a candidate read's fragments may align to. Must be at least 2. | 2 |
| `--no-post-process` | Skip the contiguity filter and keep the unfiltered result file. Post-processing assumes every candidate read has exactly `-s` matched fragments, so this is recommended whenever `-m` is set below `-s`. | off |

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Choosing `--mmf`

The default is `1 - 1/s`, which is 0.9 at the default 10 fragments. This matches
the admission gate to the fragmenter: a read is admitted once the portion
failing to match its reference approaches the length of a single fragment, which
is the smallest amount of donor sequence that fragmentation can resolve.

Recall is highly sensitive to this parameter. Because the gate keys on how well
a clean read matches its reference, the right value depends on read
accuracy and reference divergence. As such, set it from the matched-fraction distribution
of your own reads rather than relying on our default: a threshold below `1 - 1/s`
excludes reads that fragmentation could have resolved, while a threshold
approaching the matched fraction of a clean read admits reads indiscriminately.
Lower it for higher-error data. The effective value is printed at the start of
every slicing run.

## Output Format

Results are tab-separated, one row per realigned fragment:

```
<read>_frag_<n>    <genome>    AS:i:<score>
```

Rows are grouped by read and ordered by fragment number. In the unfiltered
output, blank lines separate reads.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Requirements on the Aligner

HERALD reads the realignment SAM directly, so a few aligner settings matter:

* **Primary alignments only.** Secondary (flag `0x100`) and supplementary
  (flag `0x800`) records are skipped. Configuring the aligner to report one
  alignment per fragment is recommended.
* **`AS:i:` tags.** Fragments are deduplicated by best alignment score. The tag
  is located by prefix, so its column position does not matter, but a fragment
  with no `AS:i:` tag is scored `-1` and will lose to any scored alignment.
* **Fragment names must be preserved.** The `_frag_<n>` suffix carries the
  fragment ordering that post-processing depends on. Aligners that truncate
  read names at whitespace are fine; ones that rewrite names are not.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Tests

```sh
python3 -m unittest discover -s tests -v
```

Standard library only; no test dependency beyond `tqdm`. The suite covers the
alignment-score and contiguity functions directly, and runs both stages end to
end over the small SAM fixtures in `tests/data/`. Those fixtures are also a
usable worked example if you want to see HERALD's input and output formats
without supplying your own data.

## License

Distributed under the Mozilla Public License Version 2.0. See `LICENSE.txt` for
more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Contact

##### Maintainer:
Sophia Guerra - [GitHub](https://github.com/sgguerra29)

##### Original author:
Aaron Saporito - [LinkedIn](https://www.linkedin.com/in/aaron-saporito) - [GitHub](https://github.com/aasaporito)

Project Link: [https://github.com/sgguerra29/HERALD](https://github.com/sgguerra29/HERALD)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Acknowledgments
* <a href="https://www.conncoll.edu/directories/faculty-profiles/stephen-douglass/">Dr. Stephen Douglass</a>

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
[contributors-shield]: https://img.shields.io/github/contributors/sgguerra29/HERALD.svg?style=for-the-badge
[contributors-url]: https://github.com/sgguerra29/HERALD/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/sgguerra29/HERALD.svg?style=for-the-badge
[forks-url]: https://github.com/sgguerra29/HERALD/network/members
[stars-shield]: https://img.shields.io/github/stars/sgguerra29/HERALD.svg?style=for-the-badge
[stars-url]: https://github.com/sgguerra29/HERALD/stargazers
[issues-shield]: https://img.shields.io/github/issues/sgguerra29/HERALD.svg?style=for-the-badge
[issues-url]: https://github.com/sgguerra29/HERALD/issues
[license-shield]: https://img.shields.io/badge/License-MPL_2.0-brightgreen.svg?style=for-the-badge
[license-url]: https://github.com/sgguerra29/HERALD/blob/master/LICENSE.txt
[Python-url]: https://python.org
[python-shield]: https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white

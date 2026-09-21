#!/usr/bin/env python3
"""
build_msa1003_reference.py

Build a reference for ATCC MSA-1003 from genomes downloaded from NCBI.

Reference set follows Supplementary Table 1 of:
  Hon T. et al. "Highly accurate long-read HiFi sequencing data for five
  complex genomes." Sci Data 7, 399 (2020). doi:10.1038/s41597-020-00743-4

------------
1. Scans a directory of downloaded FASTA files (any layout: an unzipped NCBI Datasets package, 
loose .fna/.fasta files, gzipped or not).
2. Identifies every sequence record by the accession in its header and assigns it to one of the 20 organisms.
3. Concatenates all replicons of an organism (chromosome + plasmids + contigs) 
into one record named for the organism, joined by runs of N so nothing aligns across a junction.
4. Writes a contig map, and a QC report comparing assembled size against the published genome size.

Usage
-----
  # print the accession list to feed to NCBI, then stop
  python3 build_msa1003_reference.py --list-accessions

  # build the reference from downloaded files
  python3 build_msa1003_reference.py --input downloads/ --output msa1003_ref.fasta

Options
-------
  --spacer N        length of the N run between replicons [200]
  --tolerance F     fractional size mismatch that triggers a warning [0.10]
  --allow-missing   write the FASTA even if some organisms are absent
"""

import argparse
import gzip
import hashlib
import os
import re
import sys

# organism, ATCC strain, % of mix, accessions (trailing * = prefix match),
# published genome size in Mb
MANIFEST = [
    ("Acinetobacter_baumannii_ATCC17978",     "ATCC 17978",   0.18, ["CP000521", "CP000522", "CP000523"], 3.97),
    ("Bacillus_cereus_ATCC10987",             "ATCC 10987",   1.80, ["AE017194", "AE017195"],             5.22),
    ("Phocaeicola_vulgatus_ATCC8482",         "ATCC 8482",    0.02, ["CP000139"],                          5.16),
    ("Bifidobacterium_adolescentis_ATCC15703","ATCC 15703",   0.02, ["AP009256"],                          2.09),
    ("Clostridium_beijerinckii_ATCC35702",    "ATCC 35702",   1.80, ["CP006777", "NZ_CP006777"],           6.00),
    ("Cutibacterium_acnes_ATCC11828",         "ATCC 11828",   0.18, ["CP003084"],                          2.49),
    ("Deinococcus_radiodurans_ATCCBAA816",    "ATCC BAA-816", 0.02, ["AE000513", "AE001825", "AE001826", "AE001827"], 3.29),
    ("Enterococcus_faecalis_ATCC47077",       "ATCC 47077",   0.02, ["NC_017316", "CP002621"],             2.74),
    ("Escherichia_coli_ATCC700926",           "ATCC 700926",  18.0, ["U00096", "NC_000913"],               4.64),
    ("Helicobacter_pylori_ATCC700392",        "ATCC 700392",  0.18, ["AE000511"],                          1.67),
    ("Lactobacillus_gasseri_ATCC33323",       "ATCC 33323",   0.18, ["CP000413"],                          1.89),
    ("Neisseria_meningitidis_ATCCBAA335",     "ATCC BAA-335", 0.18, ["AE002098"],                          2.27),
    ("Porphyromonas_gingivalis_ATCC33277",    "ATCC 33277",   18.0, ["AP009380"],                          2.35),
    ("Pseudomonas_aeruginosa_ATCC9027",       "ATCC 9027",    1.80, ["PDLX01*"],                           6.34),
    ("Cereibacter_sphaeroides_ATCC17029",     "ATCC 17029",   18.0, ["CP000577", "CP000578", "CP000579",
                                                                     "NC_009050", "NC_009049",
                                                                     "NC_009051"],            4.37),
    ("Schaalia_odontolytica_ATCC17982",       "ATCC 17982",   0.02, ["DS2645*", "AAYI01*"],                2.39),
    ("Staphylococcus_aureus_ATCCBAA1556",     "ATCC BAA-1556",1.80, ["CP000255", "CP000256", "CP000257", "CP000258"], 2.87),
    ("Staphylococcus_epidermidis_ATCC12228",  "ATCC 12228",   18.0, ["AE015929", "AE015930", "AE015931",
                                                                     "AE015932", "AE015933", "AE015934",
                                                                     "AE015935"],                          2.50),
    ("Streptococcus_agalactiae_ATCCBAA611",   "ATCC BAA-611", 1.80, ["AE009948"],                          2.16),
    ("Streptococcus_mutans_ATCC700610",       "ATCC 700610",  18.0, ["AE014133"],                          2.03),
]

ACC_RE = re.compile(r'^>\s*([A-Za-z_]{1,4}\d{5,}(?:\.\d+)?|[A-Z]{4}\d{2}\d{6}(?:\.\d+)?)')


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def build_lookup():
    """Map exact accession -> organism, and list of (prefix, organism)."""
    exact, prefixes = {}, []
    for org, _strain, _pct, accs, _mb in MANIFEST:
        for a in accs:
            if a.endswith("*"):
                prefixes.append((normalize(a[:-1]), org))
            else:
                exact[normalize(a)] = org
    return exact, prefixes


def normalize(acc):
    """Strip version suffix and RefSeq NZ_/NC_ prefix."""
    a = acc.upper().split(".")[0]
    for pre in ("NZ_", "NC_"):
        if a.startswith(pre):
            a = a[len(pre):]
    return a


def assign(acc, exact, prefixes):
    a = normalize(acc)
    if a in exact:
        return exact[a]
    for pre, org in prefixes:
        if a.startswith(pre):
            return org
    return None


def find_fastas(root):
    exts = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")
    if os.path.isfile(root):
        return [root]
    out = []
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            if f.lower().endswith(exts):
                out.append(os.path.join(dirpath, f))
    return sorted(out)


def read_records(path):
    """Yield (header, sequence) pairs."""
    header, chunks = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line.rstrip("\n"), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", help="directory containing downloaded FASTA files")
    ap.add_argument("--output", default="msa1003_ref.fasta")
    ap.add_argument("--spacer", type=int, default=200)
    ap.add_argument("--tolerance", type=float, default=0.10)
    ap.add_argument("--allow-missing", action="store_true")
    ap.add_argument("--list-accessions", action="store_true",
                    help="print accessions for NCBI Batch Entrez and exit")
    args = ap.parse_args()

    if args.list_accessions:
        for _org, _s, _p, accs, _mb in MANIFEST:
            for a in accs:
                if not a.endswith("*") and not a.startswith("NC_"):
                    print(a)
        print("\n# Note: Pseudomonas aeruginosa ATCC 9027 (PDLX01000000) and",
              file=sys.stderr)
        print("# Schaalia odontolytica ATCC 17982 are draft assemblies; download",
              file=sys.stderr)
        print("# those as whole assemblies rather than single records.",
              file=sys.stderr)
        return 0

    if not args.input:
        ap.error("--input is required unless --list-accessions is given")

    exact, prefixes = build_lookup()
    by_org = {org: [] for org, *_ in MANIFEST}
    unassigned = []
    seen_seq = {}
    dupes = []

    files = find_fastas(args.input)
    if not files:
        print(f"No FASTA files found under {args.input}", file=sys.stderr)
        return 1
    print(f"Scanning {len(files)} FASTA file(s) under {args.input}\n")

    for path in files:
        for header, seq in read_records(path):
            m = ACC_RE.match(header)
            acc = m.group(1) if m else header[1:].split()[0]
            digest = hashlib.md5(seq.upper().encode()).hexdigest()
            if digest in seen_seq:
                dupes.append((acc, seen_seq[digest], len(seq)))
                continue     # identical sequence already taken (GCF/GCA mirror)
            org = assign(acc, exact, prefixes)
            if org:
                seen_seq[digest] = acc
                by_org[org].append((acc, seq, header))
            else:
                unassigned.append((acc, len(seq), os.path.basename(path)))

    # ---- report ----
    print(f"{'organism':44s} {'reps':>4s} {'assembled':>11s} {'published':>10s} {'diff':>7s}")
    print("-" * 82)
    missing, warned = [], []
    for org, strain, pct, _accs, mb in MANIFEST:
        recs = by_org[org]
        total = sum(len(s) for _a, s, _h in recs)
        if not recs:
            missing.append(org)
            print(f"{org:44s} {0:4d} {'MISSING':>11s} {mb*1e6:10,.0f}")
            continue
        diff = (total - mb * 1e6) / (mb * 1e6)
        flag = "" if abs(diff) <= args.tolerance else "  <-- check"
        if flag:
            warned.append((org, total, mb))
        print(f"{org:44s} {len(recs):4d} {total:11,d} {mb*1e6:10,.0f} {diff:+6.1%}{flag}")

    if dupes:
        print(f"\n{len(dupes)} duplicate sequence(s) skipped "
              f"(same DNA under a second accession):")
        for acc, kept, n in dupes[:8]:
            print(f"   {acc:18s} {n:10,d} bp   identical to {kept}")
        if len(dupes) > 8:
            print(f"   ... and {len(dupes)-8} more")

    if unassigned:
        print(f"\n{len(unassigned)} record(s) did not match any organism:")
        for acc, n, src in unassigned[:20]:
            print(f"   {acc:18s} {n:10,d} bp   ({src})")
        if len(unassigned) > 20:
            print(f"   ... and {len(unassigned)-20} more")
        print("   These are ignored. If any belong in the mix, add the accession"
              "\n   to MANIFEST and rerun.")

    if missing and not args.allow_missing:
        print(f"\n{len(missing)} organism(s) missing; not writing output.")
        print("Download them, or rerun with --allow-missing.")
        return 1

    # ---- write reference ----
    spacer = "N" * args.spacer
    n_written = 0
    with open(args.output, "w") as out:
        for org, strain, pct, _accs, _mb in MANIFEST:
            recs = by_org[org]
            if not recs:
                continue
            # longest replicon first, so the chromosome leads
            recs.sort(key=lambda r: len(r[1]), reverse=True)
            seq = spacer.join(s for _a, s, _h in recs)
            out.write(f">{org}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")
            n_written += 1

    map_path = os.path.splitext(args.output)[0] + "_contig_map.tsv"
    with open(map_path, "w") as mf:
        mf.write("accession\torganism\tstrain\tpercent_of_mix\tlength_bp\n")
        for org, strain, pct, _accs, _mb in MANIFEST:
            for acc, seq, _h in by_org[org]:
                mf.write(f"{acc}\t{org}\t{strain}\t{pct}\t{len(seq)}\n")

    total_bp = sum(sum(len(s) for _a, s, _h in by_org[o])
                   for o, *_ in MANIFEST)
    print(f"\nWrote {args.output}: {n_written} records, {total_bp:,} bp "
          f"(published total 67 Mb)")
    print(f"Wrote {map_path}")
    if warned:
        print(f"\n{len(warned)} organism(s) differ from the published size by more "
              f"than {args.tolerance:.0%} — likely a missing plasmid or the wrong "
              f"assembly. Check those before using the file.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        sys.exit(0)
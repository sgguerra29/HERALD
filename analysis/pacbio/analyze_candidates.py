#!/usr/bin/env python3
"""
Analyze HERALD candidate output files.

Usage:  python3 analyze_candidates.py results_trunc*.txt

For each candidate read, determines:
  host   = genome contributing the most fragments (ties -> higher total AS)
  donor  = the other genome(s)
  pattern= run-length encoding of the genome sequence along the read
           (AB / BA = terminal junction, ABA = flanked insert, other = complex)
  insert = number of fragments in the largest non-host run

Reports per-file and pooled summaries, plus abundance-normalized candidate rates per host genome 
so reference-quality effects can be separated from coverage effects.
"""
import sys, os
from collections import Counter, defaultdict

ABUNDANCE = {
    "Escherichia_coli_ATCC700926":            (18.0,  True),
    "Porphyromonas_gingivalis_ATCC33277":     (18.0,  True),
    "Cereibacter_sphaeroides_ATCC17029":      (18.0,  True),
    "Staphylococcus_epidermidis_ATCC12228":   (18.0,  True),
    "Streptococcus_mutans_ATCC700610":        (18.0,  True),
    "Bacillus_cereus_ATCC10987":              (1.80,  True),
    "Clostridium_beijerinckii_ATCC35702":     (1.80,  True),
    "Pseudomonas_aeruginosa_ATCC9027":        (1.80,  True),
    "Staphylococcus_aureus_ATCCBAA1556":      (1.80,  True),
    "Streptococcus_agalactiae_ATCCBAA611":    (1.80,  True),
    "Acinetobacter_baumannii_ATCC17978":      (0.18,  True),
    "Cutibacterium_acnes_ATCC11828":          (0.18,  True),
    "Helicobacter_pylori_ATCC700392":         (0.18,  True),
    "Lactobacillus_gasseri_ATCC33323":        (0.18,  True),
    "Neisseria_meningitidis_ATCCBAA335":      (0.18,  True),
    "Phocaeicola_vulgatus_ATCC8482":          (0.02,  True),
    "Bifidobacterium_adolescentis_ATCC15703": (0.02,  True),
    "Deinococcus_radiodurans_ATCCBAA816":     (0.02,  True),
    "Enterococcus_faecalis_ATCC47077":        (0.02,  True),
    "Schaalia_odontolytica_ATCC17982":        (0.02,  True),
}


def frag_num(name):
    return int(name.split("_")[-1])


def parse(path):
    """Yield lists of (frag_num, genome, AS) per candidate read."""
    entry, name = [], None
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                if entry:
                    yield name, sorted(entry)
                entry, name = [], None
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            score = 0
            if len(f) > 2 and f[2].startswith("AS:i:"):
                try:
                    score = int(f[2][5:])
                except ValueError:
                    score = 0
            entry.append((frag_num(f[0]), f[1], score))
            name = f[0].split("_frag_")[0]
    if entry:
        yield name, sorted(entry)


def runs(genomes):
    """Run-length encode: [A,A,B,B,A] -> [(A,2),(B,2),(A,1)]"""
    out = []
    for g in genomes:
        if out and out[-1][0] == g:
            out[-1][1] += 1
        else:
            out.append([g, 1])
    return [(g, n) for g, n in out]


def classify(rl, host):
    """Label the architecture using letters, host = A."""
    letters, seen = [], {host: "A"}
    nxt = ord("B")
    for g, _ in rl:
        if g not in seen:
            seen[g] = chr(nxt); nxt += 1
        letters.append(seen[g])
    return "".join(letters)


def analyze(paths):
    per_file, pooled = {}, []
    for p in paths:
        reads = []
        for name, entry in parse(p):
            genomes = [g for _, g, _ in entry]
            counts = Counter(genomes)
            score_by = defaultdict(int)
            for _, g, s in entry:
                score_by[g] += s
            host = max(counts, key=lambda g: (counts[g], score_by[g]))
            rl = runs(genomes)
            non_host = [n for g, n in rl if g != host]
            donors = sorted(set(genomes) - {host})
            host_as = [s for _, g, s in entry if g == host]
            don_as = [s for _, g, s in entry if g != host]
            reads.append(dict(
                read=name, file=os.path.basename(p), nfrag=len(entry),
                host=host, donors=donors, n_genomes=len(counts),
                pattern=classify(rl, host),
                insert=max(non_host) if non_host else 0,
                host_frac=counts[host] / len(entry),
                mean_host_as=sum(host_as) / len(host_as) if host_as else 0,
                mean_donor_as=sum(don_as) / len(don_as) if don_as else 0,
            ))
        per_file[os.path.basename(p)] = reads
        pooled += reads
    return per_file, pooled


def bar(frac, width=28):
    n = int(round(frac * width))
    return "#" * n + "." * (width - n)


def main(paths):
    per_file, pooled = analyze(paths)
    N = len(pooled)
    if not N:
        print("No candidates parsed.")
        return

    print("=" * 78)
    print("PER-FILE SUMMARY")
    print("=" * 78)
    print(f"{'file':32s} {'reads':>6s} {'frags/read':>10s} {'insert=1':>9s} {'flanked':>8s}")
    for f, reads in per_file.items():
        nf = Counter(r["nfrag"] for r in reads).most_common(1)[0][0]
        one = sum(1 for r in reads if r["insert"] == 1) / len(reads)
        flank = sum(1 for r in reads if r["pattern"].startswith("A")
                    and r["pattern"].endswith("A") and len(r["pattern"]) > 1) / len(reads)
        print(f"{f:32s} {len(reads):6d} {nf:10d} {one:8.1%} {flank:7.1%}")

    print()
    print("=" * 78)
    print(f"INSERT SIZE (fragments in largest non-host run) — n={N}")
    print("=" * 78)
    ins = Counter(r["insert"] for r in pooled)
    cum = 0
    for k in sorted(ins):
        cum += ins[k]
        print(f"  {k:2d} fragment(s): {ins[k]:5d}  {ins[k]/N:6.1%}   "
              f"cumulative {cum/N:6.1%}   {bar(ins[k]/N)}")
    drop = sum(v for k, v in ins.items() if k == 1)
    print(f"\n  Requiring insert >= 2 fragments removes {drop} of {N} "
          f"({drop/N:.1%}), leaving {N-drop}.")

    print()
    print("=" * 78)
    print("ARCHITECTURE")
    print("=" * 78)
    pat = Counter(r["pattern"] for r in pooled)
    for p, c in pat.most_common(12):
        kind = ("flanked insert" if len(p) > 2 and p[0] == "A" and p[-1] == "A"
                else "terminal junction" if len(p) == 2 else "complex")
        print(f"  {p:12s} {c:5d}  {c/N:6.1%}  {kind}")
    term = sum(c for p, c in pat.items() if len(p) == 2)
    print(f"\n  Terminal junction (AB/BA): {term} ({term/N:.1%})"
          f"   Flanked (A...A): {N-term-sum(c for p,c in pat.items() if len(p)<2)} ")

    print()
    print("=" * 78)
    print("HOST GENOME vs ABUNDANCE  (candidates per unit % of input DNA)")
    print("=" * 78)
    hosts = Counter(r["host"] for r in pooled)
    rows = []
    for g, c in hosts.items():
        ab, correct = ABUNDANCE.get(g, (None, None))
        rate = c / ab if ab else float("inf")
        rows.append((rate, g, c, ab, correct))
    rows.sort(key=lambda r: (r[0] if r[0] == r[0] else -1), reverse=True)
    print(f"{'host genome':40s} {'cands':>6s} {'%DNA':>6s} {'per %DNA':>9s}  ref?")
    for rate, g, c, ab, correct in rows:
        flag = {True: "correct", False: "SUBSTITUTE", None: "not in table"}[correct]
        if ab is None:
            print(f"{g:40s} {c:6d} {'--':>6s} {'--':>9s}  {flag}")
        else:
            print(f"{g:40s} {c:6d} {ab:6.2f} {rate:9.1f}  {flag}")

    print()
    print("-" * 78)
    print("18% TIER — matched coverage, differing reference accuracy")
    print("-" * 78)
    tier = [(g, hosts.get(g, 0), ABUNDANCE[g][1]) for g in ABUNDANCE
            if ABUNDANCE[g][0] == 18.0]
    ok = [c for _, c, corr in tier if corr]
    sub = [c for _, c, corr in tier if corr is False]
    for g, c, corr in sorted(tier, key=lambda x: -x[1]):
        print(f"  {g:40s} {c:5d}  {'correct' if corr else 'SUBSTITUTE'}")
    if ok and sub:
        mo, ms = sum(ok) / len(ok), sum(sub) / len(sub)
        print(f"\n  mean candidates, correct-reference hosts:    {mo:8.1f}")
        print(f"  mean candidates, substitute-reference hosts: {ms:8.1f}")
        if mo > 0:
            print(f"  ratio (substitute / correct):                {ms/mo:8.2f}x")

    print()
    print("=" * 78)
    print("TOP HOST -> DONOR PAIRS")
    print("=" * 78)
    pairs = Counter((r["host"], d) for r in pooled for d in r["donors"])
    for (h, d), c in pairs.most_common(15):
        print(f"  {c:5d}  {h.split('_')[0][:14]:14s} -> {d.split('_')[0][:14]:14s}"
              f"  ({c/N:5.1%})")

    print()
    print("=" * 78)
    print("ALIGNMENT SCORE: host fragments vs donor fragments")
    print("=" * 78)
    hs = sum(r["mean_host_as"] for r in pooled) / N
    ds = sum(r["mean_donor_as"] for r in pooled) / N
    lower = sum(1 for r in pooled if r["mean_donor_as"] < r["mean_host_as"])
    print(f"  mean host AS  {hs:8.0f}")
    print(f"  mean donor AS {ds:8.0f}   ({ds/hs:.1%} of host)")
    print(f"  donor scores lower than host in {lower}/{N} reads ({lower/N:.1%})")

    with open("candidates_annotated.tsv", "w") as out:
        cols = ["read", "file", "nfrag", "host", "donors", "n_genomes",
                "pattern", "insert", "host_frac", "mean_host_as", "mean_donor_as"]
        out.write("\t".join(cols) + "\n")
        for r in pooled:
            r = dict(r); r["donors"] = ",".join(r["donors"])
            out.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nWrote per-read table: candidates_annotated.tsv ({N} rows)")


if __name__ == "__main__":
    main(sys.argv[1:])
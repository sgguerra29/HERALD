import re, collections
recs = {}   # acc -> (organism, length)
acc = desc = None; n = 0
for line in open("non_redundant_Gb_bac_MINUS_Corynebacterium_jeikeium_2309.fasta"):
    if line.startswith(">"):
        if acc: recs[acc] = (desc, n)
        parts = line[1:].strip().split(None, 1)
        acc, desc, n = parts[0], (parts[1] if len(parts) > 1 else ""), 0
    else:
        n += len(line.strip())
if acc: recs[acc] = (desc, n)

big = {a: v for a, v in recs.items() if v[1] >= 3_000_000}
by_sp = collections.defaultdict(list)
by_gen = collections.defaultdict(list)
for a, (d, L) in big.items():
    w = d.split()
    if len(w) >= 2:
        by_sp[" ".join(w[:2])].append((a, L))
        by_gen[w[0]].append((a, " ".join(w[:2]), L))

print("=== same species, different strain (pair 4) ===")
for sp, v in by_sp.items():
    if len(v) >= 2: print(sp, v)

print("\n=== same genus, different species (pair 3) ===")
for g, v in by_gen.items():
    if len({s for _, s, _ in v}) >= 2: print(g, v)
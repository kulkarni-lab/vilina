#!/usr/bin/env python3
"""
Synthetic 5-taxon clade for testing the comparative (rare-genomic-change) module.

True tree:  (((T1,T2),T3),(T4,T5))

Planted insertion events, each at a fixed ancestral coordinate so orthology is known:
  e1  root                 -> present in all five            (uninformative, constant)
  e2  stem of (T1,T2)      -> T1,T2                          (informative)
  e3  stem of (T1,T2,T3)   -> T1,T2,T3                       (informative)
  e4  stem of (T4,T5)      -> T4,T5                          (informative)
  e5  terminal T3          -> T3 only                        (autapomorphy)
  e6  stem of (T1,T2,T3), secondarily DELETED in T2 -> T1,T3 (tests Dollo-style loss)
  e7  SAME viral family inserted INDEPENDENTLY at two different sites, in T1 and in T4
      -> must end up as two separate loci, never one shared character
  e8  stem of (T1,T2)      -> T1,T2  (a second one, so the branch outvotes the e6 loss)

Each taxon is then diverged from the ancestor with SNPs and a few kb-scale indels, so
coordinates differ between taxa and orthology has to be recovered by flank alignment
rather than by coordinate identity.

Usage: python test/make_compare_test_data.py <outdir>
"""
import csv
import os
import random
import sys

random.seed(1729)
OUT = sys.argv[1] if len(sys.argv) > 1 else "test/compare_data"
os.makedirs(os.path.join(OUT, "genomes"), exist_ok=True)

AA = "ACDEFGHIKLMNPQRSTVWY"
CODONS = {}
_b, _aa = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, x in enumerate(_b):
    for j, y in enumerate(_b):
        for k, z in enumerate(_b):
            CODONS.setdefault(_aa[i * 16 + j * 4 + k], []).append(x + y + z)

rprot = lambda n: "M" + "".join(random.choice(AA) for _ in range(n - 1))
rdna = lambda n: "".join(random.choice("ACGT") for _ in range(n))
backtr = lambda p: "".join(random.choice(CODONS[a]) for a in p)


def mutate_prot(p, frac):
    p = list(p)
    for i in range(1, len(p)):
        if random.random() < frac:
            p[i] = random.choice(AA)
    return "".join(p)


# ---- viral proteins -------------------------------------------------------------------
VIRUSES = {
    "Borna-like virus A": ("Riboviria; Orthornavirae; Negarnaviricota; Haploviricotina; "
                           "Monjiviricetes; Mononegavirales; Bornaviridae; Orthobornavirus",
                           [("N", "nucleoprotein", 370), ("L", "RNA-dependent RNA polymerase", 1400)]),
    "Parvo-like virus B": ("Monodnaviria; Shotokuvirae; Cossaviricota; Quintoviricetes; "
                           "Piccovirales; Parvoviridae; Densovirinae",
                           [("NS1", "NS1 protein", 600), ("VP", "capsid protein VP1", 700)]),
    "Chu-like virus C": ("Riboviria; Orthornavirae; Negarnaviricota; Haploviricotina; "
                         "Monjiviricetes; Jingchuvirales; Chuviridae",
                         [("G", "glycoprotein", 650), ("N", "nucleoprotein", 450)]),
    "Toti-like virus F": ("Riboviria; Orthornavirae; Duplornaviricota; Chrymotiviricetes; "
                          "Ghabrivirales; Totiviridae",
                          [("CP", "capsid protein", 700)]),
}
prot = {}
with open(f"{OUT}/viral.faa", "w") as fa, open(f"{OUT}/viral.gpff", "w") as gp:
    n = 0
    for org, (lin, ps) in VIRUSES.items():
        for code, product, L in ps:
            n += 1
            acc = f"YP_{800000 + n}.1"
            s = rprot(L)
            prot[(org, code)] = (acc, s)
            fa.write(f">{acc} {product} [{org}]\n{s}\n")
            wrap = [lin[i:i + 66] for i in range(0, len(lin), 66)]
            gp.write(f"LOCUS       {acc.split('.')[0]}  {L} aa  linear   VRL 01-JAN-2026\n"
                     f"DEFINITION  {product} [{org}].\nACCESSION   {acc.split('.')[0]}\n"
                     f"VERSION     {acc}\nSOURCE      {org}\n  ORGANISM  {org}\n"
                     + "".join(f"            {x}\n" for x in wrap[:-1])
                     + f"            {wrap[-1]}.\nFEATURES             Location/Qualifiers\n"
                       f"ORIGIN\n//\n")
    for k in range(25):
        n += 1
        fa.write(f">YP_{800000 + n}.1\n{rprot(random.randint(150, 800))}\n")

# ---- ancestral genome -----------------------------------------------------------------
CHRS = {"chrA": 1_600_000, "chrB": 1_100_000}
anc = {c: rdna(L) for c, L in CHRS.items()}

TAXA = ["T1", "T2", "T3", "T4", "T5"]
TRUE_TREE = "(((T1,T2),T3),(T4,T5));"

# event -> (chromosome, ancestral coord, virus, protein, taxa carrying it)
EVENTS = [
    ("e1", "chrA", 200_000, "Borna-like virus A", "N",   ["T1", "T2", "T3", "T4", "T5"]),
    ("e2", "chrA", 500_000, "Chu-like virus C",   "G",   ["T1", "T2"]),
    ("e3", "chrA", 800_000, "Parvo-like virus B", "NS1", ["T1", "T2", "T3"]),
    ("e4", "chrA", 1_100_000, "Toti-like virus F", "CP", ["T4", "T5"]),
    ("e5", "chrB", 200_000, "Parvo-like virus B", "VP",  ["T3"]),
    ("e6", "chrB", 500_000, "Borna-like virus A", "L",   ["T1", "T3"]),          # lost in T2
    ("e7a", "chrB", 800_000, "Chu-like virus C",  "N",   ["T1"]),                # independent,
    ("e7b", "chrB", 950_000, "Chu-like virus C",  "N",   ["T4"]),                # different site
    # a SECOND insertion on the (T1,T2) stem: with only one, the e6 loss character ties
    # with it exactly and the cladogram collapses -- which is the real lesson that one
    # insertion per branch is not evidence.
    ("e8", "chrA", 1_400_000, "Parvo-like virus B", "VP", ["T1", "T2"]),
]

truth = []
for eid, chrom, coord, org, code, carriers in EVENTS:
    truth.append(dict(event=eid, chrom=chrom, ancestral_coord=coord, virus=org, protein=code,
                      present_in=",".join(carriers),
                      informative=int(2 <= len(carriers) < len(TAXA))))


def eve_insert(org, code, divergence):
    """A degraded but detectable integration of this viral protein."""
    acc, s = prot[(org, code)]
    p = mutate_prot(s, divergence)
    nt = list(backtr(p))
    for _ in range(random.randint(0, 2)):                      # a stop or two
        i = random.randrange(30, len(p) - 30)
        nt[3 * i:3 * i + 3] = list("TAA")
    return "".join(nt) + "TAA"


# one sequence per event, shared by all carriers (descends from one ancestral integration),
# then diverged slightly per taxon below
event_seq = {eid: eve_insert(org, code, 0.18)
             for eid, chrom, coord, org, code, carriers in EVENTS}


def diverge(seq, snp=0.02, n_indel=6):
    """SNPs plus a few kb-scale indels, so coordinates shift between taxa."""
    s = list(seq)
    for i in range(len(s)):
        if random.random() < snp:
            s[i] = random.choice("ACGT")
    s = "".join(s)
    for _ in range(n_indel):
        pos = random.randrange(20_000, len(s) - 20_000)
        if random.random() < 0.5:
            s = s[:pos] + rdna(random.randint(500, 4000)) + s[pos:]
        else:
            s = s[:pos] + s[pos + random.randint(500, 4000):]
    return s


per_taxon_coords = {}
for t in TAXA:
    genome = {}
    coords = {}
    for chrom in CHRS:
        # insert this taxon's events into the ancestral chromosome, right to left so
        # earlier coordinates stay valid
        ev = sorted([e for e in EVENTS if e[1] == chrom and t in e[5]],
                    key=lambda e: -e[2])
        s = anc[chrom]
        marks = {}
        for eid, _c, coord, org, code, _carriers in ev:
            ins = event_seq[eid]
            # small extra divergence per taxon, so they are not byte-identical
            ins = "".join(c if random.random() > 0.03 else random.choice("ACGT") for c in ins)
            s = s[:coord] + ins + s[coord:]
            marks[eid] = (coord, len(ins))
        s = diverge(s)
        genome[chrom] = s
        coords[chrom] = marks
    per_taxon_coords[t] = coords
    with open(f"{OUT}/genomes/{t}.fa", "w") as fh:
        for c, s in genome.items():
            fh.write(f">{c}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

with open(f"{OUT}/truth_events.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(truth[0]), delimiter="\t")
    w.writeheader()
    w.writerows(truth)
open(f"{OUT}/true_tree.nwk", "w").write(TRUE_TREE + "\n")

print(f"wrote {len(TAXA)} genomes to {OUT}/genomes, truth to {OUT}/truth_events.tsv")
print("true tree:", TRUE_TREE)

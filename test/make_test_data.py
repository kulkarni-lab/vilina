#!/usr/bin/env python3
"""Synthetic test set: random 'viral' proteins with RefSeq-style .gpff metadata, a host genome
with planted EVEs of known completeness / coding state, a GTF, and host proteins.
Writes truth.tsv describing what the pipeline should recover."""
import os
import random
import sys

random.seed(7)
out = sys.argv[1] if len(sys.argv) > 1 else "test_data"
os.makedirs(out, exist_ok=True)
AA = "ACDEFGHIKLMNPQRSTVWY"
CODONS = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, x in enumerate(_b):
    for j, y in enumerate(_b):
        for k, z in enumerate(_b):
            CODONS.setdefault(_aa[i * 16 + j * 4 + k], []).append(x + y + z)


def rprot(n):
    return "M" + "".join(random.choice(AA) for _ in range(n - 1))


def rdna(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def mutate(p, frac):
    p = list(p)
    for i in range(1, len(p)):
        if random.random() < frac:
            p[i] = random.choice(AA)
    return "".join(p)


def backtr(p):
    return "".join(random.choice(CODONS[a]) for a in p)


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


VIRUSES = {
    "Borna-like virus A": ("Riboviria; Orthornavirae; Negarnaviricota; Haploviricotina; Monjiviricetes; "
                           "Mononegavirales; Bornaviridae; Orthobornavirus",
                           [("N", "nucleoprotein", 370), ("P", "phosphoprotein", 200),
                            ("G", "glycoprotein", 500), ("L", "RNA-dependent RNA polymerase", 1700)]),
    "Parvo-like virus B": ("Monodnaviria; Shotokuvirae; Cossaviricota; Quintoviricetes; Piccovirales; "
                           "Parvoviridae; Densovirinae", [("NS1", "NS1 protein", 600), ("VP", "capsid protein VP1", 700)]),
    "Chu-like virus C": ("Riboviria; Orthornavirae; Negarnaviricota; Haploviricotina; Monjiviricetes; "
                         "Jingchuvirales; Chuviridae", [("L", "polymerase", 2000), ("G", "glycoprotein", 650),
                                                        ("N", "nucleoprotein", 450)]),
    "Flavi-like virus D": ("Riboviria; Orthornavirae; Kitrinoviricota; Flasuviricetes; Amarillovirales; "
                           "Flaviviridae", [("POL", "polyprotein", 3400)]),
    "Retro-like virus E": ("Riboviria; Pararnavirae; Artverviricota; Revtraviricetes; Ortervirales; "
                           "Retroviridae", [("GAG", "gag protein", 500), ("POL", "pol polyprotein", 1100)]),
    "Toti-like virus F": ("Riboviria; Orthornavirae; Duplornaviricota; Chrymotiviricetes; Ghabrivirales; "
                          "Totiviridae", [("CP", "capsid protein", 700), ("RDRP", "RNA-dependent RNA polymerase", 800)]),
    "Giant virus G": ("Varidnaviria; Bamfordvirae; Nucleocytoviricota; Megaviricetes; Imitervirales; "
                      "Mimiviridae", [("ANK", "ankyrin repeat protein", 400), ("MCP", "major capsid protein", 550)]),
}
prot = {}
with open(f"{out}/viral.faa", "w") as fa, open(f"{out}/viral.gpff", "w") as gp:
    n = 0
    for org, (lin, ps) in VIRUSES.items():
        for code, product, L in ps:
            n += 1
            acc = f"YP_{900000 + n}.1"
            s = rprot(L)
            prot[(org, code)] = (acc, s)
            fa.write(f">{acc} {product} [{org}]\n{s}\n")
            linewrap = [lin[i:i + 66] for i in range(0, len(lin), 66)]
            gp.write(f"LOCUS       {acc.split('.')[0]}  {L} aa  linear   VRL 01-JAN-2026\n"
                     f"DEFINITION  {product} [{org}].\nACCESSION   {acc.split('.')[0]}\nVERSION     {acc}\n"
                     f"SOURCE      {org}\n  ORGANISM  {org}\n"
                     + "".join(f"            {x}\n" for x in linewrap[:-1]) + f"            {linewrap[-1]}.\n"
                     f"FEATURES             Location/Qualifiers\nORIGIN\n//\n")
    # extra unrelated viral proteins (background)
    for k in range(40):
        n += 1
        acc = f"YP_{900000 + n}.1"
        fa.write(f">{acc}\n{rprot(random.randint(150, 900))}\n")

contigs = {"chr1": 2_000_000, "chr2": 1_500_000, "chr3": 1_000_000, "scaf4": 300_000,
           "scaf5": 50_000, "scaf6": 10_000}
G = {c: list(rdna(L)) for c, L in contigs.items()}
truth = []


def plant(cid, contig, pos, strand, org, code, a0, a1, div, stops=0, fs=0, insert=0, note="",
          add_start=True):
    acc, s = prot[(org, code)]
    p = mutate(s[a0:a1], div)
    if add_start and p[0] != "M":
        p = "M" + p[1:]
    nt = backtr(p)
    nt = list(nt)
    for _ in range(stops):
        i = random.randrange(30, len(p) - 30)
        nt[3 * i:3 * i + 3] = list("TAA")
    nt = "".join(nt)
    for _ in range(fs):
        i = random.randrange(90, len(nt) - 90)
        nt = nt[:i] + nt[i + 1:]
    if insert:
        mid = len(nt) // 2
        nt = nt[:mid] + rdna(insert) + nt[mid:]
    nt = nt + "TAA"
    if strand == "-":
        nt = rc(nt)
    G[contig][pos:pos + len(nt)] = list(nt)
    truth.append(dict(id=cid, contig=contig, start=pos + 1, end=pos + len(nt), strand=strand, virus=org,
                      protein=acc, aa_from=a0 + 1, aa_to=a1, cov=round((a1 - a0) / len(s), 2),
                      stops=stops, frameshifts=fs, note=note))
    return pos + 1, pos + len(nt)


T1 = plant("T1", "chr1", 200_000, "+", "Borna-like virus A", "N", 0, 370, 0.15,
           note="complete, intact, annotated host gene -> likely functional")
plant("T10", "chr1", 203_000, "+", "Borna-like virus A", "N", 0, 370, 0.25, stops=2,
      note="tandem copy, complete length but stops -> pseudogenised")
plant("T2", "chr2", 600_000, "-", "Chu-like virus C", "G", 0, 650, 0.2, note="complete, intact, unannotated")
plant("T3", "chr2", 900_000, "+", "Parvo-like virus B", "NS1", 200, 470, 0.2, note="partial (45%), intact ORF")
plant("T4", "chr1", 900_000, "+", "Borna-like virus A", "L", 0, 1700, 0.3, stops=3, fs=2,
      note="complete-length, degraded")
T5 = plant("T5", "chr3", 300_000, "+", "Flavi-like virus D", "POL", 1000, 1600, 0.25, stops=1, insert=1500,
           note="partial, split by 1.5 kb insertion, in intron")
plant("T6", "scaf5", 60, "-", "Parvo-like virus B", "VP", 140, 700, 0.2,
      note="near-complete (80%), at contig edge")
plant("T7", "chr3", 600_000, "+", "Retro-like virus E", "POL", 0, 1100, 0.2, note="EXCLUDED Retroviridae")
plant("T8", "chr2", 1_100_000, "+", "Giant virus G", "ANK", 0, 400, 0.1, note="EXCLUDED by product regex")
plant("T9", "chr2", 1_200_000, "-", "Giant virus G", "MCP", 100, 375, 0.2, stops=1, note="partial MCP, stop")
plant("T11a", "scaf4", 50_000, "+", "Toti-like virus F", "CP", 0, 700, 0.15,
      note="Toti block: CP complete intact")
plant("T11b", "scaf4", 54_000, "+", "Toti-like virus F", "RDRP", 0, 800, 0.15, stops=1,
      note="Toti block: RdRp complete, 1 stop; block = 2/2 proteins")

with open(f"{out}/genome.fa", "w") as fh:
    for c, s in G.items():
        fh.write(f">{c} synthetic\n")
        s = "".join(s)
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")

# GTF + host proteins
gtf, pep = [], []


def gene(gid, contig, strand, exons, protein=None):
    tid, pid = f"{gid}-RA", f"{gid}-PA"
    s, e = min(x[0] for x in exons), max(x[1] for x in exons)
    a = f'gene_id "{gid}"; gene_name "{gid}";'
    gtf.append(f"{contig}\tsim\tgene\t{s}\t{e}\t.\t{strand}\t.\t{a}")
    gtf.append(f"{contig}\tsim\ttranscript\t{s}\t{e}\t.\t{strand}\t.\t{a} transcript_id \"{tid}\";")
    for x0, x1 in exons:
        gtf.append(f"{contig}\tsim\texon\t{x0}\t{x1}\t.\t{strand}\t.\t{a} transcript_id \"{tid}\";")
        gtf.append(f"{contig}\tsim\tCDS\t{x0}\t{x1}\t.\t{strand}\t0\t{a} transcript_id \"{tid}\"; protein_id \"{pid}\";")
    pep.append((pid, protein or rprot(random.randint(200, 600))))


g1 = "".join(G["chr1"][T1[0] - 1:T1[1]])
tr = "".join(_aa[_b.index(g1[i]) * 16 + _b.index(g1[i + 1]) * 4 + _b.index(g1[i + 2])]
             for i in range(0, len(g1) - 2, 3)).rstrip("*")
gene("EVEgene1", "chr1", "+", [T1], protein=tr)
gene("hostFLA", "chr3", "+", [(290_000, 292_000), (310_000, 312_000)])
for k in range(20):
    c = random.choice(["chr1", "chr2", "chr3"])
    s = random.randint(1000, contigs[c] - 20000)
    if any(abs(s - t["start"]) < 15000 for t in truth if t["contig"] == c):
        continue
    gene(f"g{k:03d}", c, random.choice("+-"), [(s, s + 400), (s + 1000, s + 1600), (s + 3000, s + 3500)])
open(f"{out}/annotation.gtf", "w").write("\n".join(gtf) + "\n")
with open(f"{out}/proteins.faa", "w") as fh:
    for pid, s in pep:
        fh.write(f">{pid}\n{s}\n")
import csv
with open(f"{out}/truth.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(truth[0]), delimiter="\t")
    w.writeheader()
    w.writerows(truth)
print(f"wrote test data to {out}")

# cellular decoy: one protein closely related to the Mimiviridae MCP (host-like homolog) + background
with open(f"{out}/cellular_decoy.faa", "w") as fh:
    fh.write(f">cell_MCP_like\n{mutate(prot[('Giant virus G', 'MCP')][1], 0.05)}\n")
    for k in range(50):
        fh.write(f">cell_{k}\n{rprot(random.randint(150, 900))}\n")
# a piRNA cluster overlapping the tandem Borna N copies
open(f"{out}/pirna_clusters.bed", "w").write("chr1\t199000\t206000\tpiC1\n")

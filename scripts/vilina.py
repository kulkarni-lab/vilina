#!/usr/bin/env python3
"""
vilina.py — Vilina (विलीन, "dissolved / merged within"): endogenous viral element (EVE)
discovery & characterisation.

Sub-commands (each is one Snakemake rule; `run` chains them all):

  prepdb     viral protein FASTA (+ RefSeq .gpff or metadata TSV) -> metadata table + DIAMOND db
  split      genome FASTA(.gz) -> normalised genome + .fai + overlapping windows for parallel blastx
  hostpep    host proteome from --pep, or translated from --cds (optional input)
  call       blastx hits -> chained, de-overlapped EVE loci; extracts nucleotide seqs and ORFs
  classify   adds ORF re-alignment, reverse-search, host-gene / GTF context, completeness and
             coding-competence verdicts; writes TSV / GFF3 / BED / FASTA / clusters
  plot       figures (PNG + PDF)
  report     self-contained HTML summary
  run        whole pipeline from a YAML config without Snakemake

Only the genome and a viral protein database are required. GTF, host proteins (pep) and
CDS are optional enrichers: they add genomic context and host-gene evidence for
domestication.
"""
import argparse
import bisect
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

# ----------------------------------------------------------------------------------------
# Generic helpers
# ----------------------------------------------------------------------------------------

DEFAULTS = {
    "threads": 4,
    "diamond": "diamond",
    "diamond_sensitivity": "--sensitive",
    "diamond_extra": "",
    "window_size": 5_000_000,
    "window_overlap": 20_000,
    "n_chunks": 8,
    # hit filters
    "max_evalue": 1e-5,
    "min_bitscore": 40,
    "min_aln_aa": 30,
    # viral DB exclusions (cellular-homolog / retrotransposon noise)
    "exclude_families": ["Retroviridae", "Metaviridae", "Pseudoviridae", "Belpaoviridae"],
    "exclude_product_regex": r"ankyrin|ubiquitin|histone|BTB|BRCT|heat.?shock|actin|tubulin|"
                             r"cytochrome|serine/threonine|tyrosine.?kinase|transposase|"
                             r"ribonucleotide.?reductase|thymidylate|cyclin|bcl-2",
    # chaining
    "chain_max_gap": 3000,
    "chain_subject_tolerance": 15,
    "max_locus_overlap": 0.2,
    # completeness
    "complete_min_cov": 0.90,
    "complete_end_tolerance_aa": 15,
    "near_complete_min_cov": 0.70,
    "contig_edge_bp": 500,
    "min_host_contig": 5000,
    # coding competence
    "flank_bp": 600,
    "min_orf_aa": 80,
    "orf_retention_intact": 0.90,
    "orf_retention_truncated": 0.50,
    # clustering of EVEs into integration blocks
    "cluster_gap": 20_000,
    # optional inputs
    "gtf": None,
    "pep": None,
    "cds": None,
    "cellular_decoy_faa": None,
    "pirna_clusters_bed": None,
}


def log(msg):
    print(f"[eve] {msg}", file=sys.stderr, flush=True)


def opener(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_fasta(path):
    """Streaming FASTA reader -> (id, seq). id = first whitespace token."""
    name, chunks = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0] if line[1:].strip() else "unnamed", []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def write_fasta(fh, name, seq, width=80):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def run_cmd(cmd):
    log("$ " + cmd)
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


COMP = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")


def revcomp(s):
    return s.translate(COMP)[::-1]


CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_b):
    for _j, _c in enumerate(_b):
        for _k, _d in enumerate(_b):
            CODON[_a + _c + _d] = _aa[_i * 16 + _j * 4 + _k]


def translate(nt):
    nt = nt.upper()
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def interval_union_len(ivs):
    if not ivs:
        return 0
    ivs = sorted(ivs)
    tot, cs, ce = 0, ivs[0][0], ivs[0][1]
    for s, e in ivs[1:]:
        if s > ce + 1:
            tot += ce - cs + 1
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + ce - cs + 1


def load_cfg(path):
    import yaml
    cfg = dict(DEFAULTS)
    if path:
        with open(path) as fh:
            user = yaml.safe_load(fh) or {}
        cfg.update({k: v for k, v in user.items() if v is not None or k not in cfg})
    return cfg


# ----------------------------------------------------------------------------------------
# FASTA index (samtools-compatible .fai) for random access on very large genomes
# ----------------------------------------------------------------------------------------

class Faidx:
    def __init__(self, fasta):
        self.fasta = fasta
        fai = fasta + ".fai"
        if not os.path.exists(fai):
            raise SystemExit(f"missing {fai}; run `vilina.py split` first")
        self.idx = {}
        self.order = []
        with open(fai) as fh:
            for line in fh:
                n, ln, off, lb, lw = line.split("\t")[:5]
                self.idx[n] = (int(ln), int(off), int(lb), int(lw))
                self.order.append(n)
        self.fh = open(fasta, "rb")

    def length(self, name):
        return self.idx[name][0]

    def fetch(self, name, start, end):
        """0-based half-open."""
        ln, off, lb, lw = self.idx[name]
        start, end = max(0, start), min(ln, end)
        if end <= start:
            return ""
        b0 = off + (start // lb) * lw + start % lb
        b1 = off + ((end - 1) // lb) * lw + (end - 1) % lb + 1
        self.fh.seek(b0)
        return self.fh.read(b1 - b0).decode().replace("\n", "").replace("\r", "").upper()


# ----------------------------------------------------------------------------------------
# prepdb
# ----------------------------------------------------------------------------------------

RNA_PLUS = {"Pisuviricota", "Kitrinoviricota", "Lenarviricota"}
RNA_MINUS = {"Negarnaviricota"}
RNA_DS = {"Duplornaviricota"}


def genome_type_from_lineage(lin):
    s = set(lin)
    if "Pararnavirae" in s or "Artverviricota" in s:
        return "RT (retro/pararetro)"
    if "Riboviria" in s or "Orthornavirae" in s:
        if s & RNA_MINUS:
            return "-ssRNA"
        if s & RNA_DS:
            return "dsRNA"
        if s & RNA_PLUS:
            return "+ssRNA"
        return "RNA (unassigned)"
    if "Monodnaviria" in s:
        return "ssDNA"
    if "Duplodnaviria" in s or "Varidnaviria" in s:
        return "dsDNA"
    if "Adnaviria" in s:
        return "dsDNA"
    return "unclassified"


def parse_gpff(path):
    """Minimal GenPept flat-file parser -> dict acc -> metadata."""
    meta = {}
    rec = None
    in_org = False
    lin_buf = []

    def flush():
        if rec is None:
            return
        lin = [x.strip().rstrip(".") for x in " ".join(lin_buf).split(";") if x.strip()]
        rec["lineage"] = "; ".join(lin)
        rec["family"] = next((x for x in reversed(lin) if x.endswith("viridae")), "unassigned")
        rec["order"] = next((x for x in reversed(lin) if x.endswith("virales")), "unassigned")
        rec["realm"] = next((x for x in lin if x.endswith("viria")), "unassigned")
        rec["genome_type"] = genome_type_from_lineage(lin)
        meta[rec["acc"]] = rec

    with opener(path) as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                f = line.split()
                rec = {"acc": f[1], "slen": int(f[2]) if f[2].isdigit() else 0,
                       "product": "", "organism": ""}
                lin_buf, in_org = [], False
            elif rec is None:
                continue
            elif line.startswith("DEFINITION"):
                d = line[12:].strip()
                m = re.match(r"(.*?)\s*\[(.*)\]\.?$", d)
                rec["product"] = m.group(1) if m else d
            elif line.startswith("VERSION"):
                rec["acc"] = line.split()[1]
            elif line.startswith("  ORGANISM"):
                rec["organism"] = line[12:].strip()
                in_org = True
            elif in_org and line.startswith("            "):
                lin_buf.append(line.strip())
            elif in_org:
                in_org = False
            if line.startswith("//"):
                flush()
                rec = None
    return meta


def cmd_prepdb(a):
    os.makedirs(a.outdir, exist_ok=True)
    meta = {}
    if a.gpff:
        for g in a.gpff:
            meta.update(parse_gpff(g))
        log(f"parsed {len(meta)} records from GenPept")
    if a.meta:
        import pandas as pd
        m = pd.read_csv(a.meta, sep="\t", dtype=str).fillna("")
        for r in m.to_dict("records"):
            meta.setdefault(r["acc"], {}).update(r)
    # fallback product/organism from RefSeq-style headers: ">acc product [organism]"
    hdr = {}
    for f in a.faa:
        with opener(f) as fh:
            for line in fh:
                if line.startswith(">"):
                    t = line[1:].strip().split(None, 1)
                    if len(t) == 2:
                        m = re.match(r"(.*?)\s*\[([^\[\]]+)\]\s*$", t[1])
                        hdr[t[0]] = {"product": m.group(1) if m else t[1],
                                     "organism": m.group(2) if m else ""}
    faa = os.path.join(a.outdir, "viral_proteins.faa")
    n = 0
    rows = []
    with open(faa, "w") as out:
        for f in a.faa:
            for name, seq in read_fasta(f):
                seq = seq.rstrip("*")
                write_fasta(out, name, seq)
                n += 1
                r = dict(hdr.get(name, {}))
                r.update({k: v for k, v in (meta.get(name) or meta.get(name.split(".")[0]) or {}).items() if v})
                rows.append({
                    "acc": name, "slen": len(seq),
                    "product": r.get("product", ""), "organism": r.get("organism", ""),
                    "family": r.get("family", "unassigned") or "unassigned",
                    "order": r.get("order", "unassigned") or "unassigned",
                    "realm": r.get("realm", "unassigned") or "unassigned",
                    "genome_type": r.get("genome_type", "unclassified") or "unclassified",
                    "lineage": r.get("lineage", ""),
                })
    import pandas as pd
    df = pd.DataFrame(rows)
    df["n_proteins_in_organism"] = df.groupby("organism")["acc"].transform("count")
    df.to_csv(os.path.join(a.outdir, "viral_meta.tsv"), sep="\t", index=False)
    log(f"{n} viral proteins; {(df.family != 'unassigned').sum()} with family assignment")
    run_cmd(f"{a.diamond} makedb --in {faa} -d {os.path.join(a.outdir, 'viral')} "
            f"--threads {a.threads} --quiet")
    if a.decoy:
        comb = os.path.join(a.outdir, "reverse_combined.faa")
        with open(comb, "w") as out:
            for name, seq in read_fasta(faa):
                write_fasta(out, "VIR|" + name, seq)
            for f in a.decoy:
                for name, seq in read_fasta(f):
                    write_fasta(out, "CEL|" + name, seq.rstrip("*"))
        run_cmd(f"{a.diamond} makedb --in {comb} -d {os.path.join(a.outdir, 'reverse')} "
                f"--threads {a.threads} --quiet")
        os.remove(comb)


# ----------------------------------------------------------------------------------------
# split
# ----------------------------------------------------------------------------------------

def cmd_split(a):
    os.makedirs(a.outdir, exist_ok=True)
    gpath = os.path.join(a.outdir, "genome.fa")
    width = 80
    chunks = [open(os.path.join(a.outdir, f"chunk_{i:03d}.fa"), "w") for i in range(a.n_chunks)]
    load = [0] * a.n_chunks
    fai = []
    off = 0
    with open(gpath, "w") as g:
        for name, seq in read_fasta(a.genome):
            seq = seq.upper()
            hdr = f">{name}\n"
            g.write(hdr)
            off += len(hdr)
            fai.append((name, len(seq), off, width, width + 1))
            for i in range(0, len(seq), width):
                line = seq[i:i + width] + "\n"
                g.write(line)
                off += len(line)
            step = a.window - a.overlap
            for ws in range(0, max(1, len(seq) - a.overlap), step):
                sub = seq[ws:ws + a.window]
                if not sub:
                    break
                k = load.index(min(load))
                write_fasta(chunks[k], f"{name}@@{ws}", sub)
                load[k] += len(sub)
                if ws + a.window >= len(seq):
                    break
    for c in chunks:
        c.close()
    with open(gpath + ".fai", "w") as fh:
        for r in fai:
            fh.write("\t".join(map(str, r)) + "\n")
    log(f"genome: {len(fai)} sequences, {sum(r[1] for r in fai):,} bp; "
        f"{a.n_chunks} chunks for blastx")


# ----------------------------------------------------------------------------------------
# hostpep
# ----------------------------------------------------------------------------------------

def cmd_hostpep(a):
    with open(a.out, "w") as out:
        if a.pep:
            for n, s in read_fasta(a.pep):
                write_fasta(out, n, s.rstrip("*").replace(".", "X"))
        elif a.cds:
            for n, s in read_fasta(a.cds):
                p = translate(s).rstrip("*")
                write_fasta(out, n, p.replace("*", "X"))
        else:
            out.write("")


# ----------------------------------------------------------------------------------------
# call: hits -> EVE loci
# ----------------------------------------------------------------------------------------

BLASTX_FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                 "qend", "sstart", "send", "evalue", "bitscore", "slen", "qlen", "qframe",
                 "qseq_translated", "cigar"]


def load_meta(path):
    import pandas as pd
    m = pd.read_csv(path, sep="\t", dtype={"acc": str}).fillna("")
    return {r["acc"]: r for r in m.to_dict("records")}


def load_hits(paths, meta, cfg):
    excl_fam = set(cfg["exclude_families"] or [])
    rx = re.compile(cfg["exclude_product_regex"], re.I) if cfg["exclude_product_regex"] else None
    seen = set()
    hits = []
    nraw = nfilt = nexcl = 0
    for p in paths:
        with open(p) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 15:
                    continue
                nraw += 1
                q, s = f[0], f[1]
                ev, bs, alen = float(f[10]), float(f[11]), int(f[3])
                if ev > cfg["max_evalue"] or bs < cfg["min_bitscore"] or alen < cfg["min_aln_aa"]:
                    nfilt += 1
                    continue
                m = meta.get(s, {})
                if m.get("family") in excl_fam or (rx and rx.search(str(m.get("product", "")))):
                    nexcl += 1
                    continue
                contig, woff = q.rsplit("@@", 1) if "@@" in q else (q, "0")
                woff = int(woff)
                qs, qe = int(f[6]), int(f[7])
                strand = "+" if qs <= qe else "-"
                gs, ge = min(qs, qe) + woff, max(qs, qe) + woff
                ss, se = sorted((int(f[8]), int(f[9])))
                key = (contig, gs, ge, s, ss, se)
                if key in seen:
                    continue
                seen.add(key)
                tr = f[15] if len(f) > 15 else ""
                cig = f[16] if len(f) > 16 else ""
                hits.append(dict(contig=contig, gstart=gs, gend=ge, strand=strand, sseqid=s,
                                 sstart=ss, send=se, slen=int(f[12]), pident=float(f[2]),
                                 alen=alen, evalue=ev, bitscore=bs,
                                 stops=tr.count("*"), frameshifts=len(re.findall(r"\d*[/\\]", cig)),
                                 aa=tr))
    log(f"hits: {nraw} raw, {nfilt} below thresholds, {nexcl} excluded by family/product, "
        f"{len(hits)} kept")
    return hits


def chain_hits(hits, cfg):
    """Collinear chaining of HSPs to the same viral protein -> candidate loci."""
    gap, tol = cfg["chain_max_gap"], cfg["chain_subject_tolerance"]
    groups = defaultdict(list)
    for h in hits:
        groups[(h["contig"], h["strand"], h["sseqid"])].append(h)
    chains = []
    for (contig, strand, sid), hs in groups.items():
        hs.sort(key=lambda h: (h["gstart"], -h["bitscore"]))
        open_chains = []
        for h in hs:
            placed = False
            for c in sorted(open_chains, key=lambda c: -c["score"]):
                last = c["hsps"][-1]
                g = h["gstart"] - c["gend"]
                if g > gap:
                    continue
                # duplicate / overlapping HSP to same subject region -> absorb
                ov = min(h["gend"], c["gend"]) - max(h["gstart"], last["gstart"])
                if ov > 0.5 * (h["gend"] - h["gstart"]):
                    if h["bitscore"] > last["bitscore"] and ov > 0.5 * (last["gend"] - last["gstart"]):
                        c["hsps"][-1] = h
                    c["gend"] = max(c["gend"], h["gend"])
                    placed = True
                    break
                collinear = (h["sstart"] >= last["send"] - tol) if strand == "+" \
                    else (h["send"] <= last["sstart"] + tol)
                if collinear:
                    c["hsps"].append(h)
                    c["gend"] = max(c["gend"], h["gend"])
                    c["score"] += h["bitscore"]
                    placed = True
                    break
            if not placed:
                open_chains.append(dict(contig=contig, strand=strand, sseqid=sid,
                                        gstart=h["gstart"], gend=h["gend"],
                                        score=h["bitscore"], hsps=[h]))
        chains.extend(open_chains)
    return chains


def resolve_overlaps(chains, cfg):
    """Greedy best-score selection; overlapping lower-scoring chains become alternatives."""
    by_contig = defaultdict(list)
    for c in sorted(chains, key=lambda c: -c["score"]):
        by_contig[c["contig"]].append(c)
    loci = []
    for contig, cs in by_contig.items():
        accepted = []
        for c in cs:
            clen = c["gend"] - c["gstart"] + 1
            host = None
            for a in accepted:
                ov = min(a["gend"], c["gend"]) - max(a["gstart"], c["gstart"]) + 1
                if ov > 0 and ov / min(clen, a["gend"] - a["gstart"] + 1) > cfg["max_locus_overlap"]:
                    host = a
                    break
            if host is None:
                c["alts"] = []
                accepted.append(c)
            else:
                host["alts"].append(c)
        loci.extend(accepted)
    return loci


def best_orf(seq, core_s, core_e, min_aa):
    """ATG..stop ORF in any frame of `seq` (already in coding orientation) with the largest
    overlap with [core_s, core_e). Returns dict or None. Also returns the longest stop-free
    stretch in the best frame overlapping the core."""
    best = None
    L = len(seq)
    for fr in range(3):
        prot = translate(seq[fr:])
        # stop positions
        pos = 0
        for seg in prot.split("*"):
            seg_start_aa = pos
            pos += len(seg) + 1
            m = seg.find("M")
            if m < 0:
                continue
            orf_aa = seg[m:]
            nt_s = fr + 3 * (seg_start_aa + m)
            nt_e = nt_s + 3 * len(orf_aa) + 3  # include stop codon if present
            has_stop = nt_e <= L and seg_start_aa + len(seg) < len(prot)
            ov = min(nt_e, core_e) - max(nt_s, core_s)
            if ov <= 0:
                continue
            if best is None or ov > best["overlap"] or (ov == best["overlap"] and len(orf_aa) > best["aa_len"]):
                best = dict(frame=fr, nt_s=nt_s, nt_e=min(nt_e, L), aa=orf_aa,
                            aa_len=len(orf_aa), overlap=ov, has_stop=has_stop)
    if best and best["aa_len"] < min_aa:
        best["below_min"] = True
    return best


def cmd_call(a):
    import pandas as pd
    cfg = load_cfg(a.config)
    meta = load_meta(a.meta)
    fa = Faidx(a.genome)
    hits = load_hits(a.hits, meta, cfg)
    chains = chain_hits(hits, cfg)
    loci = resolve_overlaps(chains, cfg)
    loci.sort(key=lambda c: (fa.order.index(c["contig"]) if c["contig"] in fa.idx else 1e9, c["gstart"]))
    log(f"{len(chains)} chains -> {len(loci)} non-overlapping EVE loci")

    os.makedirs(a.outdir, exist_ok=True)
    rows, hsp_rows = [], []
    nt_out = open(os.path.join(a.outdir, "eve_nt.fa"), "w")
    orf_out = open(os.path.join(a.outdir, "eve_orfs.faa"), "w")
    aln_out = open(os.path.join(a.outdir, "eve_aligned_aa.faa"), "w")
    flank = cfg["flank_bp"]
    for i, c in enumerate(loci, 1):
        eid = f"EVE{i:05d}"
        m = meta.get(c["sseqid"], {})
        hs = c["hsps"]
        slen = hs[0]["slen"] or int(m.get("slen", 0) or 0)
        sivs = [(h["sstart"], h["send"]) for h in hs]
        covered = interval_union_len(sivs)
        smin, smax = min(s for s, _ in sivs), max(e for _, e in sivs)
        alen = sum(h["alen"] for h in hs)
        pid = sum(h["pident"] * h["alen"] for h in hs) / alen
        clen = fa.length(c["contig"]) if c["contig"] in fa.idx else 0
        # alternative assignments
        alt_fams = Counter()
        for alt in c.get("alts", []):
            alt_fams[meta.get(alt["sseqid"], {}).get("family", "unassigned")] += 1
        fam = m.get("family", "unassigned")
        other_fams = sorted(f for f in alt_fams if f != fam)
        # sequence + ORF
        gs0, ge0 = c["gstart"] - 1, c["gend"]
        fs0, fe0 = max(0, gs0 - flank), min(clen, ge0 + flank)
        seq = fa.fetch(c["contig"], fs0, fe0)
        core = fa.fetch(c["contig"], gs0, ge0)
        if c["strand"] == "-":
            seq, core = revcomp(seq), revcomp(core)
            core_s = fe0 - ge0
        else:
            core_s = gs0 - fs0
        core_e = core_s + (ge0 - gs0)
        write_fasta(nt_out, f"{eid} {c['contig']}:{c['gstart']}-{c['gend']}({c['strand']}) "
                            f"{c['sseqid']} {fam}", core)
        aa_concat = "X".join(h["aa"].replace("/", "").replace("\\", "") for h in
                             sorted(hs, key=lambda h: h["sstart"])) if any(h["aa"] for h in hs) else ""
        if aa_concat:
            write_fasta(aln_out, eid, aa_concat.replace("*", "X"))
        orf = best_orf(seq, core_s, core_e, cfg["min_orf_aa"])
        if orf:
            # genomic coordinates of ORF
            if c["strand"] == "+":
                o_gs, o_ge = fs0 + orf["nt_s"] + 1, fs0 + orf["nt_e"]
            else:
                o_gs, o_ge = fe0 - orf["nt_e"] + 1, fe0 - orf["nt_s"]
            write_fasta(orf_out, eid, orf["aa"])
        rows.append(dict(
            eve_id=eid, contig=c["contig"], start=c["gstart"], end=c["gend"], strand=c["strand"],
            length_nt=c["gend"] - c["gstart"] + 1, contig_length=clen,
            viral_protein=c["sseqid"], viral_product=m.get("product", ""),
            virus=m.get("organism", ""), family=fam, order=m.get("order", "unassigned"),
            realm=m.get("realm", "unassigned"), genome_type=m.get("genome_type", "unclassified"),
            n_hsps=len(hs), aa_identity=round(pid, 2), aligned_aa=alen,
            best_evalue=min(h["evalue"] for h in hs), sum_bitscore=round(c["score"], 1),
            viral_protein_length=slen, viral_aa_start=smin, viral_aa_end=smax,
            viral_aa_covered=covered, viral_protein_cov=round(covered / slen, 4) if slen else 0,
            stops_in_alignment=sum(h["stops"] for h in hs),
            frameshifts_in_alignment=sum(h["frameshifts"] for h in hs),
            n_alt_assignments=len(c.get("alts", [])),
            alt_families=";".join(other_fams),
            orf_aa_len=orf["aa_len"] if orf else 0,
            orf_start=o_gs if orf else "", orf_end=o_ge if orf else "",
            orf_core_overlap=round(orf["overlap"] / max(1, core_e - core_s), 3) if orf else 0,
        ))
        for j, h in enumerate(sorted(hs, key=lambda h: h["gstart"]), 1):
            hsp_rows.append(dict(eve_id=eid, hsp=j, contig=c["contig"], start=h["gstart"],
                                 end=h["gend"], strand=c["strand"], sstart=h["sstart"],
                                 send=h["send"], pident=h["pident"], evalue=h["evalue"],
                                 bitscore=h["bitscore"], stops=h["stops"],
                                 frameshifts=h["frameshifts"]))
    for fh in (nt_out, orf_out, aln_out):
        fh.close()
    pd.DataFrame(rows).to_csv(os.path.join(a.outdir, "eves_raw.tsv"), sep="\t", index=False)
    pd.DataFrame(hsp_rows).to_csv(os.path.join(a.outdir, "eve_hsps.tsv"), sep="\t", index=False)
    log(f"wrote {len(rows)} EVE loci")


# ----------------------------------------------------------------------------------------
# classify
# ----------------------------------------------------------------------------------------

def parse_attrs(s):
    d = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', s):
        d.setdefault(m.group(1), m.group(2))
    if not d:  # GFF3 style fallback
        for kv in s.split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                d.setdefault(k.strip(), v.strip())
    return d


class GeneModel:
    def __init__(self, gtf):
        self.genes = defaultdict(list)   # contig -> [(s,e,strand,gid,name)]
        self.exons = defaultdict(list)
        self.cds = defaultdict(list)
        self.alias = {}                  # protein/transcript id -> gene id
        spans = {}
        with opener(gtf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                at = parse_attrs(f[8])
                gid = at.get("gene_id") or at.get("Parent") or at.get("ID") or ""
                tid = at.get("transcript_id", "")
                s, e = int(f[3]), int(f[4])
                feat = f[2]
                key = (f[0], gid)
                if gid:
                    sp = spans.setdefault(key, [s, e, f[6], at.get("gene_name", at.get("Name", gid))])
                    sp[0], sp[1] = min(sp[0], s), max(sp[1], e)
                for k in ("protein_id", "transcript_id"):
                    if at.get(k):
                        self.alias[at[k]] = gid
                        self.alias[at[k].split(".")[0]] = gid
                if gid:
                    self.alias[gid] = gid
                if feat == "exon":
                    self.exons[f[0]].append((s, e, f[6], gid))
                elif feat == "CDS":
                    self.cds[f[0]].append((s, e, f[6], gid))
        for (contig, gid), (s, e, st, name) in spans.items():
            self.genes[contig].append((s, e, st, gid, name))
        for d in (self.genes, self.exons, self.cds):
            for k in d:
                d[k].sort()
        self.maxlen = {k: max(x[1] - x[0] for x in v) + 1 for k, v in self.genes.items()}
        self.starts = {k: [x[0] for x in v] for k, v in self.genes.items()}

    @staticmethod
    def _ov(lst, s, e):
        return [x for x in lst if x[0] <= e and x[1] >= s]

    def context(self, contig, s, e, strand):
        genes = self.genes.get(contig, [])
        if not genes:
            return dict(context="intergenic", host_gene="", host_gene_name="",
                        orientation="", nearest_gene="", nearest_gene_dist="")
        st = self.starts[contig]
        lo = bisect.bisect_left(st, s - self.maxlen[contig])
        hi = bisect.bisect_right(st, e)
        ov = [g for g in genes[lo:hi] if g[1] >= s]
        if ov:
            gids = {g[3] for g in ov}
            cds = [x for x in self._ov(self.cds.get(contig, []), s, e) if x[3] in gids]
            exo = [x for x in self._ov(self.exons.get(contig, []), s, e) if x[3] in gids]
            ctx = "CDS" if cds else ("exon (UTR/ncRNA)" if exo else "intron")
            g = (cds or exo or [(*ov[0][:4],)])[0]
            gname = next((x[4] for x in ov if x[3] == g[3]), g[3])
            return dict(context=ctx, host_gene=g[3], host_gene_name=gname,
                        orientation="sense" if g[2] == strand else "antisense",
                        nearest_gene=g[3], nearest_gene_dist=0)
        # nearest
        cand = genes[max(0, lo - 1):min(len(genes), hi + 2)]
        best, bd = None, None
        i = bisect.bisect_left(st, s)
        for g in genes[max(0, i - 50):i + 50] + cand:
            d = s - g[1] if g[1] < s else g[0] - e
            if bd is None or d < bd:
                best, bd = g, d
        return dict(context="intergenic", host_gene="", host_gene_name="", orientation="",
                    nearest_gene=best[3] if best else "", nearest_gene_dist=bd if best else "")


def read_m8(path, cols):
    import pandas as pd
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=cols)
    return pd.read_csv(path, sep="\t", header=None, names=cols)


def cmd_classify(a):
    import pandas as pd
    cfg = load_cfg(a.config)
    meta = load_meta(a.meta)
    df = pd.read_csv(os.path.join(a.indir, "eves_raw.tsv"), sep="\t")
    for col in ("virus", "viral_product", "alt_families"):
        if col in df:
            df[col] = df[col].fillna("").astype(str)
    if "virus" in df:  # no metadata for a protein -> treat the protein as its own 'virus'
        df["virus"] = df.virus.where(df.virus != "", df.viral_protein)
    out = a.outdir
    os.makedirs(out, exist_ok=True)
    if df.empty:
        log("no EVEs found")
        for f in ("eves.tsv", "eve_clusters.tsv", "host_viral_homologs.tsv"):
            pd.DataFrame().to_csv(os.path.join(out, f), sep="\t", index=False)
        open(os.path.join(out, "eves.gff3"), "w").write("##gff-version 3\n")
        open(os.path.join(out, "eves.bed"), "w").close()
        return

    # ---- ORF re-alignment to viral proteins ------------------------------------------
    cols = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send",
            "evalue", "bitscore", "qlen", "slen"]
    orfh = read_m8(a.orf_hits, cols)
    orf_cov, orf_pid = {}, {}
    for eid, g in orfh.groupby("qseqid"):
        target = df.loc[df.eve_id == eid, "viral_protein"].iloc[0]
        gg = g[g.sseqid == target]
        if gg.empty:
            continue
        cov = interval_union_len([tuple(sorted((r.sstart, r.send))) for r in gg.itertuples()])
        orf_cov[eid] = cov / gg.slen.iloc[0]
        orf_pid[eid] = gg.pident.iloc[0]
    df["orf_viral_cov"] = df.eve_id.map(orf_cov).fillna(0).round(4)
    df["orf_retention"] = (df.orf_viral_cov / df.viral_protein_cov.clip(lower=1e-9)).clip(upper=1).round(3)

    # ---- reverse search (optional) ---------------------------------------------------
    if a.reverse_hits and os.path.exists(a.reverse_hits):
        rv = read_m8(a.reverse_hits, ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"])
        best = {}
        for eid, g in rv.groupby("qseqid"):
            g = g.sort_values("bitscore", ascending=False)
            top = g.iloc[0]
            vir = g[g.sseqid.str.startswith("VIR|")]
            cel = g[g.sseqid.str.startswith("CEL|")]
            vb = vir.bitscore.max() if not vir.empty else 0
            cb = cel.bitscore.max() if not cel.empty else 0
            best[eid] = ("viral" if top.sseqid.startswith("VIR|") else "cellular",
                         top.sseqid.split("|", 1)[1], round(vb - cb, 1))
        df["reverse_best"] = df.eve_id.map(lambda x: best.get(x, ("none", "", 0))[0])
        df["reverse_best_hit"] = df.eve_id.map(lambda x: best.get(x, ("none", "", 0))[1])
        df["reverse_bitscore_margin"] = df.eve_id.map(lambda x: best.get(x, ("none", "", 0))[2])
        df["pass_reverse"] = df.reverse_best != "cellular"
    else:
        df["reverse_best"] = "not_run"
        df["pass_reverse"] = True

    # ---- completeness ------------------------------------------------------------------
    tol = cfg["complete_end_tolerance_aa"]
    edge = cfg["contig_edge_bp"]

    def completeness(r):
        n_ok = r.viral_aa_start <= 1 + tol
        c_ok = r.viral_aa_end >= r.viral_protein_length - tol
        if r.viral_protein_cov >= cfg["complete_min_cov"] and n_ok and c_ok:
            return "complete"
        if r.viral_protein_cov >= cfg["near_complete_min_cov"]:
            return "near-complete"
        return "partial"

    df["completeness"] = df.apply(completeness, axis=1)
    df["n_terminus_present"] = df.viral_aa_start <= 1 + tol
    df["c_terminus_present"] = df.viral_aa_end >= df.viral_protein_length - tol
    df["at_contig_edge"] = (df.start <= edge) | (df.contig_length - df.end <= edge)
    # integration support: host sequence on both sides; EVEs filling most of a short contig may be
    # assembled exogenous virus (contamination / active infection) rather than integrations
    df["left_flank_bp"] = df.start - 1
    df["right_flank_bp"] = df.contig_length - df.end
    df["possible_exogenous_contig"] = (df.length_nt / df.contig_length.clip(lower=1) > 0.5) | \
                                      (df.contig_length < cfg.get("min_host_contig", 5000))

    # ---- host context ------------------------------------------------------------------
    if a.gtf:
        gm = GeneModel(a.gtf)
        ctx = [gm.context(r.contig, r.start, r.end, r.strand) for r in df.itertuples()]
        df = pd.concat([df, pd.DataFrame(ctx)], axis=1)
    else:
        gm = None
        df["context"] = "no_annotation"

    # ---- host proteins with viral homology --------------------------------------------
    hv = read_m8(a.host_hits, cols)
    host_tab = pd.DataFrame()
    host_genes_viral = {}
    if not hv.empty:
        hv = hv.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
        hv["host_gene"] = hv.qseqid.map(lambda q: (gm.alias.get(q) or gm.alias.get(q.split(".")[0]) or "")
                                        if gm else "")
        hv["viral_family"] = hv.sseqid.map(lambda s: meta.get(s, {}).get("family", "unassigned"))
        hv["viral_product"] = hv.sseqid.map(lambda s: meta.get(s, {}).get("product", ""))
        hv["virus"] = hv.sseqid.map(lambda s: meta.get(s, {}).get("organism", ""))
        hv["viral_protein_cov"] = (hv.send - hv.sstart + 1).abs() / hv.slen
        host_tab = hv.rename(columns={"qseqid": "host_protein", "sseqid": "viral_protein"})
        for r in host_tab.itertuples():
            if r.host_gene:
                host_genes_viral.setdefault(r.host_gene, set()).add(r.viral_family)
    host_tab.to_csv(os.path.join(out, "host_viral_homologs.tsv"), sep="\t", index=False)
    if gm is not None:
        df["host_gene_has_viral_hit"] = [
            bool(r.host_gene) and r.family in host_genes_viral.get(r.host_gene, set())
            for r in df.itertuples()]
    else:
        df["host_gene_has_viral_hit"] = False

    # ---- piRNA clusters (optional) ----------------------------------------------------
    if a.pirna:
        pc = defaultdict(list)
        with opener(a.pirna) as fh:
            for line in fh:
                if line.startswith(("#", "track")):
                    continue
                f = line.split("\t")
                pc[f[0]].append((int(f[1]) + 1, int(f[2])))
        df["in_piRNA_cluster"] = [any(s <= r.end and e >= r.start for s, e in pc.get(r.contig, []))
                                  for r in df.itertuples()]

    # ---- coding competence & functionality verdict -----------------------------------
    def coding(r):
        intact_aln = r.stops_in_alignment == 0 and r.frameshifts_in_alignment == 0
        if r.orf_aa_len >= cfg["min_orf_aa"] and r.orf_retention >= cfg["orf_retention_intact"] and intact_aln:
            return "intact ORF"
        if r.orf_aa_len >= cfg["min_orf_aa"] and r.orf_retention >= cfg["orf_retention_truncated"]:
            return "truncated ORF"
        return "pseudogenised"

    df["coding_status"] = df.apply(coding, axis=1)

    def verdict(r):
        annotated = r.context == "CDS" and getattr(r, "orientation", "") == "sense"
        if r.coding_status == "intact ORF" and (annotated or r.host_gene_has_viral_hit):
            return "likely functional (annotated host gene)"
        if r.coding_status == "intact ORF":
            return "potentially functional (intact ORF)"
        if r.coding_status == "truncated ORF":
            return "possibly functional (truncated ORF)"
        return "unlikely coding-functional"

    df["functional_class"] = df.apply(verdict, axis=1)
    df["evidence"] = df.apply(lambda r: "; ".join(filter(None, [
        f"ORF {r.orf_aa_len} aa retains {r.orf_retention:.0%} of viral span",
        f"{r.stops_in_alignment} stop(s)" if r.stops_in_alignment else "",
        f"{r.frameshifts_in_alignment} frameshift(s)" if r.frameshifts_in_alignment else "",
        f"{r.n_hsps} HSPs" if r.n_hsps > 1 else "",
        f"in host {r.context} ({getattr(r, 'host_gene', '')})" if r.context in ("CDS", "exon (UTR/ncRNA)", "intron") else "",
        "host protein hits same viral family" if r.host_gene_has_viral_hit else "",
        "best reverse hit cellular" if not r.pass_reverse else "",
        "EVE fills most of a short contig: check for exogenous virus" if r.possible_exogenous_contig else "",
        "at contig edge (truncation may be assembly-driven)" if r.at_contig_edge else "",
        "ambiguous: also " + r.alt_families if isinstance(r.alt_families, str) and r.alt_families else "",
    ])), axis=1)

    # ---- clusters (multi-EVE integration blocks) -------------------------------------
    df = df.sort_values(["contig", "start"]).reset_index(drop=True)
    cl_id, last = [], {}
    n = 0
    for r in df.itertuples():
        key = (r.contig, r.family)
        prev = last.get(key)
        if prev and r.start - prev[1] <= cfg["cluster_gap"]:
            cid = prev[0]
        else:
            n += 1
            cid = f"CL{n:04d}"
        last[key] = (cid, max(r.end, prev[1] if prev and prev[0] == cid else r.end))
        cl_id.append(cid)
    df["cluster"] = cl_id
    clus = []
    prot_len = defaultdict(int)
    for acc, m in meta.items():
        prot_len[m.get("organism") or acc] += int(m.get("slen", 0) or 0)
    for cid, g in df.groupby("cluster"):
        org = g.virus.mode().iloc[0]
        m0 = meta.get(g.viral_protein.iloc[0], {})
        nprot = int(m0.get("n_proteins_in_organism", 0) or 0) if m0.get("organism") else 1
        go = g[g.virus == org]
        covered = 0
        for vp, gp in go.groupby("viral_protein"):
            covered += interval_union_len([(int(r.viral_aa_start), int(r.viral_aa_end)) for r in gp.itertuples()])
        total = prot_len.get(org, 0) if org else 0
        clus.append(dict(cluster=cid, contig=g.contig.iloc[0], start=g.start.min(), end=g.end.max(),
                         span_bp=g.end.max() - g.start.min() + 1, family=g.family.iloc[0],
                         dominant_virus=org, n_eves=len(g),
                         distinct_viral_proteins=g.viral_protein.nunique(),
                         virus_proteins_in_db=nprot, virus_proteins_hit=go.viral_protein.nunique(),
                         proteome_aa_covered=covered, proteome_aa_total=total,
                         genome_completeness=round(covered / total, 3) if total else None,
                         eves=",".join(g.eve_id)))
    clus = pd.DataFrame(clus)
    if not clus.empty:
        clus["genome_completeness_class"] = pd.cut(
            clus.genome_completeness.fillna(0), [-1, 0.5, cfg["complete_min_cov"] - 1e-9, 2],
            labels=["fragmentary", "partial genome", "near-complete genome"]).astype(str)
        df = df.merge(clus[["cluster", "genome_completeness", "genome_completeness_class"]].rename(
            columns={"genome_completeness": "cluster_genome_completeness",
                     "genome_completeness_class": "cluster_genome_class"}), on="cluster", how="left")
    clus.to_csv(os.path.join(out, "eve_clusters.tsv"), sep="\t", index=False)

    df = df.sort_values("eve_id").reset_index(drop=True)
    df["pass_filters"] = df.pass_reverse
    df.to_csv(os.path.join(out, "eves.tsv"), sep="\t", index=False)

    # ---- GFF3 / BED --------------------------------------------------------------------
    hsps = pd.read_csv(os.path.join(a.indir, "eve_hsps.tsv"), sep="\t")
    with open(os.path.join(out, "eves.gff3"), "w") as g, open(os.path.join(out, "eves.bed"), "w") as b:
        g.write("##gff-version 3\n")
        for r in df.itertuples():
            note = lambda s: re.sub(r"[;=,&\t]", " ", str(s))
            attrs = (f"ID={r.eve_id};Name={r.eve_id}_{note(r.family)};viral_protein={r.viral_protein};"
                     f"virus={note(r.virus)};family={note(r.family)};product={note(r.viral_product)};"
                     f"identity={r.aa_identity};viral_cov={r.viral_protein_cov};"
                     f"completeness={r.completeness};coding_status={r.coding_status};"
                     f"functional_class={note(r.functional_class)}")
            g.write(f"{r.contig}\tEVEpipe\tviral_sequence\t{r.start}\t{r.end}\t{r.sum_bitscore}\t"
                    f"{r.strand}\t.\t{attrs}\n")
            for h in hsps[hsps.eve_id == r.eve_id].itertuples():
                g.write(f"{r.contig}\tEVEpipe\tmatch_part\t{h.start}\t{h.end}\t{h.bitscore}\t"
                        f"{r.strand}\t.\tID={r.eve_id}.hsp{h.hsp};Parent={r.eve_id};"
                        f"Target={r.viral_protein} {h.sstart} {h.send}\n")
            if r.orf_aa_len and str(r.orf_start) not in ("", "nan"):
                g.write(f"{r.contig}\tEVEpipe\tORF\t{int(r.orf_start)}\t{int(r.orf_end)}\t.\t"
                        f"{r.strand}\t0\tID={r.eve_id}.orf;Parent={r.eve_id};length_aa={r.orf_aa_len}\n")
            b.write(f"{r.contig}\t{r.start - 1}\t{r.end}\t{r.eve_id}|{r.family}|{r.completeness}\t"
                    f"{min(1000, int(r.sum_bitscore))}\t{r.strand}\n")

    summ = dict(
        n_eves=int(len(df)), n_pass=int(df.pass_filters.sum()),
        completeness=df.completeness.value_counts().to_dict(),
        functional_class=df.functional_class.value_counts().to_dict(),
        families=df.family.value_counts().to_dict(),
        n_clusters_multi=int((clus.n_eves > 1).sum()) if not clus.empty else 0,
    )
    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2)
    log(json.dumps(summ))


# ----------------------------------------------------------------------------------------
# plot
# ----------------------------------------------------------------------------------------

SURFACE = "#fcfcfb"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
OTHER = "#a3a29c"
ORD3 = {"complete": "#104281", "near-complete": "#3987e5", "partial": "#86b6ef"}
FUNC = {"likely functional (annotated host gene)": "#104281",
        "potentially functional (intact ORF)": "#256abf",
        "possibly functional (truncated ORF)": "#5598e7",
        "unlikely coding-functional": "#b7d3f6"}


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
        "font.size": 9, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "legend.frameon": False, "pdf.fonttype": 42,
    })
    return plt


def _save(fig, outdir, name):
    from matplotlib.ticker import MaxNLocator
    for ax in fig.axes:
        if ax.get_xlabel() == "EVE loci":
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        if ax.get_ylabel() == "EVE loci":
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), dpi=200, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)


def _family_colors(df, k=7):
    top = df.family.value_counts().index[:k].tolist()
    cmap = {f: CAT[i] for i, f in enumerate(top)}
    return cmap, top


def _stacked_family(df, col, order, colors, title, outdir, name):
    plt = _style()
    fams = df.family.value_counts()
    keep = fams.index[:15].tolist()
    d = df.assign(fam=df.family.where(df.family.isin(keep), "Other families"))
    tab = d.groupby(["fam", col]).size().unstack(fill_value=0).reindex(columns=order, fill_value=0)
    tab = tab.loc[tab.sum(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(7.5, 0.38 * len(tab) + 1.6))
    left = [0] * len(tab)
    for c in order:
        v = tab[c].values
        ax.barh(tab.index, v, left=left, color=colors[c], height=0.62, label=c,
                edgecolor=SURFACE, linewidth=1)
        left = [l + x for l, x in zip(left, v)]
    for i, t in enumerate(tab.sum(axis=1).values):
        ax.text(t + max(left) * 0.01, i, str(int(t)), va="center", color=INK2, fontsize=8)
    ax.set_xlabel("EVE loci")
    ax.grid(axis="y", visible=False)
    ax.set_title(title)
    ax.legend(loc="lower right", fontsize=8)
    _save(fig, outdir, name)


def cmd_plot(a):
    import pandas as pd
    import numpy as np
    df = pd.read_csv(a.eves, sep="\t")
    os.makedirs(a.outdir, exist_ok=True)
    if df.empty:
        log("nothing to plot")
        return
    if "pass_filters" in df and not a.include_failed:
        df = df[df.pass_filters]
    plt = _style()

    # 1/2 family composition
    _stacked_family(df, "completeness", list(ORD3), ORD3,
                    "EVE loci by viral family and completeness", a.outdir, "01_family_completeness")
    _stacked_family(df, "functional_class", list(FUNC), FUNC,
                    "EVE loci by viral family and coding potential", a.outdir, "02_family_functionality")

    # 3 coverage vs identity
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    gts = df.genome_type.value_counts().index.tolist()
    for i, gt in enumerate(gts):
        s = df[df.genome_type == gt]
        ax.scatter(s.viral_protein_cov * 100, s.aa_identity, s=np.clip(np.sqrt(s.length_nt) * 1.2, 10, 120),
                   color=CAT[i] if i < len(CAT) else OTHER, alpha=0.85, edgecolor=SURFACE,
                   linewidth=0.8, label=f"{gt} (n={len(s)})")
    ax.axvline(90, color=INK2, lw=0.8, ls=(0, (3, 3)))
    ax.axvline(70, color=INK2, lw=0.8, ls=(0, (1, 3)))
    ax.text(90.5, ax.get_ylim()[1], "complete", fontsize=7, color=INK2, va="top")
    ax.text(70.5, ax.get_ylim()[1], "near-complete", fontsize=7, color=INK2, va="top")
    ax.set_xlabel("Viral protein covered (%)")
    ax.set_ylabel("Amino-acid identity to best viral hit (%)")
    ax.set_xlim(0, 102)
    ax.set_title("Coverage vs identity (marker size ∝ √length)")
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
    _save(fig, a.outdir, "03_coverage_identity")

    # 4 genome map
    cmap, top = _family_colors(df)
    contigs = (df.groupby("contig").agg(n=("eve_id", "size"), L=("contig_length", "first"))
               .sort_values(["n", "L"], ascending=False).head(30))
    contigs = contigs.sort_values("L")
    fig, ax = plt.subplots(figsize=(8, 0.3 * len(contigs) + 1.6))
    for i, (c, r) in enumerate(contigs.iterrows()):
        ax.plot([0, r.L / 1e6], [i, i], color=GRID, lw=5, solid_capstyle="round", zorder=1)
        s = df[df.contig == c]
        for fam, ss in s.groupby("family"):
            ax.scatter((ss.start + ss.end) / 2e6, [i] * len(ss), marker="|", s=160, linewidths=2,
                       color=cmap.get(fam, OTHER), zorder=2)
    ax.set_yticks(range(len(contigs)))
    ax.set_yticklabels(contigs.index, fontsize=7)
    ax.set_xlabel("Position (Mb)")
    ax.grid(axis="y", visible=False)
    ttl = "EVE positions on sequences carrying EVEs"
    if df.contig.nunique() > 30:
        ttl += f" (30 of {df.contig.nunique()} shown)"
    ax.set_title(ttl)
    from matplotlib.lines import Line2D
    hand = [Line2D([], [], color=cmap[f], marker="|", ls="", markersize=10, mew=2, label=f) for f in top]
    if df.family.nunique() > len(top):
        hand.append(Line2D([], [], color=OTHER, marker="|", ls="", markersize=10, mew=2, label="Other"))
    ax.legend(handles=hand, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
    _save(fig, a.outdir, "04_genome_map")

    # 5 positional retention along viral protein
    fig, ax = plt.subplots(figsize=(6.5, 4))
    bins = np.linspace(0, 1, 101)
    for i, fam in enumerate(top[:6]):
        s = df[df.family == fam]
        prof = np.zeros(100)
        for r in s.itertuples():
            a0 = (r.viral_aa_start - 1) / r.viral_protein_length
            a1 = r.viral_aa_end / r.viral_protein_length
            prof[(bins[:-1] >= a0 - 1e-9) & (bins[1:] <= a1 + 1e-9)] += 1
        ax.plot(bins[:-1] * 100 + 0.5, prof / len(s) * 100, color=CAT[i], lw=2, label=f"{fam} (n={len(s)})")
    ax.set_xlabel("Position along viral protein (% of length, N→C)")
    ax.set_ylabel("EVEs covering position (%)")
    ax.set_ylim(0, 105)
    ax.set_title("Which parts of viral proteins are retained")
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
    _save(fig, a.outdir, "05_viral_protein_retention")

    # 6 length distribution
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    L = df.length_nt.clip(lower=1)
    bins_l = np.logspace(np.log10(L.min()), np.log10(L.max() + 1), 30)
    stack = [df[df.completeness == c].length_nt for c in ORD3]
    ax.hist(stack, bins=bins_l, stacked=True, color=list(ORD3.values()), label=list(ORD3),
            edgecolor=SURFACE, linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("EVE locus length (nt, log scale)")
    ax.set_ylabel("EVE loci")
    ax.set_title("EVE length distribution")
    ax.legend(fontsize=8)
    _save(fig, a.outdir, "06_length_distribution")

    # 7 genomic context
    if "context" in df and (df.context != "no_annotation").any():
        lab = df.context + df.get("orientation", pd.Series("", index=df.index)).fillna("").map(
            lambda o: f" ({o})" if o else "")
        vc = lab.value_counts().sort_values()
        fig, ax = plt.subplots(figsize=(6.5, 0.4 * len(vc) + 1.4))
        ax.barh(vc.index, vc.values, color=CAT[0], height=0.6)
        for i, v in enumerate(vc.values):
            ax.text(v + vc.max() * 0.01, i, str(v), va="center", color=INK2, fontsize=8)
        ax.grid(axis="y", visible=False)
        ax.set_xlabel("EVE loci")
        ax.set_title("Genomic context of EVEs (relative to host genes)")
        _save(fig, a.outdir, "07_genomic_context")

    # 8 cluster genome completeness
    if a.clusters and os.path.exists(a.clusters) and os.path.getsize(a.clusters) > 1:
        cl = pd.read_csv(a.clusters, sep="\t")
        if not cl.empty and cl.genome_completeness.notna().any():
            cl = cl[cl.genome_completeness.notna()].sort_values(
                ["genome_completeness", "n_eves"], ascending=False).head(25)
            if not cl.empty:
                cl = cl.iloc[::-1]
                fig, ax = plt.subplots(figsize=(7, 0.35 * len(cl) + 1.6))
                labels = [f"{r.cluster} {r.contig}:{r.start:,} · {str(r.dominant_virus)[:35]}" for r in cl.itertuples()]
                colors = [cmap.get(f, OTHER) for f in cl.family]
                ax.barh(labels, cl.genome_completeness * 100, color=colors, height=0.6)
                for i, r in enumerate(cl.itertuples()):
                    ax.text(r.genome_completeness * 100 + 1, i,
                            f"{r.virus_proteins_hit}/{r.virus_proteins_in_db} proteins · {r.n_eves} EVE(s)", va="center",
                            fontsize=7, color=INK2)
                ax.set_xlim(0, 118)
                ax.set_xlabel("Viral proteome (aa) represented in EVE block (%)")
                ax.grid(axis="y", visible=False)
                ax.set_title("Viral genome completeness of EVE blocks")
                ax.tick_params(axis="y", labelsize=7)
                _save(fig, a.outdir, "08_cluster_genome_completeness")
    log(f"plots written to {a.outdir}")


# ----------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------

def cmd_report(a):
    import base64
    import pandas as pd
    df = pd.read_csv(a.eves, sep="\t")
    summ = json.load(open(a.summary)) if os.path.exists(a.summary) else {}
    imgs = sorted(f for f in os.listdir(a.plots) if f.endswith(".png")) if os.path.isdir(a.plots) else []
    show = ["eve_id", "contig", "start", "end", "strand", "length_nt", "family", "virus", "viral_product",
            "aa_identity", "viral_protein_cov", "completeness", "orf_aa_len", "coding_status",
            "functional_class", "context", "evidence"]
    show = [c for c in show if c in df]
    tbl = df[show].sort_values(["viral_protein_cov"], ascending=False).head(200).to_html(
        index=False, classes="t", border=0, float_format=lambda x: f"{x:.3g}") if not df.empty else "<p>No EVEs.</p>"

    def card(k, v):
        return f'<div class="c"><div class="v">{v}</div><div class="k">{k}</div></div>'

    cards = card("EVE loci", summ.get("n_eves", 0)) + card("pass reverse filter", summ.get("n_pass", 0))
    for k in ("complete", "near-complete", "partial"):
        cards += card(k, summ.get("completeness", {}).get(k, 0))
    fc = summ.get("functional_class", {})
    cards += card("intact ORF", fc.get("potentially functional (intact ORF)", 0) +
                  fc.get("likely functional (annotated host gene)", 0))
    figs = "".join(
        f'<figure><img src="data:image/png;base64,{base64.b64encode(open(os.path.join(a.plots, f), "rb").read()).decode()}"/>'
        f'<figcaption>{f[:-4]}</figcaption></figure>' for f in imgs)
    html = f"""<!doctype html><html><head><meta charset="utf-8"><title>Vilina report</title>
<style>body{{font:14px/1.5 system-ui,sans-serif;background:#fcfcfb;color:#0b0b0b;max-width:1200px;margin:24px auto;padding:0 16px}}
h1{{font-size:22px}}.sub{{color:#52514e;font-size:13px;margin:-8px 0 16px}}.cs{{display:flex;flex-wrap:wrap;gap:12px}}.c{{border:1px solid #e4e3df;border-radius:8px;padding:10px 14px;min-width:120px}}
.v{{font-size:22px;font-weight:600}}.k{{color:#52514e;font-size:12px}}figure{{margin:16px 0}}img{{max-width:100%}}
figcaption{{color:#52514e;font-size:12px}}.t{{border-collapse:collapse;font-size:12px}}.t td,.t th{{border-bottom:1px solid #e4e3df;padding:4px 6px;text-align:left;vertical-align:top}}
.wrap{{overflow-x:auto}}</style></head><body><h1>Vilina — {os.path.basename(a.title)}</h1>
<p class="sub">विलीन · endogenous viral elements</p>
<div class="cs">{cards}</div><h2>Figures</h2>{figs}<h2>EVE loci (top 200 by coverage)</h2><div class="wrap">{tbl}</div>
<p style="color:#52514e;font-size:12px">“Functional” classes describe coding competence and host-gene annotation only;
they are hypotheses to test with expression, selection (dN/dS across orthologous EVEs) or small-RNA data.</p></body></html>"""
    open(a.out, "w").write(html)
    log(f"report: {a.out}")


# ----------------------------------------------------------------------------------------
# run (no Snakemake)
# ----------------------------------------------------------------------------------------

def cmd_run(a):
    cfg = load_cfg(a.config)
    # CLI overrides let the comparative module drive one genome at a time off a shared config
    for k in ("genome", "outdir", "sample", "gtf", "pep", "cds"):
        if getattr(a, k, None):
            cfg[k] = getattr(a, k)
    me = os.path.abspath(__file__)
    py = sys.executable
    out = cfg["outdir"]
    w = os.path.join(out, "work")
    os.makedirs(w, exist_ok=True)
    dmd, t = cfg["diamond"], cfg["threads"]
    sens = cfg["diamond_sensitivity"]
    ex = cfg.get("diamond_extra") or ""
    q = lambda x: f"'{x}'"

    if getattr(a, "db", None):
        db = a.db                       # prebuilt viral DB shared across taxa
        if not os.path.exists(os.path.join(db, "viral.dmnd")):
            raise SystemExit(f"--db {db} has no viral.dmnd")
    else:
        db = os.path.join(w, "db")
    if not getattr(a, "db", None) and (not os.path.exists(os.path.join(db, "viral.dmnd"))
                                       or a.force):
        if not cfg.get("viral_faa"):
            cfg["viral_faa"], cfg["viral_gpff"] = _download_refseq_viral(
                os.path.join(w, "refseq"))
        cmd = f"{py} {me} prepdb --faa {' '.join(map(q, _aslist(cfg['viral_faa'])))} --outdir {db} " \
              f"--diamond {dmd} --threads {t}"
        if cfg.get("viral_gpff"):
            cmd += " --gpff " + " ".join(map(q, _aslist(cfg["viral_gpff"])))
        if cfg.get("viral_meta"):
            cmd += f" --meta {q(cfg['viral_meta'])}"
        if cfg.get("cellular_decoy_faa"):
            cmd += " --decoy " + " ".join(map(q, _aslist(cfg["cellular_decoy_faa"])))
        run_cmd(cmd)
    gdir = os.path.join(w, "genome")
    run_cmd(f"{py} {me} split --genome {q(cfg['genome'])} --outdir {gdir} --n-chunks {cfg['n_chunks']} "
            f"--window {cfg['window_size']} --overlap {cfg['window_overlap']}")
    hits = []
    for i in range(cfg["n_chunks"]):
        c = os.path.join(gdir, f"chunk_{i:03d}.fa")
        o = os.path.join(w, f"blastx_{i:03d}.tsv")
        if os.path.getsize(c) > 0:
            run_cmd(f"{dmd} blastx -q {c} -d {db}/viral -o {o} {sens} -F 15 --range-culling --top 10 --max-hsps 0 "
                    f"-e {cfg['max_evalue']} --threads {t} --quiet {ex} "
                    f"--outfmt 6 {' '.join(BLASTX_FIELDS)}")
            hits.append(o)
    run_cmd(f"{py} {me} call --hits {' '.join(hits)} --meta {db}/viral_meta.tsv "
            f"--genome {gdir}/genome.fa --config {q(a.config)} --outdir {w}")
    orf_hits = os.path.join(w, "orf_vs_viral.tsv")
    _blastp(dmd, f"{w}/eve_orfs.faa", f"{db}/viral", orf_hits, t, sens)
    extra = ""
    if cfg.get("cellular_decoy_faa"):
        rv = os.path.join(w, "reverse.tsv")
        if os.path.getsize(f"{w}/eve_aligned_aa.faa"):
            run_cmd(f"{dmd} blastp -q {w}/eve_aligned_aa.faa -d {db}/reverse -o {rv} {sens} -k 25 "
                    f"-e 1e-3 --threads {t} --quiet --outfmt 6 qseqid sseqid pident length evalue bitscore")
        extra += f" --reverse-hits {rv}"
    if cfg.get("pep") or cfg.get("cds"):
        hp = os.path.join(w, "host_pep.faa")
        run_cmd(f"{py} {me} hostpep --out {hp}" + (f" --pep {q(cfg['pep'])}" if cfg.get("pep") else "")
                + (f" --cds {q(cfg['cds'])}" if cfg.get("cds") and not cfg.get("pep") else ""))
        hh = os.path.join(w, "host_vs_viral.tsv")
        _blastp(dmd, hp, f"{db}/viral", hh, t, sens)
        extra += f" --host-hits {hh}"
    if cfg.get("gtf"):
        extra += f" --gtf {q(cfg['gtf'])}"
    if cfg.get("pirna_clusters_bed"):
        extra += f" --pirna {q(cfg['pirna_clusters_bed'])}"
    run_cmd(f"{py} {me} classify --indir {w} --meta {db}/viral_meta.tsv --orf-hits {orf_hits} "
            f"--config {q(a.config)} --outdir {out}{extra}")
    for f in ("eve_nt.fa", "eve_orfs.faa"):
        shutil.copy(os.path.join(w, f), os.path.join(out, f))
    if getattr(a, "skip_plots", False):
        return
    run_cmd(f"{py} {me} plot --eves {out}/eves.tsv --clusters {out}/eve_clusters.tsv --outdir {out}/plots")
    run_cmd(f"{py} {me} report --eves {out}/eves.tsv --summary {out}/summary.json --plots {out}/plots "
            f"--title {q(cfg.get('sample', 'genome'))} --out {out}/report.html")


def _download_refseq_viral(d):
    import urllib.request
    base = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
    os.makedirs(d, exist_ok=True)
    listing = urllib.request.urlopen(base).read().decode()
    files = sorted(set(re.findall(r"viral\.\d+\.protein\.(?:faa|gpff)\.gz", listing)))
    for f in files:
        dest = os.path.join(d, f)
        if not os.path.exists(dest):
            log(f"downloading {base}{f}")
            urllib.request.urlretrieve(base + f, dest)
    faa = [os.path.join(d, f) for f in files if ".faa." in f]
    gpff = [os.path.join(d, f) for f in files if ".gpff." in f]
    if not faa:
        raise SystemExit("could not list RefSeq viral files; set viral_faa/viral_gpff in the config")
    return faa, gpff


def _aslist(x):
    return x if isinstance(x, list) else [x]


def _blastp(dmd, q, db, out, t, sens):
    if os.path.exists(q) and os.path.getsize(q) > 0:
        run_cmd(f"{dmd} blastp -q {q} -d {db} -o {out} {sens} -k 25 -e 1e-5 --threads {t} --quiet "
                f"--outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen")
    else:
        open(out, "w").close()


# ----------------------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)

    s = sp.add_parser("prepdb")
    s.add_argument("--faa", nargs="+", required=True)
    s.add_argument("--gpff", nargs="*")
    s.add_argument("--meta")
    s.add_argument("--decoy", nargs="*")
    s.add_argument("--outdir", required=True)
    s.add_argument("--diamond", default="diamond")
    s.add_argument("--threads", type=int, default=4)
    s.set_defaults(fn=cmd_prepdb)

    s = sp.add_parser("split")
    s.add_argument("--genome", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--n-chunks", type=int, default=8)
    s.add_argument("--window", type=int, default=5_000_000)
    s.add_argument("--overlap", type=int, default=20_000)
    s.set_defaults(fn=cmd_split)

    s = sp.add_parser("hostpep")
    s.add_argument("--pep")
    s.add_argument("--cds")
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_hostpep)

    s = sp.add_parser("call")
    s.add_argument("--hits", nargs="+", required=True)
    s.add_argument("--meta", required=True)
    s.add_argument("--genome", required=True)
    s.add_argument("--config")
    s.add_argument("--outdir", required=True)
    s.set_defaults(fn=cmd_call)

    s = sp.add_parser("classify")
    s.add_argument("--indir", required=True)
    s.add_argument("--meta", required=True)
    s.add_argument("--orf-hits")
    s.add_argument("--reverse-hits")
    s.add_argument("--host-hits")
    s.add_argument("--gtf")
    s.add_argument("--pirna")
    s.add_argument("--config")
    s.add_argument("--outdir", required=True)
    s.set_defaults(fn=cmd_classify)

    s = sp.add_parser("plot")
    s.add_argument("--eves", required=True)
    s.add_argument("--clusters")
    s.add_argument("--outdir", required=True)
    s.add_argument("--include-failed", action="store_true")
    s.set_defaults(fn=cmd_plot)

    s = sp.add_parser("report")
    s.add_argument("--eves", required=True)
    s.add_argument("--summary", required=True)
    s.add_argument("--plots", required=True)
    s.add_argument("--title", default="genome")
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_report)

    s = sp.add_parser("run", help="whole single-genome pipeline from a YAML config")
    s.add_argument("--config", required=True)
    s.add_argument("--force", action="store_true")
    # overrides, so one config can drive many genomes (used by `compare`)
    s.add_argument("--genome", help="override config `genome`")
    s.add_argument("--outdir", help="override config `outdir`")
    s.add_argument("--sample", help="override config `sample`")
    s.add_argument("--gtf")
    s.add_argument("--pep")
    s.add_argument("--cds")
    s.add_argument("--db", help="use a prebuilt viral DIAMOND db directory instead of "
                                "building one (shares it across taxa)")
    s.add_argument("--skip-plots", action="store_true")
    s.set_defaults(fn=cmd_run)

    s = sp.add_parser("compare", help="comparative module: EVEs as rare genomic changes "
                                      "across many genomes (see vilina_compare.py)")
    s.add_argument("--config", required=True)
    s.add_argument("--genomes-dir", dest="genomes_dir")
    s.add_argument("--outdir", dest="compare_outdir")
    s.add_argument("--flank-mapper", dest="flank_mapper", choices=["minimap2", "blastn"])
    s.add_argument("--minimap2-preset", dest="minimap2_preset")
    s.add_argument("--jobs", type=int)
    s.add_argument("--skip-per-taxon", action="store_true")
    s.add_argument("--force", action="store_true")
    s.set_defaults(fn=lambda a: __import__("vilina_compare").cmd_compare(a))

    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()

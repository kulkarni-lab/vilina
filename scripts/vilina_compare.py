#!/usr/bin/env python3
"""
vilina_compare.py — Vilina comparative module: EVEs as rare genomic changes (RGCs).

An independent germline integration at a particular genomic position is effectively a
unique, unrepeatable event, so an *orthologous* integration shared by two taxa is evidence
of common ancestry — the same logic used for SINE/LINE insertion phylogenetics. The hard
part is not finding EVEs but proving orthology of the insertion SITE, and distinguishing a
true absence (ancestral state / secondary deletion) from an unassembled one.

Steps (all driven from the same config.yaml):
  1. discover    one taxon per genome FASTA in `genomes_dir`
  2. per-taxon   run the single-genome Vilina pipeline on each, sharing one viral DIAMOND DB
  3. flanks      pull the sequence flanking every EVE in every taxon
  4. map         map all flanks onto every other genome (minimap2, or blastn for deep
                 divergence) — one mapping run per target genome
  5. score       for each (EVE, other taxon): the orthologous site is located by the two
                 flanks, then scored
                    1  shared insertion  — an EVE of the same viral family sits at that site
                    0  empty site        — the flanks map ADJACENTLY, so the site exists and
                                           is demonstrably empty
                    ?  not determined    — flanks missing, inconsistent, or spanning a gap
                                           large enough to hide something
  6. orthogroup  link reciprocal shared calls into orthologous insertion loci = characters
  7. matrix      binary character matrix: TSV + NEXUS (PAUP* Dollo, MrBayes) + FASTA (IQ-TREE)
  8. splits      insertions supporting each bipartition, with four-gamete conflict check
  9. plots       pairwise EVE dotplots coloured by virus, ribbon plot in a custom taxon order,
                 presence matrix, split support, pairwise shared counts

Only step 5's "0" lets a character be informative, and only step 6's reciprocity keeps
convergent integrations of the same virus family at *different* sites from being scored as
homologous. Both are deliberately conservative: an undetermined site costs resolution,
a wrongly homologised one costs a wrong tree.
"""
import argparse
import bisect
import itertools
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vilina as V  # noqa: E402
import vilina_tree as VT  # noqa: E402

CDEFAULTS = {
    # ---- inputs -------------------------------------------------------------------
    "genomes_dir": None,          # one genome FASTA per taxon
    "taxon_order": None,          # display/ribbon order; default = alphabetical
    "taxon_names": None,          # {file stem: display name}
    "annotations_dir": None,      # optional <taxon>.gtf / <taxon>.pep.faa alongside
    "compare_outdir": None,       # default <outdir>_compare
    "reuse_per_taxon": True,      # keep existing per-taxon results instead of recomputing
    "compare_jobs": 1,            # per-taxon pipelines to run at once
    # ---- orthology ----------------------------------------------------------------
    "flank_bp_ortho": 3000,       # flank pulled either side of each EVE
    "min_flank_aln": 500,         # bp of flank that must align in the target
    "min_flank_identity": 0.75,
    "max_site_span": 50000,       # flanks further apart than this = rearranged, not a site
    "empty_site_max": 300,        # flanks this close with no EVE = scored as empty site (0)
    "site_overlap_slop": 1000,    # tolerance when asking "is a target EVE at this site"
    "allow_single_flank": True,   # one-sided anchors may support presence, never absence
    "require_reciprocal": True,   # a shared call must be found from both directions
    "same_family_required": True,
    # ---- mapping ------------------------------------------------------------------
    "flank_mapper": "minimap2",   # or "blastn" (slower, better across deep divergence)
    "minimap2": "minimap2",
    "minimap2_preset": "asm20",   # asm5/asm10/asm20 by increasing divergence
    "minimap2_extra": "",
    "blastn": "blastn",
    "makeblastdb": "makeblastdb",
    "blast_task": "dc-megablast",
    # ---- plots --------------------------------------------------------------------
    # ---- cladogram ----------------------------------------------------------------
    "build_tree": True,
    "outgroup": None,             # taxon (or list) to root on; default = last in taxon_order
    "bootstrap": 100,
    "exhaustive_max_trees": 200_000,
    "tree_search_restarts": 12,
    "max_optimal_trees": 200,
    # ---- plots ----------------------------------------------------------------------
    "ribbon_scale": "proportional",   # or "normalised"
    "dotplot_pairs": "all",           # "all", "adjacent", or [[A, B], ...]
    "max_dotplot_pairs": 30,
}

GENOME_EXT = (".fa", ".fasta", ".fna", ".fas", ".seq",
              ".fa.gz", ".fasta.gz", ".fna.gz", ".fas.gz")


def log(m):
    V.log(m)


def load_compare_cfg(path):
    import yaml
    cfg = dict(V.DEFAULTS)
    cfg.update(CDEFAULTS)
    if path:
        with open(path) as fh:
            user = yaml.safe_load(fh) or {}
        cfg.update({k: v for k, v in user.items() if v is not None or k not in cfg})
    if not cfg.get("compare_outdir"):
        cfg["compare_outdir"] = str(cfg.get("outdir", "results")).rstrip("/") + "_compare"
    return cfg


def safe_name(s):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(s))


# ----------------------------------------------------------------------------------------
# 1. discover genomes
# ----------------------------------------------------------------------------------------

def discover_genomes(cfg):
    d = cfg.get("genomes_dir")
    if not d:
        raise SystemExit("comparative module needs `genomes_dir:` in the config")
    if not os.path.isdir(d):
        raise SystemExit(f"genomes_dir not found: {d}")
    names = cfg.get("taxon_names") or {}
    found = {}
    for f in sorted(os.listdir(d)):
        low = f.lower()
        if not low.endswith(GENOME_EXT):
            continue
        stem = f
        for e in sorted(GENOME_EXT, key=len, reverse=True):
            if low.endswith(e):
                stem = f[: -len(e)]
                break
        t = safe_name(names.get(stem, names.get(f, stem)))
        if t in found:
            raise SystemExit(f"two genome files map to taxon '{t}': {found[t]} and {f}")
        found[t] = os.path.join(d, f)
    if len(found) < 2:
        raise SystemExit(f"need >= 2 genomes in {d}; found {len(found)}")
    order = cfg.get("taxon_order")
    if order:
        order = [safe_name(x) for x in order]
        missing = [t for t in order if t not in found]
        extra = [t for t in found if t not in order]
        if missing:
            raise SystemExit(f"taxon_order names with no genome file: {', '.join(missing)}")
        if extra:
            log(f"WARNING: not in taxon_order, appended alphabetically: {', '.join(sorted(extra))}")
            order = order + sorted(extra)
    else:
        order = sorted(found)
    log(f"{len(found)} taxa: {', '.join(order)}")
    return order, found


# ----------------------------------------------------------------------------------------
# 2. per-taxon EVE calling (shared viral DB)
# ----------------------------------------------------------------------------------------

def per_taxon_dir(cfg, taxon):
    return os.path.join(cfg["compare_outdir"], "per_taxon", taxon)


def run_per_taxon(cfg, order, genomes, config_path, force=False):
    out = cfg["compare_outdir"]
    work = os.path.join(out, "work")
    db = os.path.join(work, "db")
    os.makedirs(work, exist_ok=True)
    me = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vilina.py")
    py = sys.executable

    if force or not os.path.exists(os.path.join(db, "viral.dmnd")):
        faa = cfg.get("viral_faa")
        if not faa:
            faa, gpff = V._download_refseq_viral(os.path.join(work, "refseq"))
        else:
            faa, gpff = V._aslist(faa), V._aslist(cfg.get("viral_gpff") or [])
        cmd = (f"{py} {me} prepdb --faa {' '.join(map(_q, faa))} --outdir {_q(db)} "
               f"--diamond {cfg['diamond']} --threads {cfg['threads']}")
        if gpff:
            cmd += " --gpff " + " ".join(map(_q, gpff))
        if cfg.get("viral_meta"):
            cmd += f" --meta {_q(cfg['viral_meta'])}"
        if cfg.get("cellular_decoy_faa"):
            cmd += " --decoy " + " ".join(map(_q, V._aslist(cfg["cellular_decoy_faa"])))
        V.run_cmd(cmd)
    else:
        log("reusing existing shared viral DIAMOND database")

    jobs, pending = [], []
    for t in order:
        tdir = per_taxon_dir(cfg, t)
        if cfg.get("reuse_per_taxon", True) and not force and \
                os.path.exists(os.path.join(tdir, "eves.tsv")):
            log(f"reusing per-taxon results for {t}")
            continue
        cmd = (f"{py} {me} run --config {_q(config_path)} --genome {_q(genomes[t])} "
               f"--outdir {_q(tdir)} --sample {_q(t)} --db {_q(db)} --skip-plots")
        ann = cfg.get("annotations_dir")
        if ann:
            for key, exts in (("gtf", (".gtf", ".gff3", ".gff")),
                              ("pep", (".pep.faa", ".faa", ".pep.fa")),
                              ("cds", (".cds.fa", ".cds.fna"))):
                for e in exts:
                    p = os.path.join(ann, t + e)
                    if os.path.exists(p):
                        cmd += f" --{key} {_q(p)}"
                        break
        pending.append((t, cmd))

    n = max(1, int(cfg.get("compare_jobs", 1)))
    if n == 1 or len(pending) < 2:
        for t, cmd in pending:
            log(f"--- EVEs in {t} ---")
            V.run_cmd(cmd)
    else:
        log(f"running {len(pending)} per-taxon pipelines, {n} at a time")
        running = []
        for t, cmd in pending:
            while len(running) >= n:
                running = _reap(running)
            logf = open(os.path.join(cfg["compare_outdir"], "work", f"{t}.log"), "w")
            running.append((t, subprocess.Popen(cmd, shell=True, stdout=logf, stderr=logf,
                                                executable="/bin/bash"), logf))
        while running:
            running = _reap(running)
    return db


def _reap(running):
    import time
    alive = []
    for t, p, logf in running:
        rc = p.poll()
        if rc is None:
            alive.append((t, p, logf))
        else:
            logf.close()
            if rc != 0:
                raise SystemExit(f"per-taxon pipeline failed for {t} (see work/{t}.log)")
            log(f"finished {t}")
    if len(alive) == len(running):
        time.sleep(2)
    return alive


def _q(x):
    return "'" + str(x) + "'"


def taxon_genome(cfg, taxon):
    return os.path.join(per_taxon_dir(cfg, taxon), "work", "genome", "genome.fa")


def load_taxon_eves(cfg, order):
    import pandas as pd
    tabs = {}
    for t in order:
        p = os.path.join(per_taxon_dir(cfg, t), "eves.tsv")
        if not os.path.exists(p):
            raise SystemExit(f"missing {p}; run without --skip-per-taxon first")
        df = pd.read_csv(p, sep="\t")
        if df.empty:
            log(f"WARNING: no EVEs called in {t}")
            df = pd.DataFrame(columns=["eve_id", "contig", "start", "end", "strand", "family",
                                       "viral_protein", "length_nt", "completeness",
                                       "functional_class", "virus"])
        if "pass_filters" in df:
            df = df[df.pass_filters].copy()
        df["taxon"] = t
        df["key"] = t + "::" + df.eve_id.astype(str)
        tabs[t] = df.reset_index(drop=True)
    log("EVEs per taxon: " + ", ".join(f"{t}={len(tabs[t])}" for t in order))
    return tabs


# ----------------------------------------------------------------------------------------
# 3. flanks
# ----------------------------------------------------------------------------------------

def extract_flanks(cfg, order, tabs):
    """Left/right flank of every EVE. Named <taxon>::<eve_id>::L|R::<flank_len>."""
    wdir = os.path.join(cfg["compare_outdir"], "work", "flanks")
    os.makedirs(wdir, exist_ok=True)
    flank = int(cfg["flank_bp_ortho"])
    paths, lens = {}, {}
    for t in order:
        fa = V.Faidx(taxon_genome(cfg, t))
        p = os.path.join(wdir, f"{t}.fa")
        n = 0
        with open(p, "w") as out:
            for r in tabs[t].itertuples():
                clen = fa.length(r.contig)
                ls, le = max(0, r.start - 1 - flank), r.start - 1
                rs, re_ = r.end, min(clen, r.end + flank)
                for side, (s, e) in (("L", (ls, le)), ("R", (rs, re_))):
                    if e - s < int(cfg["min_flank_aln"]):
                        continue           # too close to the contig end to be usable
                    seq = fa.fetch(r.contig, s, e)
                    if len(seq) - seq.count("N") < int(cfg["min_flank_aln"]):
                        continue           # mostly assembly gap
                    V.write_fasta(out, f"{t}::{r.eve_id}::{side}::{len(seq)}", seq)
                    n += 1
                    lens[(t, r.eve_id, side)] = len(seq)
        paths[t] = p
        log(f"{t}: {n} usable flanks from {len(tabs[t])} EVEs")
    return paths, lens


# ----------------------------------------------------------------------------------------
# 4. map flanks onto every other genome — one run per target
# ----------------------------------------------------------------------------------------

def _parse_paf(path):
    hits = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            alen = int(f[10]) or 1
            hits.append(dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
                             t=f[5], ts=int(f[7]), te=int(f[8]),
                             ident=int(f[9]) / alen, alen=alen, score=int(f[9])))
    return hits


def _parse_blast(path):
    hits = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            qs, qe = int(f[2]) - 1, int(f[3])
            ss, se = int(f[6]), int(f[7])
            strand = "+" if ss <= se else "-"
            ts, te = (ss - 1, se) if strand == "+" else (se - 1, ss)
            alen = int(f[9])
            hits.append(dict(q=f[0], qlen=int(f[1]), qs=qs, qe=qe, strand=strand,
                             t=f[4], ts=ts, te=te, ident=float(f[8]) / 100.0,
                             alen=alen, score=float(f[10])))
    return hits


def map_flanks(cfg, order, flank_paths):
    """For each target taxon, map every other taxon's flanks onto it in one run."""
    wdir = os.path.join(cfg["compare_outdir"], "work", "map")
    os.makedirs(wdir, exist_ok=True)
    mapper = cfg.get("flank_mapper", "minimap2")
    out = {}
    for target in order:
        q = os.path.join(wdir, f"query_vs_{target}.fa")
        with open(q, "w") as o:
            for t in order:
                if t == target:
                    continue
                if os.path.getsize(flank_paths[t]):
                    with open(flank_paths[t]) as fh:
                        shutil.copyfileobj(fh, o)
        res = os.path.join(wdir, f"{target}." + ("paf" if mapper == "minimap2" else "tsv"))
        if os.path.getsize(q) == 0:
            open(res, "w").close()
            out[target] = []
            continue
        g = taxon_genome(cfg, target)
        if mapper == "minimap2":
            V.run_cmd(f"{cfg['minimap2']} -x {cfg['minimap2_preset']} -c --secondary=yes "
                      f"-N 5 -p 0.3 -t {cfg['threads']} {cfg.get('minimap2_extra','')} "
                      f"{_q(g)} {_q(q)} > {_q(res)} 2> {_q(res + '.log')}")
            out[target] = _parse_paf(res)
        else:
            dbp = os.path.join(wdir, f"{target}.blastdb")
            if not os.path.exists(dbp + ".nsq") and not os.path.exists(dbp + ".nin"):
                V.run_cmd(f"{cfg['makeblastdb']} -in {_q(g)} -dbtype nucl -out {_q(dbp)} "
                          f"> /dev/null")
            V.run_cmd(f"{cfg['blastn']} -task {cfg['blast_task']} -query {_q(q)} -db {_q(dbp)} "
                      f"-num_threads {cfg['threads']} -evalue 1e-10 -max_target_seqs 5 "
                      f'-outfmt "6 qseqid qlen qstart qend sseqid slen sstart send pident '
                      f'length bitscore" -out {_q(res)}')
            out[target] = _parse_blast(res)
        log(f"flank alignments onto {target}: {len(out[target])}")
    return out


# ----------------------------------------------------------------------------------------
# 5. locate and score orthologous sites
# ----------------------------------------------------------------------------------------

def _inner_edge(h, side, flank_len):
    """Target coordinate of the flank edge that touches the EVE, extrapolating over any
    unaligned tail so a partially aligned flank still points at the right place."""
    if side == "L":                      # inner edge is the 3' end of the left flank
        tail = flank_len - h["qe"]
        return h["te"] + tail if h["strand"] == "+" else h["ts"] - tail
    tail = h["qs"]                       # inner edge is the 5' start of the right flank
    return h["ts"] - tail if h["strand"] == "+" else h["te"] + tail


def build_eve_index(tabs):
    """taxon -> contig -> (sorted starts, [(start, end, eve_id, family, protein)])"""
    idx = {}
    for t, df in tabs.items():
        per = defaultdict(list)
        for r in df.itertuples():
            per[r.contig].append((int(r.start), int(r.end), r.eve_id,
                                  str(r.family), str(r.viral_protein)))
        idx[t] = {c: (sorted(v), [x[0] for x in sorted(v)]) for c, v in per.items()}
    return idx


def eves_overlapping(idx, taxon, contig, s, e):
    d = idx.get(taxon, {}).get(contig)
    if not d:
        return []
    items, starts = d
    i = bisect.bisect_left(starts, s - 10_000_000)
    return [x for x in items[max(0, i):] if x[0] <= e and x[1] >= s]


def score_sites(cfg, order, tabs, hits_by_target, flank_lens):
    """One row per (query EVE, target taxon): where the orthologous site is and what's in it."""
    idx = build_eve_index(tabs)
    fam = {}
    elen = {}
    for t, df in tabs.items():
        for r in df.itertuples():
            fam[(t, r.eve_id)] = str(r.family)
            elen[(t, r.eve_id)] = int(r.length_nt)

    min_aln, min_id = int(cfg["min_flank_aln"]), float(cfg["min_flank_identity"])
    slop = int(cfg["site_overlap_slop"])
    obs = []
    for target, hits in hits_by_target.items():
        # best hits per (query flank), keeping a few alternatives
        per_flank = defaultdict(list)
        for h in hits:
            if h["alen"] < min_aln or h["ident"] < min_id:
                continue
            q = h["q"].split("::")
            if len(q) != 4:
                continue
            per_flank[(q[0], q[1], q[2])].append(h)
        pairs = defaultdict(dict)
        for (qt, eid, side), hs in per_flank.items():
            pairs[(qt, eid)][side] = sorted(hs, key=lambda x: -x["score"])[:5]

        for (qt, eid), sides in pairs.items():
            L, R = sides.get("L", []), sides.get("R", [])
            lf = flank_lens.get((qt, eid, "L"), int(cfg["flank_bp_ortho"]))
            rf = flank_lens.get((qt, eid, "R"), int(cfg["flank_bp_ortho"]))
            best = None
            for hl, hr in itertools.product(L, R):
                if hl["t"] != hr["t"] or hl["strand"] != hr["strand"]:
                    continue
                a, b = _inner_edge(hl, "L", lf), _inner_edge(hr, "R", rf)
                lo, hi = (a, b) if hl["strand"] == "+" else (b, a)
                if hi < lo - slop:
                    continue                      # inverted / rearranged relative to query
                lo, hi = min(lo, hi), max(lo, hi)
                if hi - lo > int(cfg["max_site_span"]):
                    continue
                sc = hl["score"] + hr["score"]
                if best is None or sc > best["score"]:
                    best = dict(score=sc, contig=hl["t"], lo=int(lo), hi=int(hi),
                                strand=hl["strand"], flanks="LR",
                                ident=round(min(hl["ident"], hr["ident"]), 3))
            if best is None and cfg.get("allow_single_flank", True) and (L or R):
                h = sorted(L + R, key=lambda x: -x["score"])[0]
                side = "L" if h in L else "R"
                anchor = _inner_edge(h, side, lf if side == "L" else rf)
                w = elen.get((qt, eid), 1000) + slop
                forward = (side == "L") == (h["strand"] == "+")
                lo, hi = (anchor, anchor + w) if forward else (anchor - w, anchor)
                best = dict(score=h["score"], contig=h["t"], lo=int(lo), hi=int(hi),
                            strand=h["strand"], flanks=side, ident=round(h["ident"], 3))
            if best is None:
                obs.append(dict(query_taxon=qt, query_eve=eid, target_taxon=target,
                                status="no_anchor", target_eve="", contig="", site_start="",
                                site_end="", span="", flanks="", flank_identity=""))
                continue

            span = max(0, best["hi"] - best["lo"])
            cands = eves_overlapping(idx, target, best["contig"],
                                     best["lo"] - slop, best["hi"] + slop)
            f0 = fam.get((qt, eid), "unassigned")
            if cfg.get("same_family_required", True):
                cands = [c for c in cands if c[3] == f0 or
                         (c[3] == "unassigned" and f0 == "unassigned")]
            if cands:
                hit = max(cands, key=lambda c: min(c[1], best["hi"]) - max(c[0], best["lo"]))
                status, teve = "shared", hit[2]
            else:
                # absence is only credible when the flanks come back essentially adjacent
                el = elen.get((qt, eid), 1000)
                thresh = min(int(cfg["empty_site_max"]), max(50, int(0.25 * el)))
                if best["flanks"] == "LR" and span <= thresh:
                    status, teve = "empty", ""
                else:
                    status, teve = "ambiguous", ""
            obs.append(dict(query_taxon=qt, query_eve=eid, target_taxon=target, status=status,
                            target_eve=teve, contig=best["contig"], site_start=best["lo"],
                            site_end=best["hi"], span=span, flanks=best["flanks"],
                            flank_identity=best["ident"]))
    import pandas as pd
    df = pd.DataFrame(obs)
    if not df.empty:
        log("site calls: " + ", ".join(f"{k}={v}" for k, v in
                                       df.status.value_counts().items()))
    return df


# ----------------------------------------------------------------------------------------
# 6. orthogroups
# ----------------------------------------------------------------------------------------

class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def build_orthogroups(cfg, order, tabs, obs):
    shared = {}
    for r in obs[obs.status == "shared"].itertuples():
        shared[(r.query_taxon, r.query_eve, r.target_taxon)] = r.target_eve
    uf = UF()
    for t in order:
        for r in tabs[t].itertuples():
            uf.find((t, r.eve_id))
    n_edges = n_oneway = 0
    for (qt, qe, tt), te in shared.items():
        back = shared.get((tt, te, qt))
        if back == qe:
            uf.union((qt, qe), (tt, te))
            n_edges += 1
        elif not cfg.get("require_reciprocal", True):
            uf.union((qt, qe), (tt, te))
            n_edges += 1
        else:
            n_oneway += 1
    log(f"orthology edges: {n_edges // 2 if cfg.get('require_reciprocal', True) else n_edges} "
        f"linked, {n_oneway} one-way calls dropped")

    members = defaultdict(list)
    for t in order:
        for r in tabs[t].itertuples():
            members[uf.find((t, r.eve_id))].append((t, r.eve_id))

    # absence observations attach to the group of the querying EVE
    empt = defaultdict(set)
    for r in obs[obs.status == "empty"].itertuples():
        empt[uf.find((r.query_taxon, r.query_eve))].add(r.target_taxon)

    meta = {}
    for t in order:
        for r in tabs[t].itertuples():
            meta[(t, r.eve_id)] = r

    rows = []
    for i, (root, mem) in enumerate(
            sorted(members.items(), key=lambda kv: (-len({m[0] for m in kv[1]}), str(kv[0]))), 1):
        present = {m[0] for m in mem}
        absent = {t for t in empt.get(root, set()) if t not in present}
        fams = Counter(str(meta[m].family) for m in mem)
        rep = meta[mem[0]]
        st = {t: ("1" if t in present else ("0" if t in absent else "?")) for t in order}
        rows.append(dict(
            locus=f"EVEL{i:04d}", family=fams.most_common(1)[0][0],
            families_in_group=";".join(sorted(fams)) if len(fams) > 1 else "",
            virus=str(getattr(rep, "virus", "")), n_taxa_present=len(present),
            n_taxa_absent=len(absent), n_taxa_missing=len(order) - len(present) - len(absent),
            paralogous=int(len(mem) > len(present)),
            # no other taxon could be scored at all: usually a fragmented EVE call whose
            # flanks would not place. Harmless (contributes no signal) but not a character.
            unlinked=int(len(present) == 1 and not absent),
            representative=f"{mem[0][0]}:{rep.contig}:{rep.start}-{rep.end}",
            taxa_present=",".join(sorted(present)), taxa_absent=",".join(sorted(absent)),
            members=",".join(f"{t}:{e}" for t, e in sorted(mem)),
            pattern="".join(st[t] for t in order),
            **{f"st_{t}": st[t] for t in order}))
    import pandas as pd
    g = pd.DataFrame(rows)
    if not g.empty:
        g["informative"] = ((g.n_taxa_present >= 2) & (g.n_taxa_absent >= 1) &
                            (g.n_taxa_present < len(order))).astype(int)
        log(f"{len(g)} orthologous insertion loci; "
            f"{int(g.informative.sum())} parsimony-informative, "
            f"{int((g.n_taxa_present == 1).sum())} taxon-specific, "
            f"{int(g.unlinked.sum())} unlinked (no other taxon scorable)")
    return g


# ----------------------------------------------------------------------------------------
# 7. character matrices
# ----------------------------------------------------------------------------------------

def write_matrices(cfg, order, groups, outdir):
    chars = groups[groups.n_taxa_present >= 1].reset_index(drop=True) if not groups.empty \
        else groups
    mat = {t: "".join(chars[f"st_{t}"]) for t in order} if not chars.empty \
        else {t: "" for t in order}
    n = len(chars)

    with open(os.path.join(outdir, "eve_matrix.tsv"), "w") as fh:
        fh.write("taxon\t" + "\t".join(chars.locus) + "\n")
        for t in order:
            fh.write(t + "\t" + "\t".join(chars[f"st_{t}"]) + "\n")

    with open(os.path.join(outdir, "eve_matrix.fasta"), "w") as fh:
        for t in order:
            fh.write(f">{t}\n{mat[t]}\n")

    pad = max((len(t) for t in order), default=1) + 2
    labels = " ".join(chars.locus) if n else ""
    nex = os.path.join(outdir, "eve_matrix.nex")
    with open(nex, "w") as fh:
        fh.write("#NEXUS\n\n[ Vilina: orthologous EVE insertions as binary characters.\n"
                 "  1 = insertion present, 0 = orthologous site present and empty,\n"
                 "  ? = site not determined (flanks absent, or spanning a gap large enough\n"
                 "      to hide an insertion) -- NOT evidence of absence. ]\n\n")
        fh.write(f"BEGIN TAXA;\n  DIMENSIONS NTAX={len(order)};\n  TAXLABELS "
                 + " ".join(order) + ";\nEND;\n\n")
        fh.write(f"BEGIN CHARACTERS;\n  DIMENSIONS NCHAR={max(n, 1)};\n"
                 "  FORMAT DATATYPE=STANDARD SYMBOLS=\"01\" MISSING=? GAP=-;\n")
        if n:
            fh.write("  CHARLABELS " + labels + ";\n")
        fh.write("  MATRIX\n")
        for t in order:
            fh.write(f"    {t.ljust(pad)}{mat[t] or '?'}\n")
        fh.write("  ;\nEND;\n\n")
        fh.write("BEGIN PAUP;\n"
                 "  [ Dollo: an insertion at a given site arises once and is never regained,\n"
                 "    but may be lost by deletion. ]\n"
                 "  ctype dollo.up:all;\n"
                 "  set criterion=parsimony maxtrees=1000 increase=auto;\n"
                 "  [ hsearch; bootstrap nreps=1000 search=heuristic; ]\n"
                 "END;\n")

    with open(os.path.join(outdir, "eve_matrix_mrbayes.nex"), "w") as fh:
        fh.write(f"#NEXUS\n\nBEGIN DATA;\n  DIMENSIONS NTAX={len(order)} NCHAR={max(n, 1)};\n"
                 "  FORMAT DATATYPE=RESTRICTION MISSING=? GAP=-;\n  MATRIX\n")
        for t in order:
            fh.write(f"    {t.ljust(pad)}{mat[t] or '?'}\n")
        fh.write("  ;\nEND;\n\nBEGIN MRBAYES;\n"
                 "  [ only variable sites are scored, so correct for the missing constant ones ]\n"
                 "  lset coding=noabsencesites;\n"
                 "  [ mcmc ngen=1000000 samplefreq=1000; sump; sumt; ]\n"
                 "END;\n")
    log(f"character matrix: {len(order)} taxa x {n} loci")
    return chars


# ----------------------------------------------------------------------------------------
# 8. splits
# ----------------------------------------------------------------------------------------

def compatible(a, b, order):
    """Four-gamete test on two binary characters, ignoring '?' in either."""
    seen = set()
    for t in order:
        x, y = a[t], b[t]
        if x == "?" or y == "?":
            continue
        seen.add((x, y))
    return len(seen) < 4


def split_support(cfg, order, chars, outdir):
    import pandas as pd
    if chars.empty:
        pd.DataFrame().to_csv(os.path.join(outdir, "eve_splits.tsv"), sep="\t", index=False)
        return pd.DataFrame()
    # An insertion in {T4,T5} and one in {T1,T2,T3} are synapomorphies of opposite sides of
    # the SAME branch, so they are pooled: the key is the unordered bipartition. Polarity
    # (which side carries the insertion) is kept separately, since for a gain-only character
    # it names the derived clade.
    byp = defaultdict(list)
    for r in chars.itertuples():
        st = {t: getattr(r, f"st_{t}") for t in order}
        inn = frozenset(t for t in order if st[t] == "1")
        out = frozenset(t for t in order if st[t] == "0")
        if len(inn) < 2 or len(out) < 1:
            continue
        byp[frozenset((inn, out))].append((r.locus, st, inn))
    rows, states = [], []
    for key, items in sorted(byp.items(), key=lambda kv: (-len(kv[1]), sorted(map(sorted, kv[0])))):
        sides = sorted(key, key=lambda s: (sorted(s) != sorted(items[0][2]), sorted(s)))
        a, b = sides[0], sides[1] if len(sides) > 1 else frozenset()
        pol = Counter(",".join(sorted(i[2])) for i in items)
        rows.append(dict(branch=",".join(sorted(a)) + " | " + ",".join(sorted(b)),
                         side_a=",".join(sorted(a)), side_b=",".join(sorted(b)),
                         n_insertions=len(items),
                         derived_clades="; ".join(f"{k} ({v})" for k, v in pol.most_common()),
                         loci=",".join(i[0] for i in items)))
        states.append(items[0][1])
    sp = pd.DataFrame(rows)
    # greedy compatible set, strongest split first
    accepted, status = [], []
    for st in states:
        if all(compatible(st, acc, order) for acc in accepted):
            accepted.append(st)
            status.append("compatible")
        else:
            status.append("conflicting")
    sp["set"] = status
    sp.to_csv(os.path.join(outdir, "eve_splits.tsv"), sep="\t", index=False)
    log(f"{len(sp)} distinct splits; "
        f"{int((sp.set == 'conflicting').sum())} conflict with a better-supported one")
    return sp


def pairwise_counts(order, chars, outdir):
    import pandas as pd
    import numpy as np
    n = len(order)
    sh = np.zeros((n, n), int)
    comp = np.zeros((n, n), int)
    for r in chars.itertuples():
        st = {t: getattr(r, f"st_{t}") for t in order}
        for i, a in enumerate(order):
            for j, b in enumerate(order):
                if st[a] == "1" and st[b] == "1":
                    sh[i, j] += 1
                if st[a] != "?" and st[b] != "?":
                    comp[i, j] += 1
    pd.DataFrame(sh, index=order, columns=order).to_csv(
        os.path.join(outdir, "pairwise_shared.tsv"), sep="\t")
    return sh, comp


# ----------------------------------------------------------------------------------------
# 9. plots
# ----------------------------------------------------------------------------------------

def genome_layout(cfg, taxon):
    """Cumulative coordinates: contigs longest-first."""
    fa = V.Faidx(taxon_genome(cfg, taxon))
    cs = sorted(((c, fa.length(c)) for c in fa.order), key=lambda x: -x[1])
    off, cum = {}, 0
    for c, L in cs:
        off[c] = cum
        cum += L
    return dict(offset=off, total=cum, contigs=cs)


def family_palette(groups):
    fams = [f for f, _ in Counter(groups.family).most_common()]
    cmap = {f: V.CAT[i] for i, f in enumerate(fams[:8])}
    for f in fams[8:]:
        cmap[f] = V.OTHER
    return cmap, fams[:8]


def _legend_families(ax, cmap, top, extra=None):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=cmap[f], markersize=7,
                markeredgecolor=V.SURFACE, label=f) for f in top]
    if extra:
        h += extra
    ax.legend(handles=h, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1),
              title="viral family", title_fontsize=7)


def plot_presence_matrix(order, chars, outdir, cmap):
    import numpy as np
    plt = V._style()
    if chars.empty:
        return
    from matplotlib.colors import ListedColormap
    ch = chars.sort_values(["pattern", "family"]).reset_index(drop=True)
    code = {"1": 2, "0": 1, "?": 0}
    M = np.array([[code[getattr(r, f"st_{t}")] for r in ch.itertuples()] for t in order])
    w = max(5.5, min(14, 3.0 + 0.055 * len(ch)))
    fig, ax = plt.subplots(figsize=(w, 0.42 * len(order) + 1.8))
    ax.imshow(M, aspect="auto", interpolation="nearest",
              cmap=ListedColormap(["#ffffff", "#e4e3df", "#104281"]), vmin=0, vmax=2)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=9)
    ax.set_xticks([])
    ax.set_xlabel(f"{len(ch)} orthologous EVE insertion loci "
                  f"({int(ch.informative.sum())} parsimony-informative)")
    ax.grid(False)
    ax.set_title("EVE insertion presence / absence")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor="#104281", label="insertion present (1)"),
                       Patch(facecolor="#e4e3df", label="orthologous site empty (0)"),
                       Patch(facecolor="#ffffff", edgecolor=V.GRID, label="not determined (?)")],
              fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
    V._save(fig, outdir, "C01_presence_matrix")


def plot_split_support(sp, outdir, n_taxa):
    plt = V._style()
    if sp.empty:
        return
    d = sp.head(25).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 0.34 * len(d) + 1.8))
    colors = [V.CAT[0] if s == "compatible" else V.CAT[1] for s in d.set]
    ax.barh(range(len(d)), d.n_insertions, color=colors, height=0.62)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d.branch, fontsize=7)
    for i, v in enumerate(d.n_insertions):
        ax.text(v + max(d.n_insertions) * 0.01, i, str(v), va="center", fontsize=7, color=V.INK2)
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlim(0, max(d.n_insertions) * 1.12 + 0.3)
    ax.set_xlabel("orthologous EVE insertions supporting the split")
    ax.grid(axis="y", visible=False)
    ax.set_title("Phylogenetic signal: insertions per bipartition")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=V.CAT[0], label="compatible set"),
                       Patch(facecolor=V.CAT[1], label="conflicts with a stronger split")],
              fontsize=7, loc="lower right")
    V._save(fig, outdir, "C02_split_support")


def plot_shared_heatmap(order, sh, outdir):
    import numpy as np
    plt = V._style()
    n = len(order)
    M = sh.astype(float).copy()
    np.fill_diagonal(M, np.nan)
    fig, ax = plt.subplots(figsize=(1.0 * n + 2.4, 0.9 * n + 1.8))
    im = ax.imshow(M, cmap="Blues", interpolation="nearest")
    ax.set_xticks(range(n))
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n))
    ax.set_yticklabels(order, fontsize=8)
    mx = np.nanmax(M) if np.isfinite(M).any() else 1
    for i in range(n):
        for j in range(n):
            if i != j:
                ax.text(j, i, int(sh[i, j]), ha="center", va="center", fontsize=8,
                        color=V.SURFACE if M[i, j] > 0.6 * mx else V.INK)
    ax.grid(False)
    ax.set_title("Shared orthologous EVE insertions")
    fig.colorbar(im, ax=ax, shrink=0.75, label="shared insertions")
    V._save(fig, outdir, "C03_shared_insertions")


def plot_dotplot(cfg, a, b, chars, obs, tabs, layouts, cmap, top, outdir):
    plt = V._style()
    la, lb = layouts[a], layouts[b]
    pos = {}
    for t in (a, b):
        for r in tabs[t].itertuples():
            pos[(t, r.eve_id)] = layouts[t]["offset"][r.contig] + (r.start + r.end) / 2
    xs, ys, cs = [], [], []
    for r in chars.itertuples():
        if getattr(r, f"st_{a}") != "1" or getattr(r, f"st_{b}") != "1":
            continue
        ma = [m for m in r.members.split(",") if m.startswith(a + ":")]
        mb = [m for m in r.members.split(",") if m.startswith(b + ":")]
        if not ma or not mb:
            continue
        ka, kb = (a, ma[0].split(":", 1)[1]), (b, mb[0].split(":", 1)[1])
        if ka in pos and kb in pos:
            xs.append(pos[ka] / 1e6)
            ys.append(pos[kb] / 1e6)
            cs.append(cmap.get(r.family, V.OTHER))
    ex, ey = [], []
    for r in obs[(obs.query_taxon == a) & (obs.target_taxon == b) &
                 (obs.status == "empty")].itertuples():
        k = (a, r.query_eve)
        if k in pos and r.contig in lb["offset"]:
            ex.append(pos[k] / 1e6)
            ey.append((lb["offset"][r.contig] + (int(r.site_start) + int(r.site_end)) / 2) / 1e6)
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    for c, L in la["contigs"][:40]:
        if la["offset"][c]:
            ax.axvline(la["offset"][c] / 1e6, color=V.GRID, lw=0.6, zorder=1)
    for c, L in lb["contigs"][:40]:
        if lb["offset"][c]:
            ax.axhline(lb["offset"][c] / 1e6, color=V.GRID, lw=0.6, zorder=1)
    if ex:
        ax.scatter(ex, ey, s=26, facecolors="none", edgecolors=V.OTHER, linewidths=0.9,
                   zorder=2, label=f"site empty in {b} (n={len(ex)})")
    if xs:
        ax.scatter(xs, ys, s=46, c=cs, edgecolors=V.SURFACE, linewidths=0.8, zorder=3)
    ax.set_xlabel(f"{a} — cumulative position (Mb)")
    ax.set_ylabel(f"{b} — cumulative position (Mb)")
    ax.set_xlim(0, la["total"] / 1e6)
    ax.set_ylim(0, lb["total"] / 1e6)
    ax.set_title(f"Orthologous EVEs: {a} vs {b} (n={len(xs)} shared)")
    from matplotlib.lines import Line2D
    extra = [Line2D([], [], marker="o", ls="", markerfacecolor="none",
                    markeredgecolor=V.OTHER, markersize=6, label="site empty in target")] \
        if ex else None
    _legend_families(ax, cmap, [f for f in top if f in set(chars.family)], extra)
    V._save(fig, outdir, f"C04_dotplot_{a}_vs_{b}")


def _ribbon(ax, x0, y0, x1, y1, w, color, alpha=0.55):
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    ym = (y0 + y1) / 2
    verts = [(x0 - w, y0), (x0 - w, ym), (x1 - w, ym), (x1 - w, y1),
             (x1 + w, y1), (x1 + w, ym), (x0 + w, ym), (x0 + w, y0), (x0 - w, y0)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=color, edgecolor="none",
                           alpha=alpha, zorder=2))


def plot_ribbon(cfg, order, chars, tabs, layouts, cmap, top, outdir):
    plt = V._style()
    from matplotlib.patches import Rectangle
    if chars.empty or len(order) < 2:
        return
    norm = str(cfg.get("ribbon_scale", "proportional")).startswith("norm")
    Lmax = max(layouts[t]["total"] for t in order)
    width = {t: 1.0 if norm else layouts[t]["total"] / Lmax for t in order}
    pos = {}
    for t in order:
        for r in tabs[t].itertuples():
            pos[(t, r.eve_id)] = ((layouts[t]["offset"][r.contig] + (r.start + r.end) / 2)
                                  / layouts[t]["total"] * width[t], int(r.length_nt))
    y = {t: len(order) - 1 - i for i, t in enumerate(order)}
    th = 0.13
    fig, ax = plt.subplots(figsize=(11, 1.15 * len(order) + 1.4))
    for t in order:
        ax.add_patch(Rectangle((0, y[t] - th / 2), width[t], th, facecolor=V.GRID,
                               edgecolor="none", zorder=3))
        for c, L in layouts[t]["contigs"][:60]:
            o = layouts[t]["offset"][c]
            if o:
                x = o / layouts[t]["total"] * width[t]
                ax.plot([x, x], [y[t] - th / 2, y[t] + th / 2], color=V.SURFACE, lw=0.9, zorder=4)
        ax.text(-0.012, y[t], t, ha="right", va="center", fontsize=10)
        ax.text(width[t] + 0.008, y[t], f"{layouts[t]['total'] / 1e6:.0f} Mb · "
                f"{len(tabs[t])} EVEs", ha="left", va="center", fontsize=7, color=V.INK2)
    n_rib = 0
    for a, b in zip(order, order[1:]):
        for r in chars.itertuples():
            if getattr(r, f"st_{a}") != "1" or getattr(r, f"st_{b}") != "1":
                continue
            ma = [m for m in r.members.split(",") if m.startswith(a + ":")]
            mb = [m for m in r.members.split(",") if m.startswith(b + ":")]
            if not ma or not mb:
                continue
            ka, kb = (a, ma[0].split(":", 1)[1]), (b, mb[0].split(":", 1)[1])
            if ka not in pos or kb not in pos:
                continue
            xa, la_ = pos[ka]
            xb, lb_ = pos[kb]
            w = max(0.0022, (la_ + lb_) / 2 / Lmax)
            _ribbon(ax, xa, y[a] - th / 2, xb, y[b] + th / 2, w, cmap.get(r.family, V.OTHER))
            n_rib += 1
    ax.set_xlim(-0.16, 1.22)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.axis("off")
    ax.set_title("Shared EVE insertions between adjacent taxa"
                 + ("  (tracks normalised)" if norm else "  (tracks to scale)"),
                 loc="left")
    _legend_families(ax, cmap, [f for f in top if f in set(chars.family)])
    log(f"ribbon plot: {n_rib} ribbons between adjacent taxa")
    V._save(fig, outdir, "C05_ribbon")


# ----------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------

def write_report(cfg, order, groups, chars, sp, outdir, tree=None):
    import base64
    pdir = os.path.join(outdir, "plots")
    imgs = []
    for root, _, fs in os.walk(pdir):
        for f in sorted(fs):
            if f.endswith(".png"):
                imgs.append(os.path.join(root, f))
    imgs.sort(key=lambda p: os.path.basename(p))
    figs = "".join(
        f'<figure><img src="data:image/png;base64,'
        f'{base64.b64encode(open(p, "rb").read()).decode()}"/>'
        f'<figcaption>{os.path.basename(p)[:-4]}</figcaption></figure>' for p in imgs)
    card = lambda k, v: f'<div class="c"><div class="v">{v}</div><div class="k">{k}</div></div>'
    cards = (card("taxa", len(order)) + card("orthologous loci", len(chars)) +
             card("informative", int(chars.informative.sum()) if not chars.empty else 0) +
             card("splits", len(sp)) +
             card("conflicting", int((sp.set == "conflicting").sum()) if not sp.empty else 0))
    if tree and tree.get("built"):
        cards += (card("Dollo losses", tree["dollo_losses"]) +
                  card("optimal trees", tree["n_optimal"]))
    tree_html = ""
    if tree and tree.get("built"):
        tree_html = (f"<h2>Cladogram</h2><p>Rooted on <b>{', '.join(tree['outgroup'])}</b>. "
                     f"Dollo parsimony over {tree['n_characters']} insertion characters "
                     f"({tree['n_informative']} informative); best score "
                     f"<b>{tree['dollo_losses']}</b> loss(es) from a {tree['method']} search of "
                     f"{tree['topologies_examined']:,} topologies; {tree['n_optimal']} optimal "
                     f"tree(s); {tree['bootstrap_reps']} bootstrap replicates.</p>"
                     f"<pre style=\"font-size:12px;white-space:pre-wrap\">{tree['newick']}</pre>")
    tbl = sp.head(40).to_html(index=False, border=0, classes="t") if not sp.empty else "<p>none</p>"
    html = f"""<!doctype html><html><head><meta charset="utf-8"><title>Vilina comparative</title>
<style>body{{font:14px/1.5 system-ui,sans-serif;background:#fcfcfb;color:#0b0b0b;max-width:1200px;margin:24px auto;padding:0 16px}}
h1{{font-size:22px}}.sub{{color:#52514e;font-size:13px;margin:-8px 0 16px}}.cs{{display:flex;flex-wrap:wrap;gap:12px}}
.c{{border:1px solid #e4e3df;border-radius:8px;padding:10px 14px;min-width:120px}}.v{{font-size:22px;font-weight:600}}
.k{{color:#52514e;font-size:12px}}figure{{margin:16px 0}}img{{max-width:100%}}figcaption{{color:#52514e;font-size:12px}}
.t{{border-collapse:collapse;font-size:12px}}.t td,.t th{{border-bottom:1px solid #e4e3df;padding:4px 6px;text-align:left}}
.wrap{{overflow-x:auto}}</style></head><body>
<h1>Vilina — EVEs as rare genomic changes</h1>
<p class="sub">विलीन · {len(order)} taxa: {', '.join(order)}</p>
<div class="cs">{cards}</div>
<h2>Figures</h2>{figs}
{tree_html}
<h2>Splits</h2><div class="wrap">{tbl}</div>
<p style="color:#52514e;font-size:12px">A character is only informative where at least one taxon
scores 0 — an orthologous site located by both flanks and demonstrably empty. '?' means the site
could not be placed and is <em>not</em> evidence of absence. Run the matrix under Dollo parsimony
(<code>eve_matrix.nex</code>) rather than standard parsimony: insertions are effectively
irreversible gains, losses are deletions.</p></body></html>"""
    p = os.path.join(outdir, "compare_report.html")
    open(p, "w").write(html)
    log(f"report: {p}")


# ----------------------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------------------

def cmd_compare(a):
    import pandas as pd
    cfg = load_compare_cfg(a.config)
    for k in ("genomes_dir", "compare_outdir", "flank_mapper", "minimap2_preset"):
        if getattr(a, k.replace("compare_", "") if k == "compare_outdir" else k, None):
            cfg[k] = getattr(a, k.replace("compare_", "") if k == "compare_outdir" else k)
    if a.jobs:
        cfg["compare_jobs"] = a.jobs
    out = cfg["compare_outdir"]
    os.makedirs(os.path.join(out, "plots", "dotplots"), exist_ok=True)
    order, genomes = discover_genomes(cfg)

    if not a.skip_per_taxon:
        run_per_taxon(cfg, order, genomes, a.config, force=a.force)
    tabs = load_taxon_eves(cfg, order)
    if sum(len(v) for v in tabs.values()) == 0:
        raise SystemExit("no EVEs in any taxon — nothing to compare")

    flank_paths, flank_lens = extract_flanks(cfg, order, tabs)
    hits = map_flanks(cfg, order, flank_paths)
    obs = score_sites(cfg, order, tabs, hits, flank_lens)
    obs.to_csv(os.path.join(out, "ortho_observations.tsv"), sep="\t", index=False)

    groups = build_orthogroups(cfg, order, tabs, obs)
    groups.to_csv(os.path.join(out, "eve_orthogroups.tsv"), sep="\t", index=False)
    chars = write_matrices(cfg, order, groups, out)
    sp = split_support(cfg, order, chars, out)
    sh, comp = pairwise_counts(order, chars, out)

    layouts = {t: genome_layout(cfg, t) for t in order}
    cmap, top = family_palette(groups) if not groups.empty else ({}, [])
    pdir = os.path.join(out, "plots")
    plot_presence_matrix(order, chars, pdir, cmap)
    plot_split_support(sp, pdir, len(order))
    plot_shared_heatmap(order, sh, pdir)
    plot_ribbon(cfg, order, chars, tabs, layouts, cmap, top, pdir)

    want = cfg.get("dotplot_pairs", "all")
    if want == "adjacent":
        pairs = list(zip(order, order[1:]))
    elif isinstance(want, (list, tuple)):
        pairs = [(safe_name(x[0]), safe_name(x[1])) for x in want]
    else:
        pairs = list(itertools.combinations(order, 2))
    if len(pairs) > int(cfg["max_dotplot_pairs"]):
        log(f"{len(pairs)} pairs exceeds max_dotplot_pairs; drawing adjacent pairs only")
        pairs = list(zip(order, order[1:]))
    for x, y in pairs:
        plot_dotplot(cfg, x, y, chars, obs, tabs, layouts, cmap, top,
                     os.path.join(pdir, "dotplots"))

    tree = VT.build(cfg, order, chars, out, pdir)

    summ = dict(taxa=order, n_loci=int(len(chars)),
                n_loci_scorable=int((groups.unlinked == 0).sum()) if not groups.empty else 0,
                n_informative=int(chars.informative.sum()) if not chars.empty else 0,
                n_splits=int(len(sp)),
                n_conflicting=int((sp.set == "conflicting").sum()) if not sp.empty else 0,
                eves_per_taxon={t: int(len(tabs[t])) for t in order},
                site_calls=obs.status.value_counts().to_dict() if not obs.empty else {},
                tree=tree)
    json.dump(summ, open(os.path.join(out, "compare_summary.json"), "w"), indent=2)
    write_report(cfg, order, groups, chars, sp, out, tree)
    log(json.dumps({k: v for k, v in summ.items() if k != "taxa"}))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True)
    p.add_argument("--genomes-dir", dest="genomes_dir")
    p.add_argument("--outdir", dest="compare_outdir")
    p.add_argument("--flank-mapper", dest="flank_mapper", choices=["minimap2", "blastn"])
    p.add_argument("--minimap2-preset", dest="minimap2_preset")
    p.add_argument("--jobs", type=int, help="per-taxon pipelines to run at once")
    p.add_argument("--skip-per-taxon", action="store_true",
                   help="reuse existing per_taxon/<taxon>/eves.tsv")
    p.add_argument("--force", action="store_true")
    cmd_compare(p.parse_args())


if __name__ == "__main__":
    main()

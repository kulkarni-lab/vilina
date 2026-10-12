#!/usr/bin/env python3
"""Check the comparative module and the cladogram against the planted events and true tree.

  python test/check_compare_truth.py test/compare_data test/results_cmp
"""
import os
import re
import sys
import pandas as pd

data, res = sys.argv[1], sys.argv[2]
truth = pd.read_csv(f"{data}/truth_events.tsv", sep="\t")
g = pd.read_csv(f"{res}/eve_orthogroups.tsv", sep="\t", dtype=str)
sp = pd.read_csv(f"{res}/eve_splits.tsv", sep="\t")
taxa = sorted({t for v in truth.present_in for t in v.split(",")})
ok = True

# ---- 1. every planted event recovered as exactly one locus -----------------------------
print(f"{'event':6} {'expected taxa':22} {'locus':10} status")
recovered = {}
for r in g.itertuples():
    recovered.setdefault(r.taxa_present, []).append(r.locus)
used = set()
for r in truth.itertuples():
    want = ",".join(sorted(r.present_in.split(",")))
    hits = [x for x in recovered.get(want, []) if x not in used]
    status = "OK" if hits else "MISSING"
    if hits:
        used.add(hits[0])
    ok &= status == "OK"
    print(f"{r.event:6} {want:22} {(hits[0] if hits else '-'):10} {status}")

# ---- 2. extra loci: only tolerated when nothing else could be scored --------------------
extra = g[~g.locus.isin(used)]
bad_extra = extra[extra.unlinked != "1"] if "unlinked" in extra else extra
for r in extra.itertuples():
    kind = "unlinked fragment (harmless)" if getattr(r, "unlinked", "0") == "1" \
        else "UNEXPECTED (possible false homology)"
    print(f"\nextra locus {r.locus} [{r.taxa_present}]: {kind}")
ok &= len(bad_extra) == 0

# ---- 3. convergent same-family insertions must not merge -------------------------------
e7 = g[g.taxa_present.isin(["T1", "T4"]) & (g.family == "Chuviridae")]
print(f"\nconvergent same-family insertions kept separate: "
      f"{'OK' if len(e7) == 2 else 'FAILED'} ({len(e7)} loci, expected 2)")
ok &= len(e7) == 2

# ---- 4. cladogram topology (compared UNROOTED) -----------------------------------------
def parse_newick(s):
    s = re.sub(r"\)[^,()]*", ")", s.strip().rstrip(";"))      # drop support labels
    pos = [0]

    def node():
        if s[pos[0]] == "(":
            pos[0] += 1
            kids = [node()]
            while s[pos[0]] == ",":
                pos[0] += 1
                kids.append(node())
            pos[0] += 1                                       # ')'
            return tuple(kids)
        m = re.match(r"[^,()]+", s[pos[0]:])
        pos[0] += len(m.group(0))
        return m.group(0).strip()
    return node()


def leaves(t):
    if isinstance(t, str):
        return [t]
    out = []
    for x in t:
        out.extend(leaves(x))
    return out


def splits(t, allt):
    out = set()

    def rec(x):
        ls = frozenset(leaves(x))
        if not isinstance(x, str):
            if 1 < len(ls) < len(allt) - 1:
                out.add(frozenset({ls, frozenset(allt) - ls}))
            for k in x:
                rec(k)
        return ls
    rec(t)
    return out


tree_file = f"{res}/eve_tree.nwk"
if not os.path.exists(tree_file):
    print("\nno cladogram written — FAILED")
    ok = False
else:
    got = parse_newick(open(tree_file).read())
    true = parse_newick(open(f"{data}/true_tree.nwk").read())
    gs, ts = splits(got, taxa), splits(true, taxa)
    same = gs == ts
    ok &= same
    print(f"\ncladogram: {open(tree_file).read().strip()}")
    print(f"true tree: {open(f'{data}/true_tree.nwk').read().strip()}")
    print(f"unrooted topology matches true tree: {'OK' if same else 'FAILED'}")
    if not same:
        show = lambda S: "; ".join(sorted("|".join(sorted(",".join(sorted(h)) for h in b))
                                          for b in S))
        print(f"  recovered splits: {show(gs)}")
        print(f"  true splits:      {show(ts)}")

print("\nALL CHECKS PASSED" if ok else "\nFAILURES ABOVE")
sys.exit(0 if ok else 1)

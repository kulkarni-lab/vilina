#!/usr/bin/env python3
"""
vilina_tree.py — Dollo parsimony cladogram from orthologous EVE insertions.

An integration at a given site is a gain that happens once and can afterwards only be lost
by deletion, so the right optimality criterion is **Dollo parsimony**: every character's
1-state arises exactly once, at the MRCA of the taxa carrying it, and the tree is scored by
the number of secondary losses it forces. Standard (Fitch) parsimony would instead be free
to invent the same insertion twice, which is precisely the event the method assumes cannot
happen.

Missing data ('?') is treated as a wildcard resolved in favour of fewest losses, so an
undetermined site never manufactures a deletion.

Search: exhaustive over all topologies while that is affordable, otherwise stepwise
addition + SPR hill-climbing from several random addition orders. Outgroup monophyly is
enforced; the root is placed on the outgroup branch, which is what gives Dollo its polarity.

Runs inside `vilina.py compare`, or standalone on any binary matrix:

    python scripts/vilina_tree.py --groups eve_orthogroups.tsv --outdir . \
        --taxa T1,T2,T3,T4,T5 --outgroup T5 --bootstrap 200
"""
import argparse
import os
import random
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vilina as V  # noqa: E402

TDEFAULTS = {
    "build_tree": True,
    "outgroup": None,            # taxon name or list; default = last taxon in taxon_order
    "bootstrap": 100,
    "exhaustive_max_trees": 200_000,
    "tree_search_restarts": 12,
    "max_optimal_trees": 200,
}


# ----------------------------------------------------------------------------------------
# tree representation: nested 2-tuples, leaves are taxon strings
# ----------------------------------------------------------------------------------------

def leaves(t):
    """Works for binary search trees and for multifurcating consensus trees."""
    if isinstance(t, str):
        return [t]
    out = []
    for x in t:
        out.extend(leaves(x))
    return out


def canon(t):
    """Canonical form so topologically identical trees compare equal."""
    if isinstance(t, str):
        return t
    a, b = canon(t[0]), canon(t[1])
    return (a, b) if min(leaves(a)) <= min(leaves(b)) else (b, a)


def index_tree(t):
    """Post-order node list: (leafset, left, right); leaves have left = -1."""
    nodes = []

    def rec(x):
        if isinstance(x, str):
            nodes.append((frozenset([x]), -1, -1))
            return len(nodes) - 1
        li, ri = rec(x[0]), rec(x[1])
        nodes.append((nodes[li][0] | nodes[ri][0], li, ri))
        return len(nodes) - 1

    return nodes, rec(t)


def double_factorial_trees(m):
    """Number of rooted binary trees on m labelled taxa = (2m-3)!!"""
    n = 1
    for k in range(3, 2 * m - 2, 2):
        n *= k
    return n if m > 1 else 1


def all_rooted_trees(taxa):
    if len(taxa) == 1:
        yield taxa[0]
        return
    first, rest = taxa[0], taxa[1:]
    for t in all_rooted_trees(rest):
        yield from insert_everywhere(t, first)


def insert_everywhere(tree, sub):
    yield (sub, tree)
    if not isinstance(tree, str):
        a, b = tree
        for na in insert_everywhere(a, sub):
            yield (na, b)
        for nb in insert_everywhere(b, sub):
            yield (a, nb)


def remove_subtree(tree, target):
    if tree == target:
        return None
    if isinstance(tree, str):
        return tree
    a, b = tree
    ra = remove_subtree(a, target)
    if ra is None:
        return b
    rb = remove_subtree(b, target)
    if rb is None:
        return a
    return (ra, rb)


def subtrees(tree):
    yield tree
    if not isinstance(tree, str):
        yield from subtrees(tree[0])
        yield from subtrees(tree[1])


def spr_neighbours(tree):
    seen = {canon(tree)}
    for st in list(subtrees(tree)):
        if st == tree:
            continue
        rest = remove_subtree(tree, st)
        if rest is None:
            continue
        for nt in insert_everywhere(rest, st):
            c = canon(nt)
            if c not in seen:
                seen.add(c)
                yield nt


# ----------------------------------------------------------------------------------------
# Dollo parsimony
# ----------------------------------------------------------------------------------------

def char_losses(nodes, root, ones, zeros):
    """Losses forced by this character, and the node where the single gain is placed.
    '?' taxa are neither in `ones` nor `zeros`, so they never force a loss."""
    if not ones:
        return 0, None
    m = root
    while True:
        li, ri = nodes[m][1], nodes[m][2]
        if li < 0:
            break
        if ones <= nodes[li][0]:
            m = li
        elif ones <= nodes[ri][0]:
            m = ri
        else:
            break
    if len(ones) == 1:
        return 0, m

    def rec(i):
        ls, li, ri = nodes[i]
        if ls & ones:
            return 0 if li < 0 else rec(li) + rec(ri)
        return 1 if (ls & zeros) else 0

    return rec(m), m


def tree_cost(tree, chars):
    nodes, root = index_tree(tree)
    return sum(char_losses(nodes, root, o, z)[0] for o, z in chars)


# ----------------------------------------------------------------------------------------
# search
# ----------------------------------------------------------------------------------------

def _compose(og_tree, ig_tree):
    return (og_tree, ig_tree)


def stepwise_addition(chars, ingroup, og_tree, rng):
    order = list(ingroup)
    rng.shuffle(order)
    if len(order) == 1:
        return order[0]
    tree = (order[0], order[1])
    for taxon in order[2:]:
        best, bc = None, None
        for cand in insert_everywhere(tree, taxon):
            c = tree_cost(_compose(og_tree, cand), chars)
            if bc is None or c < bc:
                best, bc = cand, c
        tree = best
    return tree


def hill_climb(chars, ig_tree, og_tree):
    cur = ig_tree
    cost = tree_cost(_compose(og_tree, cur), chars)
    improved = True
    while improved:
        improved = False
        for nb in spr_neighbours(cur):
            c = tree_cost(_compose(og_tree, nb), chars)
            if c < cost:
                cur, cost, improved = nb, c, True
                break
    return cur, cost


def search(chars, ingroup, outgroup, cfg, quick=False):
    """Returns (best_cost, [optimal trees]) with outgroup monophyly enforced."""
    og_trees = [outgroup[0]] if len(outgroup) == 1 else list(all_rooted_trees(list(outgroup)))
    n_total = double_factorial_trees(len(ingroup)) * len(og_trees)
    cap = int(cfg.get("exhaustive_max_trees", 200_000))
    best, opt = None, []
    maxkeep = int(cfg.get("max_optimal_trees", 200))

    if n_total <= (3000 if quick else cap):
        for og in og_trees:
            for ig in all_rooted_trees(list(ingroup)):
                t = _compose(og, ig)
                c = tree_cost(t, chars)
                if best is None or c < best:
                    best, opt = c, [canon(t)]
                elif c == best and len(opt) < maxkeep:
                    ct = canon(t)
                    if ct not in opt:
                        opt.append(ct)
        return best, opt, "exhaustive", n_total

    rng = random.Random(20261011)
    restarts = 3 if quick else int(cfg.get("tree_search_restarts", 12))
    for og in og_trees:
        for _ in range(restarts):
            ig = stepwise_addition(chars, ingroup, og, rng)
            ig, c = hill_climb(chars, ig, og)
            t = canon(_compose(og, ig))
            if best is None or c < best:
                best, opt = c, [t]
            elif c == best and t not in opt and len(opt) < maxkeep:
                opt.append(t)
    return best, opt, "heuristic", n_total


# ----------------------------------------------------------------------------------------
# consensus, support, newick
# ----------------------------------------------------------------------------------------

def clades(tree):
    out = set()

    def rec(t):
        ls = frozenset(leaves(t))
        if not isinstance(t, str):
            out.add(ls)
            rec(t[0])
            rec(t[1])
        return ls

    rec(tree)
    return out


def strict_consensus(trees, taxa):
    """Clade set shared by every optimal tree, rebuilt into a (possibly multifurcating) tree."""
    common = set.intersection(*[clades(t) for t in trees])
    common = {c for c in common if 1 < len(c) < len(taxa)}
    nodes = sorted(common, key=len)
    children = {c: [] for c in nodes}
    assigned = set()
    for i, c in enumerate(nodes):
        for d in nodes[:i]:
            if d < c and not any(d < e < c for e in nodes[:i]):
                children[c].append(d)
                assigned.add(d)
    tops = [c for c in nodes if c not in assigned]

    def build(c):
        kids = [build(k) for k in children[c]]
        covered = set().union(*children[c]) if children[c] else set()
        kids += sorted(c - covered)
        return tuple(kids) if len(kids) != 1 else kids[0]

    roots = [build(c) for c in tops]
    roots += sorted(set(taxa) - set().union(*tops) if tops else set(taxa))
    return tuple(roots) if len(roots) != 1 else roots[0]


def to_newick(tree, labels=None):
    def rec(t):
        if isinstance(t, str):
            return t
        inner = ",".join(rec(x) for x in t)
        lab = str(labels.get(frozenset(leaves(t)), "") or "") if labels else ""
        return f"({inner}){lab}"
    return rec(tree) + ";"


def branch_support(tree, chars, loci, boot=None):
    """Per clade: insertions whose single gain maps there, and losses on that branch."""
    nodes, root = index_tree(tree)
    gains = defaultdict(list)
    losses = Counter()
    for (o, z), loc in zip(chars, loci):
        n_loss, m = char_losses(nodes, root, o, z)
        if m is not None:
            gains[nodes[m][0]].append(loc)
        losses[nodes[m][0] if m is not None else None] += n_loss
    rows = []
    for ls, _l, _r in nodes:
        if len(ls) < 2:
            continue
        rows.append(dict(clade=",".join(sorted(ls)), n_taxa=len(ls),
                         n_insertions=len(gains.get(ls, [])),
                         insertions=",".join(gains.get(ls, [])),
                         losses_below=losses.get(ls, 0),
                         bootstrap=(round(100 * boot.get(ls, 0.0)) if boot else "")))
    return sorted(rows, key=lambda r: (-r["n_taxa"], r["clade"])), gains


def bootstrap_clades(chars, loci, ingroup, outgroup, cfg, reps):
    rng = random.Random(7)
    n = len(chars)
    if n == 0 or reps <= 0:
        return {}
    cnt = Counter()
    for _ in range(reps):
        idx = [rng.randrange(n) for _ in range(n)]
        rs = [chars[i] for i in idx]
        _c, opt, _m, _nt = search(rs, ingroup, outgroup, cfg, quick=True)
        cl = set.intersection(*[clades(t) for t in opt]) if len(opt) > 1 else clades(opt[0])
        for c in cl:
            cnt[c] += 1
    return {c: v / reps for c, v in cnt.items()}


# ----------------------------------------------------------------------------------------
# cladogram
# ----------------------------------------------------------------------------------------

def plot_cladogram(tree, rows, outdir, title, method, cost, n_chars, name="C06_cladogram"):
    plt = V._style()
    sup = {r["clade"]: r for r in rows}
    tips = []

    def height(t):
        return 0 if isinstance(t, str) else 1 + max(height(x) for x in _kids(t))

    def _kids(t):
        return list(t) if isinstance(t, tuple) else []

    H = height(tree)
    ypos = {}
    counter = [0]

    def place(t):
        if isinstance(t, str):
            y = counter[0]
            counter[0] += 1
            tips.append(t)
            ypos[("leaf", t)] = y
            return y
        kids = sorted(_kids(t), key=lambda k: (len(leaves(k)), min(leaves(k))))
        ys = [place(k) for k in kids]
        y = sum(ys) / len(ys)
        ypos[("node", frozenset(leaves(t)))] = y
        return y

    place(tree)
    nleaf = counter[0]
    fig, ax = plt.subplots(figsize=(7.2 + 0.06 * max(len(t) for t in tips),
                                    0.46 * nleaf + 1.8))

    def xof(t):
        return H if isinstance(t, str) else H - height(t)

    def yof(t):
        return ypos[("leaf", t)] if isinstance(t, str) else ypos[("node", frozenset(leaves(t)))]

    def draw(t):
        if isinstance(t, str):
            return
        kids = sorted(_kids(t), key=lambda k: (len(leaves(k)), min(leaves(k))))
        x0 = xof(t)
        ax.plot([x0, x0], [min(yof(k) for k in kids), max(yof(k) for k in kids)],
                color=V.INK, lw=1.6, solid_capstyle="round", zorder=3)
        for k in kids:
            ax.plot([x0, xof(k)], [yof(k), yof(k)], color=V.INK, lw=1.6,
                    solid_capstyle="round", zorder=3)
            r = sup.get(",".join(sorted(leaves(k))))
            if r and r["n_taxa"] >= 2:
                xm = (x0 + xof(k)) / 2
                if r["n_insertions"]:
                    ax.text(xm, yof(k) + 0.16, f"{r['n_insertions']}", ha="center",
                            va="bottom", fontsize=9, color=V.CAT[0], fontweight="bold")
                bs = r.get("bootstrap", "")
                if bs != "":
                    ax.text(xm, yof(k) - 0.17, f"{bs}", ha="center", va="top", fontsize=7.5,
                            color=V.INK2 if bs >= 70 else V.OTHER)
            draw(k)

    draw(tree)
    for t in tips:
        ax.text(H + 0.06, yof(t), t, va="center", fontsize=10)
    ax.set_xlim(-0.25, H + 0.6 + 0.04 * max(len(t) for t in tips))
    ax.set_ylim(-0.95, nleaf - 0.45)
    ax.axis("off")
    ax.set_title(title, loc="left")
    ax.text(0, -0.72, f"Dollo parsimony · {n_chars} insertion characters · {cost} loss"
                      f"{'es' if cost != 1 else ''} · {method} search\n"
                      f"bold = insertions supporting the branch · grey = bootstrap %",
            fontsize=7.5, color=V.INK2, va="top")
    V._save(fig, outdir, name)


# ----------------------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------------------

def build(cfg, order, chars_df, outdir, plots_dir):
    """Returns a dict summary; writes newick, support table and the cladogram."""
    import pandas as pd
    res = dict(built=False)
    if not cfg.get("build_tree", True):
        return res
    if len(order) < 3:
        V.log("need >= 3 taxa for a cladogram")
        return res

    og = cfg.get("outgroup") or [order[-1]]
    og = [og] if isinstance(og, str) else list(og)
    og = [str(x) for x in og]
    bad = [t for t in og if t not in order]
    if bad:
        raise SystemExit(f"outgroup not among the taxa: {', '.join(bad)}")
    ingroup = [t for t in order if t not in og]
    if len(ingroup) < 2:
        V.log("outgroup leaves fewer than 2 ingroup taxa; skipping tree")
        return res
    if not cfg.get("outgroup"):
        V.log(f"no `outgroup:` set — rooting on '{og[0]}' (last in taxon_order); "
              f"Dollo polarity depends on this")

    chars, loci = [], []
    for r in chars_df.itertuples():
        ones = frozenset(t for t in order if getattr(r, f"st_{t}") == "1")
        zeros = frozenset(t for t in order if getattr(r, f"st_{t}") == "0")
        if not ones:
            continue
        chars.append((ones, zeros))
        loci.append(r.locus)
    inf = sum(1 for o, z in chars if 2 <= len(o) < len(order) and z)
    if inf == 0:
        V.log("WARNING: no parsimony-informative insertions — the cladogram is unresolved "
              "by the EVE data alone")

    cost, opt, method, n_total = search(chars, ingroup, og, cfg)
    V.log(f"Dollo search ({method}, {n_total:,} topologies): best = {cost} loss(es), "
          f"{len(opt)} optimal tree(s)")

    reps = int(cfg.get("bootstrap", 100))
    boot = bootstrap_clades(chars, loci, ingroup, og, cfg, reps) if reps else {}
    best = opt[0]
    cons = strict_consensus(opt, order) if len(opt) > 1 else best
    rows, gains = branch_support(best, chars, loci, boot)

    pd.DataFrame(rows).to_csv(os.path.join(outdir, "eve_tree_support.tsv"),
                              sep="\t", index=False)
    labels = {frozenset(r["clade"].split(",")): (r["bootstrap"] if r["bootstrap"] != "" else "")
              for r in rows}
    with open(os.path.join(outdir, "eve_tree.nwk"), "w") as fh:
        fh.write(to_newick(best, labels) + "\n")
    with open(os.path.join(outdir, "eve_tree_all_optimal.nwk"), "w") as fh:
        for t in opt:
            fh.write(to_newick(t) + "\n")
    if len(opt) > 1:
        with open(os.path.join(outdir, "eve_tree_consensus.nwk"), "w") as fh:
            fh.write(to_newick(cons) + "\n")

    title = "EVE insertion cladogram" + (f"  (strict consensus of {len(opt)} optimal trees)"
                                         if len(opt) > 1 else "")
    plot_cladogram(cons if len(opt) > 1 else best, rows, plots_dir, title, method, cost,
                   len(chars))
    res.update(built=True, dollo_losses=int(cost), n_optimal=len(opt), method=method,
               outgroup=og, newick=to_newick(best, labels),
               topologies_examined=int(n_total), bootstrap_reps=reps,
               n_characters=len(chars), n_informative=inf)
    return res


def main():
    import pandas as pd
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--groups", required=True, help="eve_orthogroups.tsv")
    p.add_argument("--taxa", required=True, help="comma-separated, in order")
    p.add_argument("--outgroup", help="comma-separated; default = last taxon")
    p.add_argument("--outdir", required=True)
    p.add_argument("--plots")
    p.add_argument("--bootstrap", type=int, default=100)
    a = p.parse_args()
    order = [x.strip() for x in a.taxa.split(",")]
    df = pd.read_csv(a.groups, sep="\t", dtype=str)
    cfg = dict(TDEFAULTS)
    cfg["bootstrap"] = a.bootstrap
    if a.outgroup:
        cfg["outgroup"] = [x.strip() for x in a.outgroup.split(",")]
    os.makedirs(a.outdir, exist_ok=True)
    plots = a.plots or os.path.join(a.outdir, "plots")
    os.makedirs(plots, exist_ok=True)
    r = build(cfg, order, df, a.outdir, plots)
    print(r.get("newick", "(no tree)"))


if __name__ == "__main__":
    main()

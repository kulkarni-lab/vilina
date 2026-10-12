#!/usr/bin/env python3
"""Compare pipeline calls with planted EVEs: python test/check_truth.py test/data/truth.tsv test/results/eves.tsv"""
import sys
import pandas as pd

t = pd.read_csv(sys.argv[1], sep="\t")
e = pd.read_csv(sys.argv[2], sep="\t")
rows, ok = [], 0
for r in t.itertuples():
    m = e[(e.contig == r.contig) & (e.start <= r.end) & (e.end >= r.start)]
    excluded = r.note.startswith("EXCLUDED")
    if m.empty:
        status = "correctly excluded" if excluded else "MISSED"
        rows.append((r.id, status, "", "", "", r.note))
    else:
        x = m.iloc[0]
        status = "FALSE POSITIVE (should be excluded)" if excluded else "recovered"
        rows.append((r.id, status, x.eve_id, x.completeness, x.functional_class, r.note))
    ok += status in ("recovered", "correctly excluded")
out = pd.DataFrame(rows, columns=["planted", "status", "eve_id", "completeness", "functional_class", "expected"])
print(out.to_string(index=False))
print(f"\n{ok}/{len(t)} planted elements handled correctly; "
      f"{len(e) - out.eve_id.replace('', pd.NA).dropna().nunique()} unexpected calls")
sys.exit(0 if ok == len(t) else 1)

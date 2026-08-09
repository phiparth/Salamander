#!/usr/bin/env python3
"""Transfer stage 1 - build the yeast pairs that strip2 is fine-tuned on.

Yeast Tm in this table is ORGANISM-level: every S. cerevisiae protein carries the same
value, so yeast-vs-yeast pairs have zero dTm and teach nothing. Instead we build
CROSS-ORGANISM pairs where the featurised (low-Tm) member is a yeast protein and the partner
is a hotter homolog from the same cluster - literally "how would this yeast protein look if
it came from something that lives hotter".

Output format is identical to build_pairs.py, so the same embed/train stages consume it.

    python training/build_yeast_pairs.py --csv training/data/combined_clustered_proteins.csv \
           --work_dir runs/yeast --taxid 559292
"""
import argparse, itertools, json, os
from collections import Counter

AA20 = "ACDEFGHIKLMNPQRSTVWY"
AA2IDX = {a: i for i, a in enumerate(AA20)}


def clean_seq(s):
    s = "".join(str(s).split()).upper()
    return "".join(("X" if c in "UZOB" else c) for c in s if c.isalpha())


def seq_hash(seq):
    import hashlib
    return hashlib.sha1(seq.encode()).hexdigest()[:16]


def main():
    p = argparse.ArgumentParser(description="yeast->hotter transfer pairs")
    p.add_argument("--csv", required=True)
    p.add_argument("--work_dir", required=True)
    p.add_argument("--taxid", default="559292", help="559292 = S. cerevisiae")
    p.add_argument("--cluster_col", default="Cluster ID")
    p.add_argument("--top_n_pairs", type=int, default=3)
    p.add_argument("--min_dtm", type=float, default=1.0)
    p.add_argument("--cons_thresh", type=float, default=0.90)
    a = p.parse_args()

    import pandas as pd
    from pyfamsa import Aligner, Sequence

    df = pd.read_csv(a.csv).dropna(subset=["Sequence", "Tmax_C"]).copy()
    df["Sequence"] = df["Sequence"].map(clean_seq)
    df["Organism_TaxID"] = df["Organism_TaxID"].astype(str)

    aligner = Aligner(guide_tree="upgma")
    records, unique, used = [], set(), 0

    for cid, sub in df.groupby(a.cluster_col):
        sub = sub.drop_duplicates("Sequence")
        if len(sub) < 2:
            continue
        members = list(sub.itertuples(index=False))
        if not any(str(m.Organism_TaxID) == a.taxid for m in members):
            continue                                    # cluster must contain yeast

        named = [("m%d" % i, m.Sequence) for i, m in enumerate(members)]
        tax = {"m%d" % i: str(m.Organism_TaxID) for i, m in enumerate(members)}
        tm = {"m%d" % i: float(m.Tmax_C) for i, m in enumerate(members)}
        raw = {"m%d" % i: m.Sequence for i, m in enumerate(members)}
        aln = aligner.align([Sequence(n.encode(), s.encode()) for n, s in named])
        amap = {r.id.decode(): r.sequence.decode() for r in aln}
        rows = list(amap.values())
        L = len(rows[0])

        conserved = [False] * L
        for c in range(L):
            col = [r[c] for r in rows if r[c] != "-"]
            if not col:
                continue
            _, n = Counter(col).most_common(1)[0]
            if n / len(col) >= a.cons_thresh:
                conserved[c] = True

        cand = []
        for x, y in itertools.combinations(amap, 2):
            lo, hi = (x, y) if tm[x] < tm[y] else (y, x)
            if tax[lo] != a.taxid:                      # the featurised member must be yeast
                continue
            d = tm[hi] - tm[lo]
            if d < a.min_dtm:
                continue
            cand.append((d, lo, hi))
        cand = sorted(cand)[::-1][:a.top_n_pairs]
        if cand:
            used += 1

        for dtm, lo, hi in cand:
            lo_aln, hi_aln = amap[lo], amap[hi]
            lo_pos, targets = -1, []
            for c in range(L):
                lr = lo_aln[c]
                if lr != "-":
                    lo_pos += 1
                if conserved[c] or lr == "-":
                    continue
                hr = hi_aln[c]
                if hr == "-" or hr not in AA2IDX:
                    continue
                targets.append([lo_pos, AA2IDX[hr], c])
            if not targets:
                continue
            records.append({"cluster": str(cid), "dtm": dtm, "lo_tm": tm[lo], "hi_tm": tm[hi],
                            "lo_seq": raw[lo], "lo_hash": seq_hash(raw[lo]),
                            "aln_len": L, "targets": targets})
            unique.add(raw[lo])

    os.makedirs(a.work_dir, exist_ok=True)
    json.dump(records, open(os.path.join(a.work_dir, "pairs.json"), "w"))
    with open(os.path.join(a.work_dir, "unique_seqs.txt"), "w") as fh:
        fh.write("\n".join(sorted(unique)))

    print("clusters with yeast->hotter pairs : %d" % used)
    print("pair records                      : %d" % len(records))
    print("unique low-Tm sequences (yeast)   : %d" % len(unique))
    print("training targets                  : %d" % sum(len(r["targets"]) for r in records))
    if records:
        d = sorted(r["dtm"] for r in records)
        print("dTm min/median/max                : %.0f / %.0f / %.0f"
              % (d[0], d[len(d) // 2], d[-1]))
    print("\nwrote %s/pairs.json and unique_seqs.txt" % a.work_dir)


if __name__ == "__main__":
    main()

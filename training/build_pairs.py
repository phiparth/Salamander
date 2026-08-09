#!/usr/bin/env python3
"""Stage 1 of training - turn a clustered, Tm-annotated protein table into training pairs.

The model never sees a measured ddG or dTm. Its supervision comes from evolution: within a
cluster of homologs, take a cool protein and a hotter one, and at every non-conserved
position the label is the residue the HOT protein has there. The cool protein is the one
that gets embedded; the hot one only supplies labels and is never fed to the network.

Per cluster:
  1. one FAMSA MSA over ALL members, so column indices mean the same thing for every pair
  2. a column is conserved if its modal non-gap residue covers >= --cons_thresh
  3. enumerate member pairs, require |dTm| >= --min_dtm, keep the --top_n_pairs largest
  4. emit a target at each non-conserved column where BOTH members have a residue

    python training/build_pairs.py --csv training/data/combined_clustered_proteins.csv \
           --work_dir runs/prok --exclude_taxid 559292

--exclude_taxid 559292 drops S. cerevisiae, which is how the prokaryote-only ("prok") model
was built. The distributed strip2 head is this prok model fine-tuned on yeast afterwards.
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
    p = argparse.ArgumentParser(description="build pairs.json for pairwise training")
    p.add_argument("--csv", required=True)
    p.add_argument("--work_dir", required=True)
    p.add_argument("--cluster_col", default="Cluster ID")
    p.add_argument("--cons_thresh", type=float, default=0.90,
                   help="conserved columns are excluded from the labels")
    p.add_argument("--min_dtm", type=float, default=1.0,
                   help="a pair needs at least this Tm gap to be used")
    p.add_argument("--top_n_pairs", type=int, default=3,
                   help="pairs kept per cluster, largest dTm first")
    p.add_argument("--exclude_taxid", default=None,
                   help="drop this organism entirely (559292 = S. cerevisiae)")
    p.add_argument("--max_clusters", type=int, default=0, help="0 = all; smaller for smoke tests")
    a = p.parse_args()

    import pandas as pd
    from pyfamsa import Aligner, Sequence

    df = pd.read_csv(a.csv).dropna(subset=["Sequence", "Tmax_C"]).copy()
    df["Sequence"] = df["Sequence"].map(clean_seq)
    if a.exclude_taxid:
        before = len(df)
        df = df[df["Organism_TaxID"].astype(str) != str(a.exclude_taxid)]
        print("excluded taxid %s: %d -> %d rows" % (a.exclude_taxid, before, len(df)))

    cids = sorted(df[a.cluster_col].unique())
    if a.max_clusters:
        cids = cids[:a.max_clusters]
    print("clusters: %d" % len(cids), flush=True)

    aligner = Aligner(guide_tree="upgma")
    records, unique = [], set()
    used = 0
    for k, cid in enumerate(cids, 1):
        sub = df[df[a.cluster_col] == cid].drop_duplicates("Sequence")
        if len(sub) < 2:
            continue
        named = [("m%d" % i, r.Sequence) for i, r in sub.reset_index().iterrows()]
        aln = aligner.align([Sequence(n.encode(), s.encode()) for n, s in named])
        amap = {r.id.decode(): r.sequence.decode() for r in aln}
        rows = list(amap.values())
        L = len(rows[0])
        tm = {"m%d" % i: float(r.Tmax_C) for i, r in sub.reset_index().iterrows()}
        raw = {"m%d" % i: r.Sequence for i, r in sub.reset_index().iterrows()}

        conserved = [False] * L
        for c in range(L):
            col = [r[c] for r in rows if r[c] != "-"]
            if not col:
                continue
            _, n = Counter(col).most_common(1)[0]
            if n / len(col) >= a.cons_thresh:
                conserved[c] = True

        pairs = sorted((abs(tm[x] - tm[y]), x, y)
                       for x, y in itertools.combinations(amap, 2))
        pairs = [q for q in pairs if q[0] >= a.min_dtm][::-1][:a.top_n_pairs]
        if pairs:
            used += 1

        for dtm, x, y in pairs:
            lo, hi = (x, y) if tm[x] < tm[y] else (y, x)
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
                    continue                       # missing target, not a label
                targets.append([lo_pos, AA2IDX[hr], c])
            if not targets:
                continue
            records.append({"cluster": str(cid), "dtm": dtm, "lo_tm": tm[lo], "hi_tm": tm[hi],
                            "lo_seq": raw[lo], "lo_hash": seq_hash(raw[lo]),
                            "aln_len": L, "targets": targets})
            unique.add(raw[lo])
        if k % 500 == 0:
            print("  %d/%d clusters, %d records" % (k, len(cids), len(records)), flush=True)

    os.makedirs(a.work_dir, exist_ok=True)
    json.dump(records, open(os.path.join(a.work_dir, "pairs.json"), "w"))
    with open(os.path.join(a.work_dir, "unique_seqs.txt"), "w") as fh:
        fh.write("\n".join(sorted(unique)))

    ntar = sum(len(r["targets"]) for r in records)
    print("\nclusters contributing pairs : %d" % used)
    print("pair records                : %d" % len(records))
    print("unique low-Tm sequences     : %d  (these are what stage 2 embeds)" % len(unique))
    print("training targets            : %d" % ntar)
    if records:
        d = sorted(r["dtm"] for r in records)
        print("dTm min/median/max          : %.0f / %.0f / %.0f"
              % (d[0], d[len(d) // 2], d[-1]))
    print("\nwrote %s/pairs.json and unique_seqs.txt" % a.work_dir)


if __name__ == "__main__":
    main()

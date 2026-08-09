#!/usr/bin/env python3
"""Step 2 - run strip2 over the mutable positions and apply the set1* confidence filter.

At every position that step 1 left mutable, strip2 emits a 20-way softmax. A mutation is
proposed when the top-1 residue differs from wild type. set1* then keeps only the confident
ones:

    P(mutant) > P_MIN            default 0.10   - the model must actually commit
    P(mutant) >= RATIO * P(WT)   default 1.20   - and must clearly prefer it to what is there

The ratio term is the important half. An absolute probability cut alone keeps mutations the
model rates barely above wild type; requiring a margin over P(WT) removes those.

    python pipeline/step2_design_set1star.py --work work/era \
           --model models/head_strip2.pt --p-min 0.10 --wt-ratio 1.20
"""
import argparse, json, os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from strip2.core import embed, load_head, residue_probs, AA20, AA2IDX  # noqa: E402


def main():
    p = argparse.ArgumentParser(description="strip2 inference + set1* shortlist")
    p.add_argument("--work", required=True, help="directory from step 1")
    p.add_argument("--model", default="models/head_strip2.pt")
    p.add_argument("--plm", help="ProtT5 path or HF id (else $PLM, else the HF default)")
    p.add_argument("--p-min", type=float, default=0.10, help="minimum P(mutant)")
    p.add_argument("--wt-ratio", type=float, default=1.20,
                   help="P(mutant) must be at least this multiple of P(wild type)")
    p.add_argument("--emb-cache", default=None, help="npz to cache the embedding in")
    a = p.parse_args()

    cons = json.load(open(os.path.join(a.work, "conserved.json")))
    seq, L = cons["query"], cons["L"]
    frozen = set(cons["frozen_0based"])
    print("query %s  L=%d  frozen=%d  mutable=%d"
          % (cons["query_name"], L, len(frozen), L - len(frozen)), flush=True)

    cache = a.emb_cache or os.path.join(a.work, "embedding.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        emb, pooled = z["emb"], z["pooled"]
        print("reused embedding from %s" % cache, flush=True)
    else:
        emb, pooled = embed(seq, a.plm)
        np.savez_compressed(cache, emb=emb, pooled=pooled)
        print("embedded and cached -> %s" % cache, flush=True)
    assert emb.shape[0] == L, "embedding length %d != sequence length %d" % (emb.shape[0], L)

    net, meta = load_head(a.model)
    print("model %s  in_dim=%s  hidden=%s" % (os.path.basename(a.model), meta.get("in_dim"),
                                              meta.get("hidden_dims") or meta.get("hidden")),
          flush=True)

    mutable = [r for r in range(L) if r not in frozen]
    probs = residue_probs(net, meta, emb, pooled, positions=mutable)

    baseline, kept = [], []
    for r in mutable:
        pr = probs[r]
        i1 = int(np.argmax(pr))
        a1, wt = AA20[i1], seq[r]
        if a1 == wt:
            continue
        p1 = float(pr[i1])
        pwt = float(pr[AA2IDX[wt]]) if wt in AA2IDX else 0.0
        rec = {"residue": r + 1, "wt": wt, "mut": a1,
               "p_mut": round(p1, 6), "p_wt": round(pwt, 6),
               "ratio": round(p1 / pwt, 3) if pwt > 0 else None}
        baseline.append(rec)
        if p1 > a.p_min and (pwt == 0.0 or p1 >= a.wt_ratio * pwt):
            kept.append(rec)

    design = list(seq)
    for m in kept:
        design[m["residue"] - 1] = m["mut"]
    out = {"query_name": cons["query_name"], "query": seq, "L": L,
           "p_min": a.p_min, "wt_ratio": a.wt_ratio,
           "n_frozen": len(frozen), "n_baseline": len(baseline), "n_set1star": len(kept),
           "baseline_mutations": baseline, "set1star_mutations": kept,
           "set1star_sequence": "".join(design)}
    with open(os.path.join(a.work, "set1star.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("\nconf=0 top-1 mutations : %d" % len(baseline))
    print("surviving set1*        : %d" % len(kept))
    if kept:
        print("\n%-9s %-6s %9s %9s %8s" % ("residue", "mut", "P(mut)", "P(wt)", "ratio"))
        for m in kept:
            print("%-9d %s->%s   %9.4f %9.4f %8s"
                  % (m["residue"], m["wt"], m["mut"], m["p_mut"], m["p_wt"],
                     "inf" if m["ratio"] is None else "%.2f" % m["ratio"]))
    print("\nwrote %s/set1star.json" % a.work)


if __name__ == "__main__":
    main()

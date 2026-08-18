#!/usr/bin/env python3
"""Run all four steps end to end.

    python pipeline/run_pipeline.py --query examples/era_lbd.fasta --pdb structures/era.pdb \
           --out work/era --cons-thresh 0.80 --wt-ratio 1.20 --ddg-cut 0.5 --repair

Each step writes its own JSON into --out, so you can stop after any step, change one
parameter and rerun only what follows. Embeddings and FoldX results are cached.
"""
import argparse, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def run(script, args):
    cmd = [sys.executable, os.path.join(HERE, script)] + args
    print("\n" + "-" * 70)
    print("$ " + " ".join(cmd))
    print("-" * 70, flush=True)
    r = subprocess.run(cmd)
    if r.returncode != 0:
        sys.exit("%s failed with exit code %d" % (script, r.returncode))


def skip(n, script, args, produces, resume):
    """Run a step, or report it already done.

    run_pipeline drives every step, including the ortholog search - so re-running it after
    a later step failed repeats step 1, and repeats whatever was wrong with it. --resume
    reuses finished outputs instead.
    """
    if resume and os.path.exists(produces):
        print("\nstep %d: %s already exists, skipping (--resume)"
              % (n, os.path.basename(produces)), flush=True)
        return
    run(script, args)


def main():
    p = argparse.ArgumentParser(description="strip2 thermostable design, all four steps")
    p.add_argument("--query", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--pdb", help="structure for the FoldX steps; omitted = stop after step 2")
    p.add_argument("--model", default=os.path.join(ROOT, "models", "head_strip2.pt"))
    p.add_argument("--plm")
    # step 1
    p.add_argument("--source", default="orthodb", choices=["orthodb", "blast", "local"])
    p.add_argument("--family")
    p.add_argument("--n-orthologs", type=int, default=10)
    p.add_argument("--cons-thresh", type=float, default=0.80)
    p.add_argument("--min-identity", type=float, default=0.30)
    p.add_argument("--length-band", type=float, default=2.0)
    p.add_argument("--extra-freeze", default="")
    p.add_argument("--email")
    p.add_argument("--gene", help="gene name; avoids OrthoDB /blast (see step 1)")
    p.add_argument("--og", help="OrthoDB orthogroup id; avoids OrthoDB /blast")
    p.add_argument("--alignment", help="a finished MSA; skips the ortholog search AND the "
                                      "aligner (see step 1)")
    p.add_argument("--aligner", default="auto", choices=["auto", "famsa", "biopython"])
    p.add_argument("--resume", action="store_true",
                   help="skip steps whose output already exists in --out")
    # step 2
    p.add_argument("--p-min", type=float, default=0.10)
    p.add_argument("--wt-ratio", type=float, default=1.20)
    # steps 3-4
    p.add_argument("--foldx")
    p.add_argument("--chain", default="A")
    p.add_argument("--ddg-cut", type=float, default=0.5)
    p.add_argument("--repair", action="store_true")
    p.add_argument("--jobs", type=int, default=6)
    p.add_argument("--no-combined", action="store_true")
    a = p.parse_args()

    os.makedirs(a.out, exist_ok=True)

    s1 = ["--query", a.query, "--out", a.out, "--source", a.source,
          "--n-orthologs", str(a.n_orthologs), "--cons-thresh", str(a.cons_thresh),
          "--min-identity", str(a.min_identity), "--length-band", str(a.length_band)]
    if a.family:
        s1 += ["--family", a.family]
    if a.email:
        s1 += ["--email", a.email]
    if a.extra_freeze:
        s1 += ["--extra-freeze", a.extra_freeze]
    if a.gene:
        s1 += ["--gene", a.gene]
    if a.og:
        s1 += ["--og", a.og]
    if a.alignment:
        s1 += ["--alignment", a.alignment]
    if a.aligner != "auto":
        s1 += ["--aligner", a.aligner]
    skip(1, "step1_conserved_sites.py", s1, os.path.join(a.out, "conserved.json"), a.resume)

    s2 = ["--work", a.out, "--model", a.model, "--p-min", str(a.p_min),
          "--wt-ratio", str(a.wt_ratio)]
    if a.plm:
        s2 += ["--plm", a.plm]
    skip(2, "step2_design_set1star.py", s2, os.path.join(a.out, "set1star.json"), a.resume)

    if not a.pdb:
        print("\nno --pdb given, stopping after step 2. The set1* design is in "
              "%s/set1star.json" % a.out)
        return

    s3 = ["--work", a.out, "--pdb", a.pdb, "--chain", a.chain,
          "--ddg-cut", str(a.ddg_cut), "--jobs", str(a.jobs)]
    if a.foldx:
        s3 += ["--foldx", a.foldx]
    if a.repair:
        s3 += ["--repair"]
    run("step3_foldx_filter.py", s3)

    s4 = ["--work", a.out]
    if not a.no_combined:
        s4 += ["--combined"]
    if a.foldx:
        s4 += ["--foldx", a.foldx]
    run("step4_final_sequence.py", s4)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Step 4 - assemble the final design and write a report.

Optionally scores the whole surviving set as ONE multi-mutant in FoldX. That number is not
the sum of the singles: mutations interact, and on our benchmark the combined value ran
3 kcal/mol below additive on hAChE and 1.7 above it on TPH1. If you intend to quote a ddG
for the design, quote this one, not the sum.

    python pipeline/step4_final_sequence.py --work work/era --combined
"""
import argparse, glob, json, os, shutil, subprocess, sys, time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from strip2 import foldx as foldx_find  # noqa: E402  (`fx` is a local below)


def dif_values(d):
    f = [x for x in os.listdir(d) if x.startswith("Dif_") and x.endswith(".fxout")]
    if not f:
        return None
    vals = []
    for ln in open(os.path.join(d, f[0])):
        q = ln.rstrip("\n").split("\t")
        if len(q) > 2 and q[0].endswith(".pdb"):
            try:
                vals.append(float(q[1]))
            except ValueError:
                pass
    return vals


def main():
    p = argparse.ArgumentParser(description="assemble the final strip2 design")
    p.add_argument("--work", required=True)
    p.add_argument("--combined", action="store_true",
                   help="also score the full design as one multi-mutant in FoldX")
    p.add_argument("--foldx", help="FoldX executable; else $FOLDX, else found automatically (see find_foldx.py)")
    p.add_argument("--name", default=None, help="name for the output FASTA record")
    a = p.parse_args()

    cons = json.load(open(os.path.join(a.work, "conserved.json")))
    s1 = json.load(open(os.path.join(a.work, "set1star.json")))
    fx_path = os.path.join(a.work, "foldx_filter.json")
    have_fx = os.path.exists(fx_path)
    fx = json.load(open(fx_path)) if have_fx else None

    seq = cons["query"]
    muts = fx["kept"] if have_fx else s1["set1star_mutations"]
    final = list(seq)
    for m in muts:
        final[m["residue"] - 1] = m["mut"]
    final = "".join(final)
    name = a.name or (cons["query_name"] + "_strip2_design")

    combined = None
    if a.combined and have_fx and muts:
        exe = foldx_find.resolve(a.foldx, required=False)
        if not exe:
            print("FoldX not found - skipping the combined score", flush=True)
        else:
            root = os.path.join(a.work, "foldx")
            base = os.path.join(root, "start_Repair.pdb")
            pdbname = "start_Repair.pdb"
            if not os.path.exists(base):
                base, pdbname = os.path.join(root, "start.pdb"), "start.pdb"
            d = os.path.join(root, "combined")
            os.makedirs(d, exist_ok=True)
            shutil.copy(base, os.path.join(d, pdbname))
            # one comma-separated line = a single variant carrying every mutation
            with open(os.path.join(d, "individual_list.txt"), "w", newline="\n") as fh:
                fh.write(",".join("%s%s%d%s" % (m["wt"], fx["chain"], m["residue"], m["mut"])
                                  for m in muts) + ";\n")
            if dif_values(d) is None:
                print("scoring the combined design (%d mutations) ..." % len(muts), flush=True)
                t0 = time.time()
                subprocess.run([exe, "--command=BuildModel", "--pdb=" + pdbname,
                                "--mutant-file=individual_list.txt", "--numberOfRuns=1"],
                               cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                print("  finished in %ds" % (time.time() - t0), flush=True)
            v = dif_values(d)
            if v:
                combined = round(sum(v) / len(v), 2)

    n_mut = len(muts)
    additive = round(sum(m.get("ddg", 0.0) for m in muts), 2) if have_fx else None
    rep = {"name": name, "query": seq, "final_sequence": final, "length": len(seq),
           "n_frozen": cons["n_frozen"], "cons_thresh": cons["cons_thresh"],
           "n_orthologs": cons["n_orthologs"],
           "n_baseline": s1["n_baseline"], "n_set1star": s1["n_set1star"],
           "ddg_cut": fx["ddg_cut"] if have_fx else None,
           "n_final": n_mut, "mutations": muts,
           "sum_single_ddg": additive, "combined_ddg": combined}
    with open(os.path.join(a.work, "final.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    with open(os.path.join(a.work, "final.fasta"), "w") as fh:
        fh.write(">%s\n" % name)
        for i in range(0, len(final), 60):
            fh.write(final[i:i + 60] + "\n")

    bar = "=" * 66
    print("\n" + bar)
    print("%s   %d residues" % (name, len(seq)))
    print(bar)
    print("%-44s %d" % ("orthologs used", cons["n_orthologs"]))
    print("%-44s %d" % ("frozen at conservation >= %.2f" % cons["cons_thresh"], cons["n_frozen"]))
    print("%-44s %d" % ("mutable positions", cons["L"] - cons["n_frozen"]))
    print("%-44s %d" % ("strip2 top-1 mutations", s1["n_baseline"]))
    print("%-44s %d" % ("after set1* (P>%.2f, ratio>=%.2f)" % (s1["p_min"], s1["wt_ratio"]),
                        s1["n_set1star"]))
    if have_fx:
        print("%-44s %d" % ("after FoldX ddG <= %.3f" % fx["ddg_cut"], n_mut))
        print("%-44s %.2f kcal/mol" % ("sum of surviving single ddG", additive))
    if combined is not None:
        print("%-44s %.2f kcal/mol" % ("combined ddG of the whole design", combined))
        print("%-44s %+.2f kcal/mol" % ("epistasis (combined - additive)", combined - additive))
    print(bar)
    if muts:
        print("\n%-9s %-8s %9s %9s" % ("residue", "change", "P(mut)", "ddG"))
        for m in sorted(muts, key=lambda x: x["residue"]):
            print("%-9d %s -> %s   %9.4f %9s"
                  % (m["residue"], m["wt"], m["mut"], m["p_mut"],
                     "%.3f" % m["ddg"] if "ddg" in m else "-"))
    else:
        print("\nno mutations survived the filters - try a lower --cons-thresh, a lower")
        print("--wt-ratio, or a more permissive --ddg-cut")
    print("\nwrote %s/final.fasta and final.json" % a.work)


if __name__ == "__main__":
    main()

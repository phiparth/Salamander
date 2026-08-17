#!/usr/bin/env python3
"""Step 3 - score every shortlisted mutation as an isolated point mutation in FoldX and
drop the destabilising ones.

strip2 is trained on evolutionary signal and knows nothing about structure. It reaches for
prolines because real thermophiles carry them - but not for whether the local backbone can
accommodate one. In our benchmark the same substitution type appeared at both extremes:
Gly->Pro at -2.85 kcal/mol on ERa, Gln->Pro at +4.60 on hDnmt3a. This step is what separates
those two cases.

    python pipeline/step3_foldx_filter.py --work work/era --pdb structures/era.pdb \
           --foldx /path/to/foldx --ddg-cut 0.5 --repair --jobs 6

The PDB must contain the same sequence as the query (verified before anything runs).
--repair runs FoldX RepairPDB first; skip it only if your structure is already repaired.
"""
import argparse, glob, json, os, shutil, subprocess, sys, time, collections

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from strip2 import foldx as fx  # noqa: E402

THREE2ONE = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
             "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
             "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
             "TRP": "W", "TYR": "Y"}


def pdb_residues(path):
    """-> {(chain, resseq): one_letter} from ATOM records, first model only."""
    out = {}
    for ln in open(path, encoding="utf-8", errors="ignore"):
        if ln.startswith("ENDMDL"):
            break
        if ln.startswith("ATOM"):
            rn, ch, num = ln[17:20].strip(), ln[21], ln[22:26].strip()
            if rn in THREE2ONE and num:
                out[(ch, int(num))] = THREE2ONE[rn]
    return out


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


def run_foldx(exe, cwd, args):
    return subprocess.run([exe] + args, cwd=cwd, stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL)


def main():
    p = argparse.ArgumentParser(description="FoldX single-point ddG filter")
    p.add_argument("--work", required=True)
    p.add_argument("--pdb", required=True, help="structure matching the query sequence")
    p.add_argument("--foldx", help="FoldX executable; else $FOLDX, else found automatically (see find_foldx.py)")
    p.add_argument("--chain", default="A")
    p.add_argument("--ddg-cut", type=float, default=0.5,
                   help="keep mutations with single-point ddG <= this (kcal/mol)")
    p.add_argument("--repair", action="store_true", help="run RepairPDB first")
    p.add_argument("--jobs", type=int, default=6, help="parallel FoldX processes")
    a = p.parse_args()

    exe = fx.resolve(a.foldx)
    print("FoldX: %s" % exe, flush=True)
    comp = fx.companions(exe)
    if not comp["rotabase.txt"] and not comp["molecules"]:
        print("  warning: neither rotabase.txt nor molecules/ is beside the executable;"
              " FoldX 4 will not run", flush=True)

    s1 = json.load(open(os.path.join(a.work, "set1star.json")))
    muts = s1["set1star_mutations"]
    seq = s1["query"]
    if not muts:
        sys.exit("no set1* mutations to score - nothing to do")

    # ---- verify the structure matches the design frame BEFORE spending any compute ----
    res = pdb_residues(a.pdb)
    problems = []
    for m in muts:
        key = (a.chain, m["residue"])
        if key not in res:
            problems.append("%s%d absent from chain %s" % (m["wt"], m["residue"], a.chain))
        elif res[key] != m["wt"]:
            problems.append("%s%d is %s in the PDB, design says %s"
                            % (m["wt"], m["residue"], res[key], m["wt"]))
    if problems:
        print("STRUCTURE MISMATCH - the PDB numbering does not match the query sequence:")
        for q in problems[:10]:
            print("   " + q)
        sys.exit("aborting; renumber the PDB or pass the right --chain")
    print("verified all %d mutations against %s chain %s"
          % (len(muts), os.path.basename(a.pdb), a.chain), flush=True)

    root = os.path.join(a.work, "foldx")
    os.makedirs(root, exist_ok=True)
    start = os.path.join(root, "start.pdb")
    shutil.copy(a.pdb, start)

    if a.repair:
        rep = os.path.join(root, "start_Repair.pdb")
        if not os.path.exists(rep):
            print("RepairPDB (this is the slow part) ...", flush=True)
            t0 = time.time()
            run_foldx(exe, root, ["--command=RepairPDB", "--pdb=start.pdb"])
            print("  repaired in %ds" % (time.time() - t0), flush=True)
        if not os.path.exists(rep):
            sys.exit("RepairPDB produced nothing - check the FoldX licence and rotabase")
        base, pdbname = rep, "start_Repair.pdb"
    else:
        base, pdbname = start, "start.pdb"

    # ---- one chunk per worker, each holding several single mutations ----
    nch = max(1, min(a.jobs, len(muts)))
    chunks = [[] for _ in range(nch)]
    for i, m in enumerate(muts):
        chunks[i % nch].append(m)
    dirs = []
    for ci, cm in enumerate(chunks):
        d = os.path.join(root, "chunk_%d" % ci)
        os.makedirs(d, exist_ok=True)
        shutil.copy(base, os.path.join(d, pdbname))
        with open(os.path.join(d, "individual_list.txt"), "w", newline="\n") as fh:
            fh.write("\n".join("%s%s%d%s;" % (m["wt"], a.chain, m["residue"], m["mut"])
                               for m in cm) + "\n")
        dirs.append((d, cm))

    todo = [(d, cm) for d, cm in dirs if dif_values(d) is None]
    if todo:
        print("scoring %d mutations in %d parallel jobs ..." % (len(muts), len(todo)), flush=True)
        t0 = time.time()
        procs = [(subprocess.Popen([exe, "--command=BuildModel", "--pdb=" + pdbname,
                                    "--mutant-file=individual_list.txt", "--numberOfRuns=1"],
                                   cwd=d, stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL), d) for d, _ in todo]
        while any(pr.poll() is None for pr, _ in procs):
            time.sleep(10)
        print("  finished in %ds" % (time.time() - t0), flush=True)

    single = {}
    for d, cm in dirs:
        v = dif_values(d)
        if v is None or len(v) != len(cm):
            sys.exit("FoldX gave %s rows for %s, expected %d"
                     % ("no" if v is None else len(v), d, len(cm)))
        for m, x in zip(cm, v):
            single["%s%d%s" % (m["wt"], m["residue"], m["mut"])] = round(x, 3)

    kept, dropped = [], []
    for m in muts:
        k = "%s%d%s" % (m["wt"], m["residue"], m["mut"])
        m = dict(m, ddg=single[k])
        (kept if single[k] <= a.ddg_cut else dropped).append(m)

    final = list(seq)
    for m in kept:
        final[m["residue"] - 1] = m["mut"]
    out = {"query_name": s1["query_name"], "query": seq, "ddg_cut": a.ddg_cut,
           "chain": a.chain, "pdb": os.path.abspath(a.pdb), "repaired": bool(a.repair),
           "n_in": len(muts), "n_kept": len(kept), "n_dropped": len(dropped),
           "kept": kept, "dropped": dropped, "single_ddg": single,
           "filtered_sequence": "".join(final)}
    with open(os.path.join(a.work, "foldx_filter.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("\n%-9s %-6s %10s %s" % ("residue", "mut", "ddG", "verdict"))
    for m in sorted(muts, key=lambda x: single["%s%d%s" % (x["wt"], x["residue"], x["mut"])]):
        k = single["%s%d%s" % (m["wt"], m["residue"], m["mut"])]
        print("%-9d %s->%s   %10.3f %s" % (m["residue"], m["wt"], m["mut"], k,
                                           "keep" if k <= a.ddg_cut else "drop"))
    print("\nkept %d of %d at ddG <= %.3f kcal/mol" % (len(kept), len(muts), a.ddg_cut))
    print("wrote %s/foldx_filter.json" % a.work)


if __name__ == "__main__":
    main()

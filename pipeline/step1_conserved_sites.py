#!/usr/bin/env python3
"""Step 1 - find orthologs, align them, and freeze the conserved sites.

Positions that are conserved across the ortholog family are almost always doing something
(catalysis, binding, fold nucleation). Freezing them before design is what keeps the active
site intact - in the ERa/hAChE benchmarks PROSS never mutated inside this band either.

    python pipeline/step1_conserved_sites.py --query examples/era_lbd.fasta \
           --out work/era --cons-thresh 0.80 --n-orthologs 10

Ortholog source (--source):
    orthodb   OrthoDB orthologs via the UniProt/OrthoDB REST route (default)
    blast     NCBI blastp, one hit per genus
    local     skip search, use --family (a FASTA you supply)

Everything that changes the frozen set is a flag: --cons-thresh, --n-orthologs, --source,
--min-identity, --length-band.
"""
import argparse, json, os, sys, collections

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from strip2.core import read_fasta, write_fasta, clean_seq, AA20  # noqa: E402


def fetch_orthodb(seq, n, log, min_identity, length_band, og=None, gene=None):
    """OrthoDB orthologs, one sequence per genus."""
    from strip2.orthologs import orthodb_orthologs
    return orthodb_orthologs(seq, n=n, min_identity=min_identity,
                             length_band=length_band, log=log, og=og, gene=gene)


def fetch_blast(seq, n, log, min_identity, length_band, email):
    from strip2.orthologs import blast_orthologs
    return blast_orthologs(seq, n=n, min_identity=min_identity,
                           length_band=length_band, email=email, log=log)


def align(query_name, query, homologs, log):
    """One MSA over query + orthologs. FAMSA, falling back to Biopython if unavailable."""
    try:
        from pyfamsa import Aligner, Sequence
        recs = [Sequence(query_name.encode(), query.encode())]
        recs += [Sequence(n.encode(), s.encode()) for n, s in homologs]
        aln = Aligner(guide_tree="upgma").align(recs)
        return [(r.id.decode(), r.sequence.decode()) for r in aln]
    except ImportError:
        log("pyfamsa not installed - falling back to pairwise alignment")
        from Bio import pairwise2
        rows = [(query_name, query)]
        for nm, s in homologs:
            a = pairwise2.align.globalms(query, s, 2, -1, -10, -0.5, one_alignment_only=True)[0]
            rows.append((nm, a.seqB))
        return rows


def conserved_columns(rows, cons_thresh):
    """A column is conserved if its dominant non-gap residue covers >= cons_thresh."""
    L = len(rows[0][1])
    flags, detail = [], []
    for c in range(L):
        col = [r[1][c] for r in rows if r[1][c] != "-"]
        if not col:
            flags.append(False); detail.append(("-", 0.0, 0)); continue
        modal, k = collections.Counter(col).most_common(1)[0]
        frac = k / len(col)
        flags.append(frac >= cons_thresh)
        detail.append((modal, round(frac, 3), len(col)))
    return flags, detail


def main():
    p = argparse.ArgumentParser(description="freeze conserved sites from an ortholog MSA")
    p.add_argument("--query", required=True, help="FASTA with a single query sequence")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--source", default="orthodb", choices=["orthodb", "blast", "local"])
    p.add_argument("--family", help="FASTA of orthologs (required for --source local)")
    p.add_argument("--n-orthologs", type=int, default=10)
    p.add_argument("--cons-thresh", type=float, default=0.80,
                   help="fraction of non-gap residues the modal residue must reach (0-1)")
    p.add_argument("--min-identity", type=float, default=0.30,
                   help="reject orthologs below this identity to the query")
    p.add_argument("--length-band", type=float, default=2.0,
                   help="keep orthologs whose length is within [1/x, x] times the query")
    p.add_argument("--email", help="contact address for NCBI blastp (--source blast)")
    p.add_argument("--og", help="OrthoDB orthogroup id, e.g. 4385266at2759; skips the "
                               "sequence search, which is the endpoint that breaks")
    p.add_argument("--gene", help="gene name, e.g. ESR1; resolves the family through "
                                 "OrthoDB /search instead of /blast")
    p.add_argument("--extra-freeze", default="",
                   help="space-separated 1-based residue numbers to freeze regardless")
    a = p.parse_args()

    os.makedirs(a.out, exist_ok=True)
    log = lambda *m: print(*m, flush=True)

    recs = read_fasta(a.query)
    if len(recs) != 1:
        sys.exit("--query must contain exactly one sequence, found %d" % len(recs))
    qname, qseq = recs[0]
    log("query %s, %d residues" % (qname, len(qseq)))

    if a.source == "local":
        if not a.family:
            sys.exit("--source local requires --family")
        homologs = [(n, s) for n, s in read_fasta(a.family) if s != qseq]
    elif a.source == "orthodb":
        homologs = fetch_orthodb(qseq, a.n_orthologs, log, a.min_identity, a.length_band,
                                 og=a.og, gene=a.gene)
    else:
        homologs = fetch_blast(qseq, a.n_orthologs, log, a.min_identity, a.length_band, a.email)
    if not homologs:
        sys.exit("no orthologs found.\n"
                 "  If the message above says OrthoDB /blast is not responding, that is an\n"
                 "  outage on their side, not a problem with your install. Use any of:\n"
                 "    --gene <GENE NAME>            resolve the family by name\n"
                 "    --og <ORTHOGROUP ID>          e.g. --og 4385266at2759\n"
                 "    --source blast --email you@institution.edu\n"
                 "    --source local --family my_orthologs.fasta")
    log("using %d orthologs" % len(homologs))

    rows = align(qname, qseq, homologs, log)
    flags, detail = conserved_columns(rows, a.cons_thresh)

    # map alignment columns back onto query residue indices
    qrow = dict(rows)[qname]
    frozen0, qpos = [], -1
    percol = []
    for c, ch in enumerate(qrow):
        if ch == "-":
            continue
        qpos += 1
        modal, frac, depth = detail[c]
        percol.append({"residue": qpos + 1, "aa": ch, "column": c,
                       "modal": modal, "conservation": frac, "depth": depth,
                       "frozen": bool(flags[c])})
        if flags[c]:
            frozen0.append(qpos)

    for tok in a.extra_freeze.split():
        try:
            r = int(tok) - 1
        except ValueError:
            continue
        if 0 <= r < len(qseq) and r not in frozen0:
            frozen0.append(r)
            percol[r]["frozen"] = True
    frozen0 = sorted(set(frozen0))

    write_fasta(os.path.join(a.out, "orthologs.fasta"), homologs)
    write_fasta(os.path.join(a.out, "alignment.fasta"), rows)
    out = {"query_name": qname, "query": qseq, "L": len(qseq),
           "cons_thresh": a.cons_thresh, "n_orthologs": len(homologs),
           "source": a.source, "frozen_0based": frozen0, "n_frozen": len(frozen0),
           "mutable": len(qseq) - len(frozen0), "per_residue": percol}
    with open(os.path.join(a.out, "conserved.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    log("\nconserved at >= %.2f : %d of %d residues frozen, %d mutable"
        % (a.cons_thresh, len(frozen0), len(qseq), len(qseq) - len(frozen0)))
    log("wrote %s/conserved.json, alignment.fasta, orthologs.fasta" % a.out)


if __name__ == "__main__":
    main()

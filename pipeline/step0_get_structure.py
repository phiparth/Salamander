#!/usr/bin/env python3
"""Step 0 (optional) - obtain a structure for your sequence, so you do not have to hunt for
one by hand.

FoldX needs atomic coordinates; it cannot work from a sequence. But the structure does not
have to be experimental - a prediction is fine, and is what the published benchmark used.

Two sources, tried in order:

    alphafold   AlphaFold DB, by UniProt accession. Instant download, already computed,
                numbered 1..L with no gaps, so it always matches a full-length FASTA.
                Requires the protein to be in UniProt.

    esmfold     ESMFold via the public ESM Atlas API. Takes any sequence, no account, no
                GPU, a few seconds. Practical limit is roughly 400 residues; longer
                sequences are usually rejected by the server.

    python pipeline/step0_get_structure.py --query my.fasta --out structures/mine.pdb
    python pipeline/step0_get_structure.py --uniprot P03372 --out structures/era.pdb
    python pipeline/step0_get_structure.py --query my.fasta --source esmfold --out s.pdb

IMPORTANT - numbering. AlphaFold models cover the FULL UniProt sequence. If your query is a
domain (say residues 300-550), the downloaded model is numbered 300-550 in UniProt space
while your FASTA is 1-250. Those do not match and step 3 will refuse to run. Use
--renumber to renumber the output 1..L against your query, or fold the domain with ESMFold
instead, which numbers from 1 by construction.
"""
import argparse, json, os, sys, urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from strip2.core import read_fasta, clean_seq  # noqa: E402

AFDB = "https://alphafold.ebi.ac.uk/api/prediction/%s"
ESM = "https://api.esmatlas.com/foldSequence/v1/pdb/"
UA = {"User-Agent": "Mozilla/5.0 (salamander-structure)"}
THREE2ONE = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
             "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
             "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
             "TRP": "W", "TYR": "Y"}


def get(url, data=None, timeout=300):
    req = urllib.request.Request(url, data=data, headers=UA)
    return urllib.request.urlopen(req, timeout=timeout).read()


def from_alphafold(acc, log):
    log("querying AlphaFold DB for %s ..." % acc)
    meta = json.loads(get(AFDB % acc, timeout=90).decode())
    if not meta:
        raise SystemExit("AlphaFold DB has no entry for %s" % acc)
    m = meta[0]
    log("  %s  (%s, model version %s)" % (m.get("entryId"),
                                          m.get("uniprotDescription", "")[:40],
                                          m.get("latestVersion")))
    return get(m["pdbUrl"], timeout=300).decode()


def from_esmfold(seq, log):
    log("folding %d residues with ESMFold ..." % len(seq))
    if len(seq) > 400:
        log("  WARNING: %d residues; the public API usually refuses more than ~400" % len(seq))
    return get(ESM, data=seq.encode(), timeout=600).decode()


def pdb_seq(text, chain=None):
    """(chain, resseq) -> one letter, first model only."""
    out = {}
    for ln in text.splitlines():
        if ln.startswith("ENDMDL"):
            break
        if ln.startswith("ATOM"):
            rn, ch, num = ln[17:20].strip(), ln[21], ln[22:26].strip()
            if rn in THREE2ONE and num and (chain is None or ch == chain):
                out[(ch, int(num))] = THREE2ONE[rn]
    return out


def renumber(text, offset, keep=None, chain=None):
    """Shift residue numbers by -offset. If `keep` is given, drop residues outside it."""
    out = []
    for ln in text.splitlines():
        if ln.startswith(("ATOM", "HETATM", "TER")) and len(ln) > 26:
            try:
                n = int(ln[22:26])
            except ValueError:
                out.append(ln); continue
            if keep is not None and n not in keep:
                continue
            if chain is not None and ln[21] != chain:
                continue
            out.append(ln[:22] + ("%4d" % (n - offset)) + ln[26:])
        elif ln.startswith("END"):
            out.append(ln)
    out.append("END")
    return "\n".join(out) + "\n"


def extract_query(text, res, chain, seq, log):
    """Locate the query inside a longer structure and cut it out, renumbered 1..L.

    AlphaFold models cover the whole UniProt entry, so a domain query never lines up.
    Finding the query as a substring of the modelled sequence turns a full-length model
    into a usable one.
    """
    nums = sorted(n for c, n in res if c == chain)
    struct_seq = "".join(res[(chain, n)] for n in nums)
    i = struct_seq.find(seq)
    if i < 0:
        return None
    span = nums[i:i + len(seq)]
    log("  found the query at structure residues %d-%d; extracting and renumbering to 1-%d"
        % (span[0], span[-1], len(seq)))
    return renumber(text, span[0] - 1, keep=set(span), chain=chain)


def main():
    p = argparse.ArgumentParser(description="fetch or predict a structure for FoldX")
    p.add_argument("--query", help="FASTA with one sequence (needed for esmfold and checking)")
    p.add_argument("--uniprot", help="UniProt accession, e.g. P03372 (alphafold source)")
    p.add_argument("--out", required=True, help="where to write the .pdb")
    p.add_argument("--source", default="auto", choices=["auto", "alphafold", "esmfold"])
    p.add_argument("--renumber", action="store_true",
                   help="renumber so the first residue is 1 (fixes domain offsets)")
    p.add_argument("--chain", default="A")
    a = p.parse_args()
    log = lambda *m: print(*m, flush=True)

    if not a.query and not a.uniprot:
        sys.exit("give --query (a FASTA) or --uniprot (an accession), or both")

    seq = None
    if a.query:
        recs = read_fasta(a.query)
        if len(recs) != 1:
            sys.exit("--query must hold exactly one sequence, found %d" % len(recs))
        seq = clean_seq(recs[0][1])
        log("query: %s, %d residues" % (recs[0][0], len(seq)))

    src = a.source
    if src == "auto":
        src = "alphafold" if a.uniprot else "esmfold"
    log("source: %s" % src)

    if src == "alphafold":
        if not a.uniprot:
            sys.exit("--source alphafold needs --uniprot")
        text = from_alphafold(a.uniprot, log)
    else:
        if not seq:
            sys.exit("--source esmfold needs --query")
        text = from_esmfold(seq, log)

    res = pdb_seq(text, a.chain)
    if not res:
        sys.exit("no ATOM records for chain %s in the downloaded structure" % a.chain)
    nums = sorted(n for _, n in res)
    log("structure: chain %s, %d residues, numbered %d-%d"
        % (a.chain, len(res), nums[0], nums[-1]))

    # a full-length model against a domain query: cut the query out and renumber it
    if seq and len(res) > len(seq):
        cut = extract_query(text, res, a.chain, seq, log)
        if cut:
            text = cut
            res = pdb_seq(text, a.chain)
            nums = sorted(n for _, n in res)
        else:
            log("  note: the query is not an exact substring of the modelled sequence")

    if a.renumber and nums and nums[0] != 1:
        log("renumbering: %d..%d -> 1..%d" % (nums[0], nums[-1], nums[-1] - nums[0] + 1))
        text = renumber(text, nums[0] - 1)
        res = pdb_seq(text, a.chain)
        nums = sorted(n for _, n in res)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")
    log("wrote %s" % a.out)

    # ---- the check that step 3 will run anyway, reported here while it is cheap to fix ----
    if seq:
        bad = [(i, seq[i - 1], res[(a.chain, i)])
               for i in range(1, len(seq) + 1)
               if (a.chain, i) in res and res[(a.chain, i)] != seq[i - 1]]
        missing = [i for i in range(1, len(seq) + 1) if (a.chain, i) not in res]
        if not bad and not missing:
            log("\nOK: every residue 1-%d matches your query. Ready for step 3." % len(seq))
        else:
            log("\nNUMBERING MISMATCH - step 3 would refuse this structure:")
            if missing:
                log("  %d query residues absent from the structure, first few: %s"
                    % (len(missing), missing[:8]))
            for i, q, s in bad[:8]:
                log("  residue %d: query has %s, structure has %s" % (i, q, s))
            log("  Fix: rerun with --renumber, or use --source esmfold, or trim the FASTA")


if __name__ == "__main__":
    main()

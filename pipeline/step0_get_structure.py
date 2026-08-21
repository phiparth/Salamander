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
import argparse, json, os, re, sys, urllib.request, urllib.error

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


# UniProt accession grammar, plus AlphaFold's own AF-<acc>-F1 form.
ACC_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$",
                    re.I)
AFID_RE = re.compile(r"^AF-[A-Z0-9]+-F[0-9]+$", re.I)


def check_accession(acc):
    """AlphaFold indexes by accession, not by gene symbol.

    Step 1 takes --gene ACHE, so it is natural to hand --uniprot the same thing here. That
    returns HTTP 400 'Invalid identifier format', which does not say what is wrong. Catching
    it before the call, and naming the difference, is worth the six lines.
    """
    acc = acc.strip()
    if ACC_RE.match(acc) or AFID_RE.match(acc):
        return acc

    # the README writes YOUR_ACCESSION where your value goes; it gets pasted as-is
    if acc.upper() in ("YOUR_ACCESSION", "YOUR-ACCESSION", "ACCESSION",
                       "YOUR_UNIPROT", "UNIPROT_ID", "YOUR_GENE"):
        sys.exit(
            "--uniprot %s is the placeholder from the README, not a real value.\n"
            "\n"
            "  Replace it with your protein's UniProt accession, e.g.:\n"
            "    --uniprot P22303      acetylcholinesterase (ACHE)\n"
            "    --uniprot P03372      estrogen receptor alpha (ESR1)\n"
            "\n"
            "  Look yours up by gene symbol at https://www.uniprot.org/ - the accession is\n"
            "  the short code like P22303 shown as 'Entry' in the results.\n"
            "\n"
            "  Benchmark proteins:\n"
            "    ESR1 P03372   TPH1 P17752   ACHE P22303   SIRT6 Q8N6T7\n"
            "    DNMT3A Q9Y6K1   PRSS1 P07477   VPS26A O75436   STXBP1 P61764" % acc)
    looks_like_gene = acc.isalnum() and acc.upper() == acc and not any(c.isdigit() for c in acc[:1])
    msg = ["%r is not a UniProt accession." % acc,
           "",
           "  AlphaFold DB is indexed by ACCESSION (P22303, Q8N6T7, O75436) - a letter,",
           "  then digits and letters, 6 or 10 characters. It is not the gene symbol.",
           ""]
    if looks_like_gene:
        msg += ["  That looks like a gene symbol. Note the two flags differ:",
                "    step 1  --gene ACHE       the gene symbol, for OrthoDB",
                "    step 0  --uniprot P22303  the accession, for AlphaFold",
                ""]
    msg += ["  Find the accession by searching the gene symbol at uniprot.org, or:",
            "    https://rest.uniprot.org/uniprotkb/search?query=gene:%s+AND+organism_id:9606&fields=accession"
            % acc,
            "",
            "  Accessions for the benchmark proteins:",
            "    ESR1 P03372   TPH1 P17752   ACHE P22303   SIRT6 Q8N6T7",
            "    DNMT3A Q9Y6K1   PRSS1 P07477   VPS26A O75436   STXBP1 P61764"]
    sys.exit("\n".join(msg))


def from_alphafold(acc, log):
    acc = check_accession(acc)
    log("querying AlphaFold DB for %s ..." % acc)
    try:
        raw = get(AFDB % acc, timeout=90)
    except urllib.error.HTTPError as e:
        body = ""
        try:
            body = e.read().decode("utf-8", "replace")[:200]
        except Exception:  # noqa: BLE001
            pass
        if e.code == 404:
            raise SystemExit("AlphaFold DB has no model for %s.\n"
                             "  Not every UniProt entry has one. Use --source esmfold to\n"
                             "  predict it instead, or fetch an experimental structure from\n"
                             "  rcsb.org and pass it to step 3 directly." % acc)
        raise SystemExit("AlphaFold DB refused the request for %s: HTTP %d\n  %s"
                         % (acc, e.code, body))
    meta = json.loads(raw.decode())
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


def extract_query(text, res, chain, seq, log, max_mismatch=0):
    """Locate the query inside a longer structure and cut it out, renumbered 1..L.

    AlphaFold models cover the whole UniProt entry, so a domain query never lines up.
    Finding the query inside the modelled sequence turns a full-length model into a usable
    one.

    An exact substring is tried first. With --max-mismatch N, the best-scoring window is
    accepted if it differs at no more than N positions. That covers the common case of an
    engineered construct: catalytically dead point mutants (a PTP catalytic Cys->Ser, say)
    differ from the UniProt entry at a handful of residues but are otherwise the same
    protein, and the wild-type backbone is a perfectly good template for them.
    """
    nums = sorted(n for c, n in res if c == chain)
    struct_seq = "".join(res[(chain, n)] for n in nums)
    i = struct_seq.find(seq)
    mism = []
    if i < 0:
        if max_mismatch <= 0 or len(struct_seq) < len(seq):
            return None
        best = (-1, -1)
        for off in range(len(struct_seq) - len(seq) + 1):
            m = sum(1 for k in range(len(seq)) if struct_seq[off + k] == seq[k])
            if m > best[0]:
                best = (m, off)
        m, i = best
        n_bad = len(seq) - m
        if n_bad > max_mismatch:
            log("  best window differs at %d positions, more than --max-mismatch %d"
                % (n_bad, max_mismatch))
            return None
        mism = [(k + 1, seq[k], struct_seq[i + k])
                for k in range(len(seq)) if struct_seq[i + k] != seq[k]]
    span = nums[i:i + len(seq)]
    log("  found the query at structure residues %d-%d; extracting and renumbering to 1-%d"
        % (span[0], span[-1], len(seq)))
    for k, want, got in mism:
        log("    NOTE position %d: your sequence has %s, the structure has %s"
            % (k, want, got))
    if mism:
        log("    -> %d position(s) differ. FoldX will run, but any mutation AT these" % len(mism))
        log("       positions is scored against the structure's residue, not yours.")
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
    p.add_argument("--max-mismatch", type=int, default=0,
                   help="accept a near-match window differing at up to N positions")
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
        cut = extract_query(text, res, a.chain, seq, log, a.max_mismatch)
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
        elif missing or len(bad) > a.max_mismatch:
            # step 3 verifies the residues it is about to mutate; wholesale disagreement
            # here means the frames do not correspond at all
            log("\nNUMBERING MISMATCH - step 3 would refuse this structure:")
            if missing:
                log("  %d query residues absent from the structure, first few: %s"
                    % (len(missing), missing[:8]))
            for i, q, s in bad[:8]:
                log("  residue %d: query has %s, structure has %s" % (i, q, s))
            log("  Fix: rerun with --renumber, or use --source esmfold, or trim the FASTA")
        else:
            # a handful of accepted differences: usable, with one caveat
            log("\nUSABLE, with %d accepted difference(s):" % len(bad))
            for i, q, s in bad[:8]:
                log("  residue %d: query has %s, structure has %s" % (i, q, s))
            log("  Step 3 checks only the residues it mutates, so this structure is fine")
            log("  UNLESS a proposed mutation lands on one of the positions above - it")
            log("  would then be scored against the structure's residue, not yours. Step 3")
            log("  verifies that and stops if so. Common cause: your sequence is the mature")
            log("  chain and the model is the full precursor, or vice versa.")


if __name__ == "__main__":
    main()

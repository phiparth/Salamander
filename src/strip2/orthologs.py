# -*- coding: utf-8 -*-
"""Ortholog retrieval - OrthoDB (default) and NCBI blastp (fallback).

Both return [(name, sequence)] with at most one sequence per genus. One-per-genus matters:
a family dominated by twelve near-identical mammalian sequences produces a conservation
profile that reflects sampling bias rather than real constraint.
"""
import re
import time
import urllib.parse
import urllib.request

from .core import clean_seq

ORTHODB = "https://data.orthodb.org/v12"
UA = {"User-Agent": "Mozilla/5.0 (strip2-thermostable-designer)"}


def _retry(fn, what, log, tries=3):
    for i in range(tries):
        try:
            return fn()
        except Exception as e:
            if i == tries - 1:
                raise
            log("  %s failed (%s); retry %d/%d" % (what, e, i + 1, tries - 1))
            time.sleep(3 * (i + 1))


def _identity(a, b):
    """Crude ungapped identity on the overlap - only used for coarse filtering."""
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    return sum(1 for i in range(n) if a[i] == b[i]) / float(n)


def _one_per_genus(members, query, n, min_identity, length_band, log):
    L = len(query)
    seen, out = set(), []
    for org, gid, seq in members:
        seq = clean_seq(seq)
        genus = org.split()[0].lower() if org else ""
        if not seq or not genus or genus in seen:
            continue
        if not (L / length_band <= len(seq) <= L * length_band):
            continue
        seen.add(genus)
        out.append(("%s (%s)" % (org, gid), seq))
    if min_identity > 0:
        keep = [h for h in out if _identity(query, h[1]) >= min_identity]
        # only enforce identity if it leaves a usable family
        if len(keep) >= 3:
            out = keep
        elif keep != out:
            log("  identity filter would leave %d orthologs; keeping the unfiltered set"
                % len(keep))
    return out[:n]


def orthodb_orthologs(query, n=10, min_identity=0.30, length_band=2.0, log=print,
                      og=None, gene=None):
    """OrthoDB v12 REST: /search -> orthogroups, /fasta -> members.

    The family is named, not inferred:

        gene  a gene name such as ESR1 - resolved through /search
        og    an orthogroup id such as 4385266at2759 - goes straight to /fasta

    There used to be a third route: hand /blast the query sequence and let OrthoDB work out
    which gene it is. That endpoint is the only one that runs an alignment on their server,
    and it returns HTTP 500 - its own crash, "rapsearch: unhandled exception during
    asyncio.run() shutdown" - on every API version while /search, /genesearch and /fasta
    answer normally. It has been removed rather than left as a default that always fails.
    """
    import json as _json
    query = clean_seq(query)

    def api(path):
        return _retry(lambda: urllib.request.urlopen(
            urllib.request.Request(ORTHODB + path, headers=UA), timeout=60).read().decode(),
            "OrthoDB API", log)

    if og:
        ogs = [og]
        log("[ortho] using orthogroup %s directly (skipping the sequence search)" % og)
    elif gene:
        log("[ortho] OrthoDB gene-name search for %r ..." % gene)
        try:
            d = _json.loads(api("/search?query=" + urllib.parse.quote(gene)))
        except Exception as e:  # noqa: BLE001 - reported, not raised
            log("[ortho] OrthoDB /search failed: %s" % e)
            return []
        ogs = [x for x in (d.get("data") or []) if re.match(r"^\d+at\d+$", str(x))]
        if not ogs:
            log("[ortho] no orthogroups matched the gene name %r" % gene)
            return []
        log("[ortho] %d orthogroup(s) matched: %s" % (len(ogs), ", ".join(ogs[:4])))
    else:
        log("[ortho] no family given. OrthoDB needs --gene: it cannot identify the family")
        log("[ortho] from the sequence alone (that endpoint is down on their side).")
        log("[ortho]   --gene ESR1     the gene symbol for your protein")
        return []

    best = None
    for og in ogs[:4]:                                 # try a few levels, keep the most diverse
        time.sleep(1.1)
        try:
            fa = api("/fasta?id=" + og)
        except Exception:
            continue
        members, meta, cur = [], None, ""
        for ln in fa.splitlines():
            if ln.startswith(">"):
                if meta is not None and cur:
                    members.append((meta.get("organism_name", ""),
                                    meta.get("pub_gene_id", "?"), cur))
                cur = ""
                try:
                    meta = _json.loads(ln[ln.find("{"):])
                except Exception:
                    meta = {}
            else:
                cur += ln.strip()
        if meta is not None and cur:
            members.append((meta.get("organism_name", ""), meta.get("pub_gene_id", "?"), cur))
        genera = {m[0].split()[0] for m in members if m[0]}
        log("[ortho]   %s (level %s): %d members, %d genera"
            % (og, og.split("at")[1], len(members), len(genera)))
        if best is None or len(genera) > best[1]:
            best = (og, len(genera), members)
        if len(genera) >= n * 2:
            break

    if best is None:
        return []
    og, _, members = best
    out = _one_per_genus(members, query, n, min_identity, length_band, log)
    log("[ortho] OrthoDB: %d orthologs from %s (one per genus)" % (len(out), og))
    return out


def blast_orthologs(query, n=10, min_identity=0.30, length_band=2.0, email=None,
                    hitlist=150, log=print):
    """NCBI blastp against nr, one hit per genus. Slow (1-7 min) and needs internet."""
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import Entrez
    if email:
        Entrez.email = email
    query = clean_seq(query)
    log("[ortho] NCBI blastp (%d aa), this takes a few minutes ..." % len(query))
    handle = _retry(lambda: NCBIWWW.qblast("blastp", "nr", query, hitlist_size=hitlist),
                    "NCBI blastp", log)
    rec = NCBIXML.read(handle)
    members = []
    for aln in rec.alignments:
        m = re.search(r"\[([^\]]+)\]\s*$", aln.hit_def)
        org = m.group(1) if m else ""
        seq = "".join(hsp.sbjct.replace("-", "") for hsp in aln.hsps[:1])
        members.append((org, aln.accession, seq))
    out = _one_per_genus(members, query, n, min_identity, length_band, log)
    log("[ortho] blastp: %d orthologs (one per genus)" % len(out))
    return out

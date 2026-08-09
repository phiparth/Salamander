# -*- coding: utf-8 -*-
"""Core primitives shared by every pipeline step.

Only what the strip2 pipeline needs: ProtT5 embedding, checkpoint loading, the MLP head,
and the positional encoding that must match training exactly.
"""
import os
import numpy as np

AA20 = "ACDEFGHIKLMNPQRSTVWY"
AA2IDX = {a: i for i, a in enumerate(AA20)}
IDX2AA = {i: a for a, i in AA2IDX.items()}
DEFAULT_PLM = os.environ.get("PLM", "Rostlab/prot_t5_xl_half_uniref50-enc")


def clean_seq(seq):
    """Uppercase, strip whitespace, map the non-standard residues ProtT5 does not know."""
    s = "".join(seq.split()).upper()
    return "".join(("X" if c in "UZOB" else c) for c in s if c.isalpha())


def read_fasta(path):
    """-> [(name, sequence)] ; tolerates blank lines and wrapped sequences."""
    out, name, buf = [], None, []
    for ln in open(path, encoding="utf-8"):
        ln = ln.rstrip("\n\r")
        if ln.startswith(">"):
            if name is not None:
                out.append((name, clean_seq("".join(buf))))
            name, buf = ln[1:].strip(), []
        elif ln.strip():
            buf.append(ln.strip())
    if name is not None:
        out.append((name, clean_seq("".join(buf))))
    return out


def write_fasta(path, records, width=60):
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(">%s\n" % name)
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


# --------------------------------------------------------------------------------------
#  positional encoding - MUST match training (train_prok.py assemble())
#  sinusoid of the residue index divided by the sequence length.
# --------------------------------------------------------------------------------------
def pos_enc(pos, length, dim=16):
    p = pos / max(length, 1)
    freqs = np.exp(np.linspace(0, np.log(50.0), dim // 2))
    return np.concatenate([np.sin(p * freqs), np.cos(p * freqs)]).astype(np.float32)


# --------------------------------------------------------------------------------------
#  ProtT5
# --------------------------------------------------------------------------------------
_PLM_CACHE = {}


def load_plm(plm=None, log=print):
    """Load T5 encoder + tokenizer once per process. `plm` may be a local directory."""
    import torch
    from transformers import T5EncoderModel, T5Tokenizer
    plm = plm or DEFAULT_PLM
    if plm in _PLM_CACHE:
        return _PLM_CACHE[plm]
    log("loading ProtT5 from %s ..." % plm)
    tok = T5Tokenizer.from_pretrained(plm, do_lower_case=False, legacy=True)
    mdl = T5EncoderModel.from_pretrained(plm)
    mdl = mdl.eval().to("cuda" if torch.cuda.is_available() else "cpu")
    _PLM_CACHE[plm] = (tok, mdl)
    return tok, mdl


def embed(seq, plm=None, log=print):
    """-> (per_residue [L,1024] float32, pooled [1024] float32)

    Pooling is norm-weighted, exactly as in training: each residue vector is weighted by
    its own L2 norm before averaging. A plain mean would not reproduce the training input.
    """
    import torch
    tok, mdl = load_plm(plm, log)
    s = clean_seq(seq)
    spaced = " ".join(list(s))
    enc = tok(spaced, return_tensors="pt", add_special_tokens=True)
    enc = {k: v.to(mdl.device) for k, v in enc.items()}
    with torch.no_grad():
        out = mdl(**enc).last_hidden_state[0].cpu().numpy().astype(np.float32)
    emb = out[:len(s)]
    w = np.linalg.norm(emb, axis=1, keepdims=True) + 1e-8
    pooled = (emb * w).sum(0) / w.sum()
    return emb, pooled.astype(np.float32)


# --------------------------------------------------------------------------------------
#  the head
# --------------------------------------------------------------------------------------
def _make_head(in_dim, hidden, dropout, out=20):
    """`hidden` may be an int (two equal layers) or a list of widths."""
    import torch.nn as nn
    dims = [hidden, hidden] if isinstance(hidden, int) else list(hidden)
    layers, prev = [], in_dim
    for h in dims:
        layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(dropout)]
        prev = h
    layers += [nn.Linear(prev, out)]
    return nn.Sequential(*layers)


def load_head(path):
    """-> (net, meta). Checkpoints store their own geometry, so nothing is hardcoded."""
    import torch
    ck = torch.load(path, map_location="cpu", weights_only=False)
    hidden = ck.get("hidden_dims") or ck["hidden"]
    net = _make_head(ck["in_dim"], hidden, ck["dropout"])
    sd = ck["state_dict"]
    # training wrapped the Sequential in self.net, so the keys carry a "net." prefix
    if all(k.startswith("net.") for k in sd):
        sd = {k[4:]: v for k, v in sd.items()}
    net.load_state_dict(sd)
    net.eval()
    meta = {k: v for k, v in ck.items() if k != "state_dict"}
    return net, meta


def residue_probs(net, meta, emb, pooled, positions=None):
    """20-way softmax at each requested residue index (0-based). -> {pos: np.array[20]}"""
    import torch
    L = emb.shape[0]
    pos_dim = meta.get("pos_dim", 16)
    mode = meta.get("input_mode", "pooled_plus_residue")
    cond = meta.get("condition_dtm", False)
    out = {}
    idx = range(L) if positions is None else positions
    with torch.no_grad():
        for r in idx:
            pe = pos_enc(r, L, pos_dim)
            if mode == "pooled":
                x = np.concatenate([pooled, pe])
            elif mode == "per_residue":
                x = np.concatenate([emb[r], pe])
            else:
                x = np.concatenate([pooled, emb[r], pe])
            if cond:
                x = np.concatenate([x, [np.float32(30.0 / 100.0)]])
            logits = net(torch.from_numpy(x[None].astype(np.float32)))
            out[r] = torch.softmax(logits, 1)[0].numpy()
    return out

#!/usr/bin/env python3
"""Stage 3 of training - the streaming pairwise trainer that produced the prok head.

One training example = one (low-Tm protein, position) pair, with the hot homolog's residue
as the label:

    input   [ pooled(lo) 1024 | emb(lo)[r] 1024 | pos_enc(r/L) 16 ]  =  2064
    model   MLP 2064 -> h -> h -> 20
    label   the residue the HOT partner has at that column
    loss    class-weighted cross-entropy, each example scaled by sqrt(dTm)

Streaming rather than eager: the naive version materialises the whole feature matrix and
needs ~100 GB at full scale. This holds one [N,1024] per-residue array plus the pooled
vectors and assembles each batch on the fly.

    python training/train_prok.py --work_dir runs/prok --epochs 30

The train/validation split is by CLUSTER, never by example. Alignment errors are correlated
within a cluster, so a random split leaks homologs across the boundary and inflates accuracy.
"""
import argparse, gc, json, os
import numpy as np


def main():
    p = argparse.ArgumentParser(description="train the pairwise head")
    p.add_argument("--work_dir", required=True)
    p.add_argument("--out", default=None, help="checkpoint path (default work_dir/head.pt)")
    p.add_argument("--plm", default="prot_t5_xl", help="recorded in the checkpoint only")
    p.add_argument("--pos_dim", type=int, default=16)
    p.add_argument("--hidden", type=int, default=512)
    p.add_argument("--dropout", type=float, default=0.3)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--batch", type=int, default=256)
    p.add_argument("--epochs", type=int, default=30)
    p.add_argument("--val_frac", type=float, default=0.15)
    p.add_argument("--dtm_weighting", default="sqrt", choices=["none", "linear", "sqrt"])
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()

    import torch
    import torch.nn as nn

    wd = a.work_dir
    cache = os.path.join(wd, "emb")
    records = json.load(open(os.path.join(wd, "pairs.json")))

    hashes, h2i = [], {}
    for rec in records:
        h = rec["lo_hash"]
        if h not in h2i:
            h2i[h] = len(hashes); hashes.append(h)
    print("[pair] %d records, %d unique low-Tm sequences" % (len(records), len(hashes)),
          flush=True)

    emb_store = [None] * len(hashes)
    pooled = np.zeros((len(hashes), 1024), np.float32)
    lolen = np.zeros(len(hashes), np.int32)
    miss = 0
    for h, i in h2i.items():
        f = os.path.join(cache, h + ".npy")
        if not os.path.exists(f):
            miss += 1; continue
        e = np.load(f)
        emb_store[i] = e
        lolen[i] = e.shape[0]
        w = np.linalg.norm(e, axis=1, keepdims=True) + 1e-8   # norm-weighted pooling
        pooled[i] = (e * w).sum(0) / w.sum()
    print("[pair] embeddings loaded (%d missing from cache)" % miss, flush=True)
    if miss:
        print("[pair] run training/embed_sequences.py first if that number is large", flush=True)

    N = sum(1 for rec in records if emb_store[h2i[rec["lo_hash"]]] is not None
            for r, _, _ in rec["targets"] if r < lolen[h2i[rec["lo_hash"]]])
    print("[pair] %d targets; per-residue matrix %.1f GB" % (N, N * 1024 * 4 / 1e9), flush=True)

    feat = np.zeros((N, 1024), np.float32)
    lo_idx = np.zeros(N, np.int32); rr = np.zeros(N, np.int32)
    yy = np.zeros(N, np.int8); dtm = np.zeros(N, np.float32); clu = np.zeros(N, np.int32)
    cl_names, cl2i = [], {}
    k = 0
    for rec in records:
        i = h2i[rec["lo_hash"]]; e = emb_store[i]
        if e is None:
            continue
        L = e.shape[0]; d = float(rec["dtm"]); c = str(rec["cluster"])
        if c not in cl2i:
            cl2i[c] = len(cl_names); cl_names.append(c)
        ci = cl2i[c]
        for r, ta, _ in rec["targets"]:
            if r >= L:
                continue
            feat[k] = e[r]; lo_idx[k] = i; rr[k] = r; yy[k] = ta; dtm[k] = d; clu[k] = ci
            k += 1
    emb_store = None
    gc.collect()
    lo_len_per_target = lolen[lo_idx]

    if a.dtm_weighting == "none":
        sw = np.ones(N, np.float32)
    elif a.dtm_weighting == "linear":
        sw = np.clip(dtm, 0, None) + 1e-6
    else:
        sw = np.sqrt(np.clip(dtm, 0, None)) + 1e-6
    sw = (sw / sw.mean()).astype(np.float32)

    rng = np.random.default_rng(a.seed)
    uniq = np.arange(len(cl_names)); rng.shuffle(uniq)
    n_val = max(1, int(len(uniq) * a.val_frac))
    is_val = np.isin(clu, uniq[:n_val])
    tr_idx = np.where(~is_val)[0]; va_idx = np.where(is_val)[0]
    print("[pair] train %d / val %d  (%d train clusters / %d val clusters)"
          % (len(tr_idx), len(va_idx), len(uniq) - n_val, n_val), flush=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    in_dim = 1024 + 1024 + a.pos_dim
    freqs = np.exp(np.linspace(0, np.log(50.0), a.pos_dim // 2))

    def assemble(idx):
        pl = pooled[lo_idx[idx]]
        q = rr[idx] / np.maximum(lo_len_per_target[idx], 1)
        ang = q[:, None] * freqs[None, :]
        pe = np.concatenate([np.sin(ang), np.cos(ang)], 1).astype(np.float32)
        return np.concatenate([pl, feat[idx], pe], 1)

    class Head(nn.Module):
        def __init__(self, d, h, dr):
            super().__init__()
            self.net = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(dr),
                                     nn.Linear(h, h), nn.ReLU(), nn.Dropout(dr),
                                     nn.Linear(h, 20))

        def forward(self, x):
            return self.net(x)

    model = Head(in_dim, a.hidden, a.dropout).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=1e-4)
    counts = np.bincount(yy[tr_idx].astype(np.int64), minlength=20).astype(np.float32)
    cw = torch.tensor(counts.sum() / (counts + 1.0), device=device)
    cw = cw / cw.mean()
    lossf = nn.CrossEntropyLoss(weight=cw, reduction="none")

    base = (yy[va_idx] == np.bincount(yy[tr_idx].astype(np.int64)).argmax()).mean()
    print("[pair] majority-class baseline val_acc = %.3f" % base, flush=True)

    def evaluate():
        model.eval(); ok = tot = 0
        with torch.no_grad():
            for st in range(0, len(va_idx), a.batch):
                b = va_idx[st:st + a.batch]
                xb = torch.from_numpy(assemble(b)).to(device)
                ok += int((model(xb).argmax(1).cpu().numpy() == yy[b]).sum()); tot += len(b)
        return ok / max(tot, 1)

    out_path = a.out or os.path.join(wd, "head.pt")
    best = 0.0
    for ep in range(1, a.epochs + 1):
        model.train()
        order = rng.permutation(len(tr_idx)); run = 0.0
        for st in range(0, len(tr_idx), a.batch):
            b = tr_idx[order[st:st + a.batch]]
            xb = torch.from_numpy(assemble(b)).to(device)
            yb = torch.from_numpy(yy[b].astype(np.int64)).to(device)
            wb = torch.from_numpy(sw[b]).to(device)
            opt.zero_grad()
            loss = (lossf(model(xb), yb) * wb).mean()
            loss.backward(); opt.step()
            run += loss.item() * len(b)
        acc = evaluate()
        print("[pair] epoch %02d  loss=%.3f  val_acc=%.3f" % (ep, run / len(tr_idx), acc),
              flush=True)
        if acc > best:
            best = acc
            torch.save({"state_dict": model.state_dict(), "in_dim": in_dim,
                        "hidden": a.hidden, "dropout": a.dropout,
                        "input_mode": "pooled_plus_residue", "pos_dim": a.pos_dim,
                        "plm": a.plm, "condition_dtm": False,
                        "dtm_weighting": a.dtm_weighting, "val_acc": acc}, out_path)
    print("[pair] best val_acc=%.3f (baseline %.3f); saved %s" % (best, base, out_path))


if __name__ == "__main__":
    main()

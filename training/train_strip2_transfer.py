#!/usr/bin/env python3
"""Transfer stage 2 - fine-tune the prok head on yeast pairs to produce strip2.

The prok head is  Linear(2064,512) ReLU Drop | Linear(512,512) ReLU Drop | Linear(512,20).
strip2 is the "strip2_block+clf" configuration:

    block 1  Linear(2064,512)   FROZEN    - keeps the general ProtT5 -> residue mapping
    block 2  Linear(512,512)    retrained
    classifier Linear(512,20)   retrained

Freezing block 1 is what makes this transfer rather than retraining. The yeast pair set is
much smaller than the prokaryote one, and letting the whole network move overfits it; keeping
the first projection fixed preserves what the large set taught and adapts only the head.

    python training/train_strip2_transfer.py --work_dir runs/yeast \
           --pretrained runs/prok/head.pt --epochs 60 --out models/head_strip2.pt

--frozen-blocks lets you try the alternatives: 0 = full fine-tune, 2 = classifier only
(that is the "strip1" linear-probe variant, which scored lower).
"""
import argparse, gc, json, os
import numpy as np


def main():
    p = argparse.ArgumentParser(description="fine-tune prok -> strip2 on yeast pairs")
    p.add_argument("--work_dir", required=True, help="yeast pairs.json + emb/")
    p.add_argument("--pretrained", required=True, help="prok head.pt")
    p.add_argument("--out", default=None, help="default work_dir/head_strip2.pt")
    p.add_argument("--frozen-blocks", type=int, default=1,
                   help="1 = strip2 (default), 0 = full fine-tune, 2 = classifier only")
    p.add_argument("--plm", default="prot_t5_xl")
    p.add_argument("--pos_dim", type=int, default=16)
    p.add_argument("--hidden", type=int, default=512)
    p.add_argument("--dropout", type=float, default=0.3)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--batch", type=int, default=256)
    p.add_argument("--epochs", type=int, default=60)
    p.add_argument("--val_frac", type=float, default=0.15)
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
    print("[strip2] %d yeast pair records, %d unique sequences" % (len(records), len(hashes)),
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
        emb_store[i] = e; lolen[i] = e.shape[0]
        w = np.linalg.norm(e, axis=1, keepdims=True) + 1e-8
        pooled[i] = (e * w).sum(0) / w.sum()
    print("[strip2] embeddings loaded (%d missing)" % miss, flush=True)

    N = sum(1 for rec in records if emb_store[h2i[rec["lo_hash"]]] is not None
            for r, _, _ in rec["targets"] if r < lolen[h2i[rec["lo_hash"]]])
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
        for r, ta, _ in rec["targets"]:
            if r >= L:
                continue
            feat[k] = e[r]; lo_idx[k] = i; rr[k] = r; yy[k] = ta
            dtm[k] = d; clu[k] = cl2i[c]
            k += 1
    emb_store = None
    gc.collect()
    lo_len_per_target = lolen[lo_idx]
    print("[strip2] %d targets over %d clusters" % (N, len(cl_names)), flush=True)

    sw = np.sqrt(np.clip(dtm, 0, None)) + 1e-6
    sw = (sw / sw.mean()).astype(np.float32)
    rng = np.random.default_rng(a.seed)
    uniq = np.arange(len(cl_names)); rng.shuffle(uniq)
    n_val = max(1, int(len(uniq) * a.val_frac))
    is_val = np.isin(clu, uniq[:n_val])
    tr_idx = np.where(~is_val)[0]; va_idx = np.where(is_val)[0]
    print("[strip2] train %d / val %d (%d/%d clusters)"
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

    # ---- rebuild the prok architecture and load its weights ----
    def fresh():
        return nn.Sequential(nn.Linear(in_dim, a.hidden), nn.ReLU(), nn.Dropout(a.dropout),
                             nn.Linear(a.hidden, a.hidden), nn.ReLU(), nn.Dropout(a.dropout),
                             nn.Linear(a.hidden, 20))

    ck = torch.load(a.pretrained, map_location="cpu", weights_only=False)
    sd = ck["state_dict"]
    if all(key.startswith("net.") for key in sd):
        sd = {key[4:]: v for key, v in sd.items()}
    model = fresh()
    model.load_state_dict(sd)
    print("[strip2] loaded prok weights from %s (val_acc %s)"
          % (os.path.basename(a.pretrained), ck.get("val_acc")), flush=True)

    # a "block" is Linear+ReLU+Dropout, so 3 modules each
    for i, mod in enumerate(model):
        if i < a.frozen_blocks * 3:
            for prm in mod.parameters():
                prm.requires_grad = False
    frozen = sum(1 for prm in model.parameters() if not prm.requires_grad)
    train_n = sum(1 for prm in model.parameters() if prm.requires_grad)
    print("[strip2] frozen_blocks=%d -> %d frozen tensors, %d trainable"
          % (a.frozen_blocks, frozen, train_n), flush=True)
    model.to(device)

    params = [prm for prm in model.parameters() if prm.requires_grad]
    opt = torch.optim.AdamW(params, lr=a.lr, weight_decay=1e-4)
    counts = np.bincount(yy[tr_idx].astype(np.int64), minlength=20).astype(np.float32)
    cw = torch.tensor(counts.sum() / (counts + 1.0), device=device)
    cw = cw / cw.mean()
    lossf = nn.CrossEntropyLoss(weight=cw, reduction="none")

    def evaluate(m):
        m.eval(); ok = tot = 0
        with torch.no_grad():
            for st in range(0, len(va_idx), a.batch):
                b = va_idx[st:st + a.batch]
                xb = torch.from_numpy(assemble(b)).to(device)
                ok += int((m(xb).argmax(1).cpu().numpy() == yy[b]).sum()); tot += len(b)
        return ok / max(tot, 1)

    base = (yy[va_idx] == np.bincount(yy[tr_idx].astype(np.int64)).argmax()).mean()
    print("[strip2] majority-class baseline      val_acc=%.4f" % base, flush=True)
    print("[strip2] prok model, no fine-tuning   val_acc=%.4f" % evaluate(model), flush=True)

    out_path = a.out or os.path.join(wd, "head_strip2.pt")
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
        acc = evaluate(model)
        print("[strip2] epoch %02d  loss=%.3f  val_acc=%.4f" % (ep, run / len(tr_idx), acc),
              flush=True)
        if acc > best:
            best = acc
            torch.save({"state_dict": model.state_dict(), "in_dim": in_dim,
                        "hidden": a.hidden, "dropout": a.dropout,
                        "input_mode": "pooled_plus_residue", "pos_dim": a.pos_dim,
                        "plm": a.plm, "condition_dtm": False, "dtm_weighting": "sqrt",
                        "config": "strip2_block+clf" if a.frozen_blocks == 1
                                  else "frozen_blocks=%d" % a.frozen_blocks,
                        "transfer": "yeast", "val_acc": acc}, out_path)
    print("[strip2] best val_acc=%.4f (baseline %.4f); saved %s" % (best, base, out_path))


if __name__ == "__main__":
    main()

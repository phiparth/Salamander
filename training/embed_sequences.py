#!/usr/bin/env python3
"""Stage 2 of training - embed every unique low-Tm sequence with ProtT5.

One .npy per sequence, named by its hash, so stage 3 can stream them. Writes are atomic
(temp file + rename): a truncated .npy from an interrupted run silently corrupts training.

    python training/embed_sequences.py --work_dir runs/prok --plm /path/to/prot_t5_xl

Batching is length-aware because ProtT5 attention is O(L^2) - a fixed batch size OOMs on
the longest sequences. Anything that still fails on GPU is retried on CPU.
"""
import argparse, hashlib, os, sys
import numpy as np


def seq_hash(seq):
    return hashlib.sha1(seq.encode()).hexdigest()[:16]


def main():
    p = argparse.ArgumentParser(description="cache ProtT5 embeddings for training")
    p.add_argument("--work_dir", required=True)
    p.add_argument("--plm", default=os.environ.get("PLM",
                                                   "Rostlab/prot_t5_xl_half_uniref50-enc"))
    p.add_argument("--max_tokens", type=int, default=6000,
                   help="approximate residues per batch; lower it if you hit OOM")
    a = p.parse_args()

    import torch
    from transformers import T5EncoderModel, T5Tokenizer

    cache = os.path.join(a.work_dir, "emb")
    os.makedirs(cache, exist_ok=True)
    seqs = [s.strip() for s in
            open(os.path.join(a.work_dir, "unique_seqs.txt"), encoding="utf-8") if s.strip()]
    todo = [s for s in seqs if not os.path.exists(os.path.join(cache, seq_hash(s) + ".npy"))]
    print("%d unique sequences, %d still to embed" % (len(seqs), len(todo)), flush=True)
    if not todo:
        return

    tok = T5Tokenizer.from_pretrained(a.plm, do_lower_case=False, legacy=True)
    mdl = T5EncoderModel.from_pretrained(a.plm).eval()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    mdl = mdl.to(dev)
    print("device %s" % dev, flush=True)

    todo.sort(key=len)                                   # group similar lengths together
    batches, cur = [], []
    for s in todo:
        if cur and (len(cur) + 1) * len(s) > a.max_tokens:
            batches.append(cur); cur = []
        cur.append(s)
    if cur:
        batches.append(cur)
    print("%d batches" % len(batches), flush=True)

    def embed_batch(bs, device):
        spaced = [" ".join(list(s)) for s in bs]
        enc = tok(spaced, return_tensors="pt", padding=True, add_special_tokens=True)
        enc = {k: v.to(device) for k, v in enc.items()}
        with torch.no_grad():
            out = mdl.to(device)(**enc).last_hidden_state.cpu().numpy().astype(np.float32)
        return [out[i][:len(s)] for i, s in enumerate(bs)]

    done = 0
    for bi, bs in enumerate(batches, 1):
        try:
            embs = embed_batch(bs, dev)
        except RuntimeError as e:
            print("  batch %d failed on %s (%s); retrying one-by-one on CPU"
                  % (bi, dev, str(e)[:60]), flush=True)
            embs = []
            for s in bs:
                embs += embed_batch([s], "cpu")
            mdl.to(dev)
        for s, e in zip(bs, embs):
            h = seq_hash(s)
            tmp = os.path.join(cache, h + ".tmp.npy")     # np.save appends .npy, so the
            np.save(tmp, e)                               # temp name must already end in it
            os.replace(tmp, os.path.join(cache, h + ".npy"))
            done += 1
        if bi % 20 == 0 or bi == len(batches):
            print("  %d/%d batches, %d sequences" % (bi, len(batches), done), flush=True)

    print("\nembeddings in %s" % cache)


if __name__ == "__main__":
    main()

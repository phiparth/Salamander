# strip2 — thermostable protein design

Redesign a protein for higher melting temperature, without touching its active site.

`strip2` is a per-residue classifier on frozen [ProtT5](https://github.com/agemagician/ProtTrans)
embeddings. It learns from evolution rather than from thermodynamics: given a family of
homologs where some members come from hot organisms and some from cold, it learns which
residue a *hotter* relative would carry at each position.

Because it knows sequence but not structure, the model alone proposes too many mutations and
some of them are structurally destructive. The pipeline therefore runs three filters around
it — conservation, model confidence, and FoldX stability — and only mutations that survive
all three reach the final sequence.

```
query sequence
      │
      ▼  step 1   orthologs → MSA → freeze conserved sites          (protects the active site)
      ▼  step 2   strip2 → 20-way softmax → set1* confidence cut    (keeps confident calls)
      ▼  step 3   FoldX single-point ΔΔG → drop destabilising       (keeps the fold intact)
      ▼  step 4   assemble, score the whole design, report
final sequence
```

---

## Requirements

You must supply two things yourself; neither is redistributable:

| | |
|---|---|
| **ProtT5** | `Rostlab/prot_t5_xl_half_uniref50-enc`, or a local copy. Set `PLM=/path/to/model` or pass `--plm`. ~11 GB for the full checkpoint; only the encoder is loaded. |
| **FoldX 5** | Academic licence from [foldxsuite.crg.eu](https://foldxsuite.crg.eu). Set `FOLDX=/path/to/foldx` or pass `--foldx`. Keep its `molecules/` directory beside the executable. |

Then:

```bash
pip install -r requirements.txt
```

A GPU makes embedding faster but is not required — everything runs on CPU.

---

## Quick start

```bash
python pipeline/run_pipeline.py \
    --query examples/era_lbd.fasta \
    --pdb   structures/era_lbd.pdb \
    --out   work/era \
    --cons-thresh 0.80 \
    --wt-ratio 1.20 \
    --ddg-cut 0.05 \
    --repair
```

Results land in `work/era/`: `final.fasta`, `final.json`, and one JSON per step so you can
inspect or restart anywhere.

Without `--pdb` the run stops after step 2 and gives you the set1\* design — useful when you
have no structure, but the stability filter is the step that catches the model's worst ideas,
so treat that output as provisional.

---

## The steps

### Step 1 — freeze the conserved sites

```bash
python pipeline/step1_conserved_sites.py --query Q.fasta --out work/Q \
       --source orthodb --n-orthologs 10 --cons-thresh 0.80
```

Finds orthologs (OrthoDB by default; `--source blast` for NCBI blastp, `--source local
--family F.fasta` to supply your own), aligns them with FAMSA, and marks a column conserved
when its modal non-gap residue covers at least `--cons-thresh` of the non-gap residues.
Those positions are frozen for the rest of the run.

At most one sequence per genus is kept. A family of twelve near-identical mammalian
sequences would otherwise produce a conservation profile that reflects sampling, not
constraint.

| flag | default | effect |
|---|---|---|
| `--cons-thresh` | 0.80 | higher → fewer frozen, more aggressive design |
| `--n-orthologs` | 10 | more orthologs → better-estimated conservation, slower |
| `--min-identity` | 0.30 | rejects distant hits that align badly |
| `--length-band` | 2.0 | keeps orthologs within ×2 of the query length |
| `--extra-freeze` | — | force-freeze residues you know matter, e.g. `"25 159 175"` |

On ERα-LBD at 0.80 this freezes 119 of 250 residues. **PROSS never mutated inside this band
either** in our benchmark, which is a useful sanity check that the threshold is sane.

### Step 2 — design and shortlist

```bash
python pipeline/step2_design_set1star.py --work work/Q --p-min 0.10 --wt-ratio 1.20
```

strip2 emits a 20-way softmax at each mutable position. A mutation is proposed when the
top-1 residue differs from wild type. **set1\*** then keeps only:

```
P(mutant) > 0.10                    the model must actually commit
P(mutant) ≥ 1.20 × P(wild type)     and must clearly prefer it to what is there
```

The ratio is the important half — an absolute cut alone keeps mutations the model rates
barely above wild type. Raise `--wt-ratio` for a smaller, more conservative design.

### Step 3 — FoldX stability filter

```bash
python pipeline/step3_foldx_filter.py --work work/Q --pdb Q.pdb --ddg-cut 0.05 --repair --jobs 6
```

Each surviving mutation is scored as an **isolated point mutation**. Anything above
`--ddg-cut` is reverted to wild type.

This step exists because strip2 has no structural sense. It has learned that proline
substitutions are a thermostabilising move — real thermophiles do carry them — but not
*where* the backbone tolerates one. In our benchmark the identical substitution type sat at
both extremes: **Gly→Pro at −2.85 kcal/mol** on ERα, **Gln→Pro at +4.60** on hDnmt3a. A
cheap single-point screen separates them.

The PDB is checked against the query sequence before any compute is spent; a numbering
mismatch aborts with the offending residues listed.

> **On `--ddg-cut 0.05`**: this is strict. FoldX's own error on a single-point ΔΔG is around
> 1 kcal/mol and its run-to-run noise is 0.2–0.5, so a 0.05 threshold sits well inside the
> uncertainty — which mutations fall either side is partly determined by noise. It yields a
> small, safe design. We validated **0.5** empirically across five proteins (it kept ~53% of
> set1\* mutations and every design came out net stabilising). Use 0.05 when you want to be
> conservative, 0.5 when you want more mutations, and treat neither as a physical constant.

### Step 4 — final sequence

```bash
python pipeline/step4_final_sequence.py --work work/Q --combined
```

Writes `final.fasta` and a report. With `--combined` it also scores the whole surviving set
as **one multi-mutant**, which is not the sum of the singles — mutations interact. On our
benchmark the combined value ran 3.2 kcal/mol *below* additive on hAChE and 1.7 *above* it on
TPH1. If you quote a ΔΔG for the design, quote this one.

---

## Training

The distributed head is `models/head_strip2.pt`. To reproduce the base model from scratch:

```bash
# 1. pairs from the clustered, Tm-annotated table  (excludes yeast → the "prok" model)
python training/build_pairs.py --csv training/data/combined_clustered_proteins.csv \
       --work_dir runs/prok --exclude_taxid 559292

# 2. cache ProtT5 embeddings for every unique low-Tm sequence
python training/embed_sequences.py --work_dir runs/prok --plm /path/to/prot_t5_xl

# 3. train
python training/train_prok.py --work_dir runs/prok --epochs 30
```

**How the supervision works.** No measured ΔΔG or ΔTm is ever used as a label. Within a
cluster of homologs, a cool protein and a hotter one are paired; at every non-conserved
alignment column the label is the residue the *hot* protein carries. The cool protein is
embedded; the hot one only supplies labels and never enters the network.

One example is:

```
input   [ pooled(lo) 1024 | emb(lo)[r] 1024 | pos_enc(r/L) 16 ]  =  2064
model   MLP 2064 → 512 → 512 → 20
label   the hot partner's residue at that column
loss    class-weighted cross-entropy × sqrt(ΔTm), normalised to mean 1
```

Three details that matter:

- **One MSA per cluster, built from all members at once** — not pairwise. Conserved-column
  detection needs a common frame, and a register agreed by seven sequences is far more
  reliable than one agreed by two.
- **Roughly a quarter of labels say "keep"** — the hot partner often has the same residue.
  The model is a 20-way classifier, not a mutation proposer, which is why it frequently
  returns wild type and emits nothing.
- **The split is by cluster, never by example.** Alignment errors are correlated within a
  cluster; a random split leaks homologs across the boundary and inflates accuracy.

Validation accuracy is ~0.24–0.26 against a majority-class baseline of ~0.10. That looks
modest because it is genuine held-out-family generalisation.

`training/data/combined_clustered_proteins.csv` (40 MB) is the clustered, Tm-annotated
source table. `--exclude_taxid 559292` drops *S. cerevisiae* and gives the prokaryote-only
base model; the released strip2 head is that model fine-tuned on yeast afterwards.

---

## Repository layout

```
models/head_strip2.pt         the released model
src/strip2/core.py            embedding, checkpoint loading, positional encoding
src/strip2/orthologs.py       OrthoDB and blastp ortholog retrieval
pipeline/step1..step4         the four stages
pipeline/run_pipeline.py      all four in sequence
training/build_pairs.py       stage 1 — pairs.json
training/embed_sequences.py   stage 2 — ProtT5 cache
training/train_prok.py        stage 3 — the trainer
training/data/                the training table
docs/method.md                design rationale and benchmark results
```

---

## Limitations

Worth reading before trusting an output.

- **Positions are chosen independently.** The head predicts each residue from its own
  context, so it cannot design cooperativity between mutations the way a combinatorial
  Rosetta protocol can. On ERα, PROSS's 24 mutations sum to +5.28 kcal/mol as isolated
  singles yet the assembled design scores −0.60 — about 5.9 kcal/mol gained purely from
  combination. strip2 has no mechanism for that.
- **ΔΔG is a proxy for the wrong observable.** These filters optimise folding stability;
  the model predicts thermostability. They correlate but are not the same, and no step here
  verifies that predicted ΔTm survives the filtering.
- **Multi-mutant ΔΔG degrades with mutation count.** Both FoldX and Rosetta relax from the
  wild-type backbone, so a design changing >10% of the sequence is increasingly measuring
  "how badly do these residues fit the original backbone". Treat large designs as
  directional.
- **N-terminal bias.** The positional encoding cannot distinguish residue 1 of a 250-mer from
  residue 1 of a 600-mer, and initiator methionine dominates that position in training. Check
  the first few residues of any design.
- **Nothing here is experimentally validated.** Every number in this repository is
  computational.

---

## Licence

MIT — see `LICENSE`. ProtT5 and FoldX carry their own licences and are not included.

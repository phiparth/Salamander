# Salamander

**Evolution-guided thermostable protein design.** Give it a protein sequence; it gives you back
a redesigned sequence predicted to melt at a higher temperature, with the active site left
alone.

*Named for the creature that was said to live in fire.*

---

## Contents

1. [What you need before you start](#1-what-you-need-before-you-start)
2. [Install](#2-install)
3. [Get ProtT5](#3-get-prott5-required)
4. [Get FoldX](#4-get-foldx-required-for-steps-3-and-4)
5. [Get a structure for your protein](#5-get-a-structure-for-your-protein)
6. [Run the example](#6-run-the-example-end-to-end)
7. [Run it on your own protein](#7-run-it-on-your-own-protein)
8. [What every file means](#8-what-every-file-means)
9. [Every option, explained](#9-every-option-explained)
10. [Training from scratch](#10-training-from-scratch-optional)
11. [Troubleshooting](#11-troubleshooting)
12. [How it works](#12-how-it-works)

---

## 1. What you need before you start

| | What | Where to get it | Size | Needed for |
|---|---|---|---|---|
| ✅ | **Python 3.9 – 3.12** | [python.org/downloads](https://www.python.org/downloads/) | — | everything |
| ✅ | **Git** | [git-scm.com/downloads](https://git-scm.com/downloads) | — | cloning this repo |
| ✅ | **ProtT5 model** | [huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc](https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc) | ~2.5 GB | steps 2 and all training |
| ✅ | **FoldX 5** | [foldxsuite.crg.eu](https://foldxsuite.crg.eu/) — free academic licence, registration required | ~50 MB | steps 3 and 4 |
| ✅ | **A PDB structure** of your protein | [rcsb.org](https://www.rcsb.org/) or [alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk/) | ~1 MB | steps 3 and 4 |
| ⬜ | GPU with CUDA | — | — | optional, ~10× faster embedding |
| ⬜ | Internet connection | — | — | step 1 only (ortholog search) |

**Disk space**: about 3 GB total, almost all of it ProtT5.

**Time for one protein**: 2–5 minutes for steps 1–2, then 10–60 minutes for steps 3–4
depending on protein size and how many mutations survive. `RepairPDB` in step 3 is the slow
part.

---

## 2. Install

Open a terminal (Command Prompt or PowerShell on Windows, Terminal on macOS/Linux) and run:

```bash
git clone https://github.com/phiparth/Salamander.git
cd Salamander
```

Create an isolated Python environment so this does not disturb your other projects:

**Windows**
```bash
python -m venv .venv
.venv\Scripts\activate
```

**macOS / Linux**
```bash
python3 -m venv .venv
source .venv/bin/activate
```

Your prompt should now start with `(.venv)`. Install the dependencies:

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

That installs `torch`, `transformers`, `sentencepiece`, `protobuf`, `numpy`, `pandas`,
`pyfamsa` and `biopython`. It takes a few minutes, mostly downloading PyTorch.

Check it worked:

```bash
python pipeline/step1_conserved_sites.py --help
```

You should see the usage text. If you see `ModuleNotFoundError`, the virtual environment is
not active — rerun the `activate` line above.

> **Every command in this README assumes you are in the `Salamander` folder with `(.venv)`
> active.** If you close the terminal, `cd` back and re-run `activate`.

---

## 3. Get ProtT5 (required)

ProtT5 turns a protein sequence into the numeric embedding the model reads. It is **not**
included here — it is 2.5 GB and has its own licence.

### Option A — let it download automatically (simplest)

Do nothing. The first time you run step 2, `transformers` downloads the model from Hugging
Face and caches it in your home folder (`~/.cache/huggingface`). You need internet for that
first run only.

### Option B — download it yourself (recommended if you will run offline or repeatedly)

```bash
pip install huggingface_hub
python -c "from huggingface_hub import snapshot_download; snapshot_download('Rostlab/prot_t5_xl_half_uniref50-enc', local_dir='models/prot_t5')"
```

That writes the model into `models/prot_t5/`. Then tell Salamander where it is, either per
command with `--plm models/prot_t5`, or once for the whole session:

**Windows (PowerShell)**
```powershell
$env:PLM = "models/prot_t5"
```

**macOS / Linux**
```bash
export PLM=models/prot_t5
```

**Check it loads:**
```bash
python -c "from transformers import T5EncoderModel; T5EncoderModel.from_pretrained('models/prot_t5'); print('ProtT5 OK')"
```

---

## 4. Get FoldX (required for steps 3 and 4)

FoldX computes the stability cost of each proposed mutation. It is commercial software, free
for academic use, and cannot be redistributed here.

1. Go to **[foldxsuite.crg.eu](https://foldxsuite.crg.eu/)** and register for an academic
   licence. Approval usually takes a day or two.
2. Download the build for your operating system and unzip it.
3. **Keep the `molecules/` folder next to the executable.** FoldX will not run without it.
   A correct install looks like:

```
C:/FoldX/                      (or ~/foldx/ on macOS/Linux)
├── foldx_YYYYMMDD.exe         the executable (name varies by release)
├── molecules/                 REQUIRED - do not move or rename
└── rotabase.txt               FoldX 4 only; FoldX 5 does not use it
```

4. Tell Salamander where it is:

**Windows (PowerShell)**
```powershell
$env:FOLDX = "C:/FoldX/foldx_20251231.exe"
```

**macOS / Linux**
```bash
export FOLDX=~/foldx/foldx
```

Or pass `--foldx C:/FoldX/foldx_20251231.exe` on each command. On Windows, if you skip both,
Salamander looks for the newest `.exe` in `C:/FoldX/` automatically.

**Check it runs:**
```bash
"C:/FoldX/foldx_20251231.exe" --version
```

> **If you have FoldX 4 and it reports an expired licence**, that build's licence ran out on
> 2025-12-31. You need FoldX 5. If you have both, delete or rename `rotabase.txt` — FoldX 5
> rejects the FoldX 4 rotabase file.

---

## 5. Get a structure for your protein

Steps 3 and 4 need a 3D structure whose **residue numbering matches your sequence exactly**.

**If your protein is in the PDB**: search [rcsb.org](https://www.rcsb.org/), download the
`.pdb` file. Be aware that crystal structures often start at residue 20-something and have
gaps — see the numbering note below.

**If it is not**: get a predicted structure from
[AlphaFold DB](https://alphafold.ebi.ac.uk/) by UniProt accession. AlphaFold models are
numbered 1..L with no gaps, which matches a plain FASTA sequence, so they are usually the
easier choice.

Put it anywhere, e.g. `structures/my_protein.pdb`. (The `structures/` folder is in
`.gitignore`, so nothing you put there gets committed.)

> **Numbering must match.** If your FASTA residue 45 is `V`, the PDB must have `VAL` at
> residue 45 of the chain you pass with `--chain`. Step 3 checks this before running anything
> and aborts with the exact mismatches listed if it fails. **This is the single most common
> reason the pipeline stops.** Fix it by using an AlphaFold model, or by trimming your FASTA
> to match the construct in the crystal structure.

---

## 6. Run the example end to end

The repo ships with the ligand-binding domain of human estrogen receptor α at
[`examples/era_lbd.fasta`](examples/era_lbd.fasta) (250 residues).

### 6a. Without a structure (steps 1–2 only)

This works immediately after install — no FoldX, no PDB needed:

```bash
python pipeline/step1_conserved_sites.py --query examples/era_lbd.fasta --out work/era
python pipeline/step2_design_set1star.py --work work/era
```

**What you get:**

```
work/era/
├── orthologs.fasta      the ortholog sequences OrthoDB returned
├── alignment.fasta      the multiple sequence alignment
├── conserved.json       which residues are frozen, and why
├── embedding.npz        cached ProtT5 embedding (delete to recompute)
└── set1star.json        the proposed mutations  <-- your result
```

Step 1 prints something like `conserved at >= 0.80 : 119 of 250 residues frozen, 131
mutable`. Step 2 prints a table of surviving mutations with their probabilities.

### 6b. Full run with FoldX

Download the ERα LBD structure — [PDB 1ERE](https://www.rcsb.org/structure/1ERE) or the
[AlphaFold model for P03372](https://alphafold.ebi.ac.uk/entry/P03372) — save it as
`structures/era.pdb`, then:

```bash
python pipeline/run_pipeline.py --query examples/era_lbd.fasta --pdb structures/era.pdb --out work/era --repair
```

That runs all four steps in order. Add `--foldx <path>` if you did not set `$FOLDX`.

**Final result:** `work/era/final.fasta` — the redesigned sequence.

---

## 7. Run it on your own protein

**Step 0.** Put your sequence in a FASTA file, one sequence only:

```
>MyProtein
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ...
```

Save as `my_protein.fasta`.

**Step 1 — freeze the conserved sites.**

```bash
python pipeline/step1_conserved_sites.py --query my_protein.fasta --out work/mine --cons-thresh 0.80 --n-orthologs 10
```

Finds orthologs via OrthoDB, aligns them, and marks conserved positions as untouchable.
Reads: your FASTA. Writes: `work/mine/conserved.json`, `alignment.fasta`, `orthologs.fasta`.

*No internet, or OrthoDB finds nothing?* Supply your own ortholog FASTA:

```bash
python pipeline/step1_conserved_sites.py --query my_protein.fasta --out work/mine --source local --family my_orthologs.fasta
```

**Step 2 — propose mutations.**

```bash
python pipeline/step2_design_set1star.py --work work/mine
```

Reads: `work/mine/conserved.json`. Writes: `work/mine/set1star.json`, `embedding.npz`.
Add `--plm models/prot_t5` if you downloaded ProtT5 manually and did not set `$PLM`.

**Step 3 — drop the destabilising ones.**

```bash
python pipeline/step3_foldx_filter.py --work work/mine --pdb structures/mine.pdb --repair --jobs 6
```

Reads: `work/mine/set1star.json` + your PDB. Writes: `work/mine/foldx_filter.json`.
Use `--repair` the first time. On reruns you can drop it — the repaired structure is cached
in `work/mine/foldx/`.

**Step 4 — assemble the final sequence.**

```bash
python pipeline/step4_final_sequence.py --work work/mine --combined
```

Reads: everything above. Writes: **`work/mine/final.fasta`** and `final.json`.

**Or do all four at once:**

```bash
python pipeline/run_pipeline.py --query my_protein.fasta --pdb structures/mine.pdb --out work/mine --repair
```

---

## 8. What every file means

### Files in this repository

| Path | What it is |
|---|---|
| [`models/head_strip2.pt`](models/head_strip2.pt) | **The trained model.** 5 MB. Loaded by step 2. |
| [`pipeline/step1_conserved_sites.py`](pipeline/step1_conserved_sites.py) | Ortholog search → MSA → frozen sites |
| [`pipeline/step2_design_set1star.py`](pipeline/step2_design_set1star.py) | Model inference → confident mutations |
| [`pipeline/step3_foldx_filter.py`](pipeline/step3_foldx_filter.py) | FoldX per-mutation stability filter |
| [`pipeline/step4_final_sequence.py`](pipeline/step4_final_sequence.py) | Assemble + score the whole design |
| [`pipeline/run_pipeline.py`](pipeline/run_pipeline.py) | Runs steps 1–4 in sequence |
| [`src/strip2/core.py`](src/strip2/core.py) | Embedding, checkpoint loading, positional encoding |
| [`src/strip2/orthologs.py`](src/strip2/orthologs.py) | OrthoDB and NCBI blastp retrieval |
| [`training/`](training/) | Scripts to retrain from scratch — see section 10 |
| [`training/data/combined_clustered_proteins.csv`](training/data/combined_clustered_proteins.csv) | 40 MB clustered, Tm-annotated training table |
| [`docs/method.md`](docs/method.md) | Benchmark results, the proline finding, caveats |
| [`examples/era_lbd.fasta`](examples/era_lbd.fasta) | Worked example sequence |

### Files the pipeline creates, in your `--out` folder

| File | Written by | Contains |
|---|---|---|
| `orthologs.fasta` | step 1 | the ortholog sequences that were found |
| `alignment.fasta` | step 1 | the MSA of query + orthologs |
| `conserved.json` | step 1 | per-residue conservation and the frozen list |
| `embedding.npz` | step 2 | cached ProtT5 embedding — delete to force recompute |
| `set1star.json` | step 2 | proposed mutations with probabilities |
| `foldx/` | step 3 | FoldX working directories (large, safe to delete) |
| `foldx_filter.json` | step 3 | per-mutation ΔΔG, kept vs dropped |
| **`final.fasta`** | step 4 | **your redesigned sequence** |
| `final.json` | step 4 | full report: every mutation, ΔΔG, and the summary |

---

## 9. Every option, explained

### Step 1 — `step1_conserved_sites.py`

| Option | Default | What it does |
|---|---|---|
| `--query` | *required* | FASTA with exactly one sequence |
| `--out` | *required* | output folder; created if missing |
| `--source` | `orthodb` | `orthodb`, `blast` (NCBI, slow), or `local` |
| `--family` | — | your own ortholog FASTA; required with `--source local` |
| `--n-orthologs` | `10` | how many orthologs to use |
| `--cons-thresh` | `0.80` | fraction agreeing before a column is frozen. **Higher = fewer frozen = more aggressive design** |
| `--min-identity` | `0.30` | reject orthologs below this identity |
| `--length-band` | `2.0` | keep orthologs within ×2 of query length |
| `--extra-freeze` | — | force-freeze residues you know matter, e.g. `"25 159 175"` (1-based) |
| `--email` | — | contact address, required by NCBI when using `--source blast` |

> With the default 10 orthologs, conservation is quantised to steps of 0.1 — so
> `--cons-thresh 0.80` and `0.85` give **identical** frozen sets. Raise `--n-orthologs` to
> 20–30 if you want finer control.

### Step 2 — `step2_design_set1star.py`

| Option | Default | What it does |
|---|---|---|
| `--work` | *required* | the folder from step 1 |
| `--model` | `models/head_strip2.pt` | model checkpoint |
| `--plm` | `$PLM`, else Hugging Face | ProtT5 path or HF id |
| `--p-min` | `0.10` | minimum P(mutant) to keep a mutation |
| `--wt-ratio` | `1.20` | P(mutant) must be ≥ this × P(wild type). **Raise for fewer, safer mutations** |
| `--emb-cache` | `<work>/embedding.npz` | where to cache the embedding |

### Step 3 — `step3_foldx_filter.py`

| Option | Default | What it does |
|---|---|---|
| `--work` | *required* | the folder from steps 1–2 |
| `--pdb` | *required* | structure matching the query numbering |
| `--foldx` | `$FOLDX`, else `C:/FoldX/*.exe` | FoldX executable |
| `--chain` | `A` | which chain in the PDB |
| `--ddg-cut` | `0.5` | keep mutations with ΔΔG ≤ this (kcal/mol) |
| `--repair` | off | run `RepairPDB` first — **use it on the first run** |
| `--jobs` | `6` | parallel FoldX processes; set to your core count |

> **On `--ddg-cut 0.5`**: this default is the one we validated. Across five proteins it kept
> ~53% of mutations and every design came out net stabilising (see
> [`docs/method.md`](docs/method.md)). FoldX's own error is around 1 kcal/mol, so a threshold
> much below 0.5 is inside the noise. Use `1.0` for a larger, bolder design.

### Step 4 — `step4_final_sequence.py`

| Option | Default | What it does |
|---|---|---|
| `--work` | *required* | the folder from steps 1–3 |
| `--combined` | off | also score the whole design as one multi-mutant |
| `--foldx` | `$FOLDX` | FoldX executable |
| `--name` | `<query>_strip2_design` | name for the FASTA record |

### `run_pipeline.py`

Accepts every option above and passes it through. Omit `--pdb` to stop after step 2.
`--no-combined` skips the final multi-mutant score.

---

## 10. Training from scratch (optional)

**You do not need this to use Salamander** — `models/head_strip2.pt` is ready to go. This
section is for reproducing or retraining the model.

Training runs in two phases and needs a GPU (a CPU run takes days).

### Phase A — the prokaryote base model

```bash
# A1  build training pairs, excluding yeast          (~30 min, CPU)
python training/build_pairs.py --csv training/data/combined_clustered_proteins.csv --work_dir runs/prok --exclude_taxid 559292

# A2  cache ProtT5 embeddings                        (hours, GPU strongly recommended)
python training/embed_sequences.py --work_dir runs/prok --plm models/prot_t5

# A3  train                                          (~2 h on a V100)
python training/train_prok.py --work_dir runs/prok --epochs 30
```

Produces `runs/prok/head.pt`.

### Phase B — the yeast transfer that makes strip2

```bash
# B1  yeast -> hotter-homolog pairs
python training/build_yeast_pairs.py --csv training/data/combined_clustered_proteins.csv --work_dir runs/yeast

# B2  embed the yeast sequences
python training/embed_sequences.py --work_dir runs/yeast --plm models/prot_t5

# B3  fine-tune prok -> strip2
python training/train_strip2_transfer.py --work_dir runs/yeast --pretrained runs/prok/head.pt --epochs 60 --out models/head_strip2.pt
```

**Why two phases.** Yeast Tm in this dataset is *organism-level* — every *S. cerevisiae*
protein carries the same value — so yeast-vs-yeast pairs have zero ΔTm and teach nothing.
`build_yeast_pairs.py` therefore pairs each yeast protein with a **hotter homolog from
another organism** in the same cluster.

**Why only one block is frozen.** strip2 is the `strip2_block+clf` configuration:

```
block 1     Linear(2064, 512)   FROZEN      keeps the general ProtT5 -> residue mapping
block 2     Linear(512, 512)    retrained
classifier  Linear(512, 20)     retrained
```

The yeast set is far smaller than the prokaryote one; letting the whole network move overfits
it. `--frozen-blocks 0` gives a full fine-tune, `--frozen-blocks 2` a classifier-only linear
probe (the "strip1" variant, val_acc 0.215 vs strip2's 0.256).

Expect validation accuracy around **0.24–0.26** against a majority-class baseline of ~0.10.
That looks low because it is genuine held-out-*family* generalisation — the split is by
cluster, so no homolog of a training protein appears in validation.

---

## 11. Troubleshooting

**`ModuleNotFoundError: No module named 'torch'`**
The virtual environment is not active. Run `.venv\Scripts\activate` (Windows) or
`source .venv/bin/activate`, then `pip install -r requirements.txt`.

**`FoldX not found - pass --foldx or set $FOLDX`**
Set the environment variable (section 4) or pass `--foldx /full/path/to/foldx`. Use the full
path including the filename, not the folder.

**`RepairPDB produced nothing - check the FoldX licence and rotabase`**
Either the licence expired (FoldX 4 licences ended 2025-12-31 — get FoldX 5), or the
`molecules/` folder is missing from beside the executable, or a FoldX 4 `rotabase.txt` is
present and FoldX 5 is rejecting it. Rename it to `rotabase.txt.bak`.

**`STRUCTURE MISMATCH - the PDB numbering does not match the query sequence`**
The most common failure. Step 3 lists the offending residues. Either use an AlphaFold model
(numbered 1..L with no gaps) or trim your FASTA to match the crystal construct. Also check
`--chain` — the default is `A`.

**`no set1* mutations to score - nothing to do`**
The filters removed everything. Lower `--wt-ratio` (try `1.10`), lower `--p-min` (try `0.05`),
or raise `--cons-thresh` in step 1 so fewer positions are frozen.

**`no orthologs found`**
OrthoDB had no match. Supply your own family with `--source local --family your_file.fasta`,
or try `--source blast --email you@example.com` (slow, several minutes).

**Step 2 is very slow / runs out of memory**
It is embedding on CPU. Either accept it (a few minutes for most proteins) or install a
CUDA-enabled PyTorch build from [pytorch.org](https://pytorch.org/get-started/locally/).

**`ImportError: pyfamsa`, or alignment falls back to pairwise**
`pip install pyfamsa`. The fallback works but produces a poorer MSA and therefore a less
reliable frozen set.

**Everything is frozen — `250 of 250 residues frozen`**
Your "orthologs" are nearly identical to the query, so every column looks conserved. Use more
distant orthologs, raise `--min-identity` to reject near-duplicates, or lower
`--cons-thresh`.

---

## 12. How it works

```
your sequence
      │
      ▼  step 1   orthologs → MSA → freeze conserved sites      protects the active site
      ▼  step 2   strip2 → 20-way softmax → set1* cut           keeps confident calls
      ▼  step 3   FoldX single-point ΔΔG → drop destabilising   keeps the fold intact
      ▼  step 4   assemble, score the whole design
redesigned sequence
```

**strip2** is a per-residue classifier on frozen ProtT5 embeddings. It learns from evolution,
not thermodynamics: given families of homologs where some members come from hot organisms and
some from cold, it learns which residue a *hotter* relative would carry at each position. No
measured ΔΔG or ΔTm is ever used as a training label.

Because it sees sequence but not structure, the model alone proposes too many mutations and
some are structurally destructive — it has learned that **proline substitutions are a
thermostabilising move** (real thermophiles carry them) but not *where* the backbone tolerates
one. In our benchmark the same substitution type sat at both extremes: **Gly→Pro at −2.85
kcal/mol** on ERα, **Gln→Pro at +4.60** on hDnmt3a. Step 3 is what separates those cases.

Full benchmark results, method comparisons and limitations are in
**[`docs/method.md`](docs/method.md)**.

### Limitations worth knowing

- **Positions are chosen independently**, so the model cannot design cooperativity between
  mutations the way a combinatorial Rosetta protocol can. On ERα, PROSS's 24 mutations sum to
  +5.28 kcal/mol as isolated singles yet score −0.60 assembled — about 5.9 kcal/mol gained
  purely from combination. Salamander has no mechanism for that.
- **ΔΔG is a proxy for the wrong observable.** The filters optimise folding stability; the
  model predicts thermostability. They correlate but are not the same, and nothing here
  verifies that predicted ΔTm survives filtering.
- **Multi-mutant ΔΔG degrades with mutation count.** FoldX relaxes from the wild-type
  backbone, so a design changing more than ~10% of the sequence increasingly measures "how
  badly do these residues fit the original backbone". Treat large designs as directional.
- **N-terminal bias.** The positional encoding cannot tell residue 1 of a 250-mer from residue
  1 of a 600-mer, and initiator methionine dominates that position in training. Check the
  first few residues of any design.
- **Nothing here is experimentally validated.** Every number in this repository is
  computational.

---

## Licence

MIT — see [`LICENSE`](LICENSE). ProtT5 and FoldX carry their own licences and are not
included in this repository.

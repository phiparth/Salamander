# Salamander

**Evolution-guided thermostable protein design.** Give it a protein sequence; it gives you back
a redesigned sequence predicted to melt at a higher temperature, with the active site left
alone.

*Named for the creature that was said to live in fire.*

> ### The structure step is not optional
>
> Every published sequence-based thermostability tool — TemStaPro, DeepSTABp, ProLaTherm,
> TemBERTure, ESMStabP, NOMELT — reads sequence and nothing else. **Salamander's contribution
> is what happens after the model speaks: every proposed mutation is scored against a real 3D
> structure in FoldX, and the destabilising ones are thrown away.**
>
> That is the whole point, and it is measurable. Before the structure filter, our designs are
> destabilising on every protein we tested. After it, **all five flip to net stabilising and
> four of the five beat PROSS** — see [Why the structure filter is the point](#why-the-structure-filter-is-the-point).
>
> So: **you need FoldX and you need a PDB.** Steps 1–2 will run without them, but what they
> emit is a candidate list, not a design. Do not treat it as a result.

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

> **Every code block below holds exactly one command.** Copy them one at a time, in order.

---

## 1. What you need before you start

| | What | Where to get it | Size | Needed for |
|---|---|---|---|---|
| ✅ | **Python 3.9 – 3.12** | [python.org/downloads](https://www.python.org/downloads/) | — | everything |
| ✅ | **Git** | [git-scm.com/downloads](https://git-scm.com/downloads) | — | cloning this repo |
| ✅ | **ProtT5 model** | [huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc](https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc) | ~2.5 GB | step 2, all training |
| ✅ | **FoldX 5** | [foldxsuite.crg.eu](https://foldxsuite.crg.eu/) — free academic licence | ~50 MB | steps 3 and 4 |
| ✅ | **A PDB structure** | step 0 fetches it for you, or [rcsb.org](https://www.rcsb.org/) | ~1 MB | steps 3 and 4 |
| ⬜ | GPU with CUDA | — | — | optional, ~10× faster embedding |
| ⬜ | Internet | — | — | steps 0 and 1 only |

**Disk space**: about 3 GB, almost all ProtT5.
**Time for one protein**: 2–5 min for steps 1–2, then 10–60 min for steps 3–4. `RepairPDB` in
step 3 is the slow part.

---

## 2. Install

Open a terminal — Command Prompt or PowerShell on Windows, Terminal on macOS/Linux.

Clone the repository:

```bash
git clone https://github.com/phiparth/Salamander.git
```

Enter the folder:

```bash
cd Salamander
```

Create an isolated Python environment so this does not disturb your other projects.

> **You need Python 3.9 or newer.** On Windows, do **not** type `python` for this step. Other
> software (MGLTools, ArcGIS, Cygwin, Anaconda) puts its own `python` on your PATH, and it is
> often Python 2.7 — which cannot run any of this. `py -3` asks Windows for a real Python 3.
> Check with `py -3 --version` before continuing.

> **The next four blocks are two pairs. Run the pair for YOUR operating system — not all four.**

**Windows** — create it:

```bash
py -3 -m venv .venv
```

**Windows** — then activate it:

```bash
.venv\Scripts\activate
```

**macOS / Linux** — create it:

```bash
python3 -m venv .venv
```

**macOS / Linux** — then activate it:

```bash
source .venv/bin/activate
```

Your prompt should now start with `(.venv)`. From here on, `python` means the one inside
`.venv`, which is the correct one. Upgrade pip:

```bash
python -m pip install --upgrade pip
```

Install the dependencies (a few minutes, mostly PyTorch):

```bash
pip install -r requirements.txt
```

Confirm it worked:

```bash
python check_env.py
```

This checks your Python version, your working directory, every required package, and the
ProtT5 weights, then prints `OK` or a `FAIL` line naming the fix. **Run it whenever anything
below goes wrong** — it reports the cause instead of leaving you with a traceback.

At this stage it is expected to pass sections 1–3 and report ProtT5 as missing; you install
that next.

If it reports `Python 2.7 is too old`, the virtual environment is not active: your prompt is
missing `(.venv)`. Run the `activate` line again, in this same terminal.

> **Every command from here on assumes you are inside the `Salamander` folder with `(.venv)`
> active.** If you close the terminal, `cd` back and re-run `activate`.

---

## 3. Get ProtT5 (required)

ProtT5 turns a sequence into the numeric embedding the model reads. It is **not** bundled here
— 2.5 GB, separate licence.

### Option A — automatic (simplest)

Do nothing. The first run of step 2 downloads it from Hugging Face into
`~/.cache/huggingface`. Internet needed for that first run only.

### Option B — download it yourself (better if you will run offline or often)

Install the downloader:

```bash
pip install huggingface_hub
```

Download the model into `models/prot_t5/`:

```bash
python -c "from huggingface_hub import snapshot_download; snapshot_download('Rostlab/prot_t5_xl_half_uniref50-enc', local_dir='models/prot_t5')"
```

Now tell Salamander where it is.

> **Run ONE of the next three blocks — whichever matches your shell. Not all three.**
> Prompt starts with `PS C:\>` → PowerShell. Prompt starts with `C:\>` → Command Prompt.

**Either — Windows, PowerShell:**

```powershell
$env:PLM = "models/prot_t5"
```

**Or — Windows, Command Prompt:**

```bash
set PLM=models/prot_t5
```

**Or — macOS / Linux:**

```bash
export PLM=models/prot_t5
```

Confirm it loads:

```bash
python check_env.py
```

The last section of the report lists every ProtT5 file with its size and then actually loads
the model. `pytorch_model.bin` must be roughly **2.2–2.5 GB**. If it is a few kilobytes the
download produced a Git LFS pointer instead of the weights — delete `models/prot_t5` and run
the download command again.

Alternatively skip the variable and pass `--plm models/prot_t5` on each command.

### Two version traps

ProtT5 is distributed as a `pytorch_model.bin` file — there is no `model.safetensors` in that
repository. That makes it sensitive to library versions in two ways, and `requirements.txt`
already pins around both. Only read this if you installed `transformers` or `torch` by hand.

- **`transformers` 5.x removed support for `.bin` checkpoints.** ProtT5 will not load on it.
  Stay on 4.x: `pip install "transformers>=4.30,<5"`
- **`transformers` 4.56+ refuses `.bin` files when `torch` is older than 2.6**, as a
  `torch.load` security measure. Either upgrade torch, or hold transformers below 4.56.

`python check_env.py` detects both combinations and tells you which one you have.

---

## 4. Get FoldX (required for steps 3 and 4)

FoldX computes the stability cost of each proposed mutation. Commercial software, free for
academic use, cannot be redistributed here.

### 4a. Register

Create an account: **[foldxsuite.crg.eu/user/register](https://foldxsuite.crg.eu/user/register)**

Use your **institutional email address**. Personal Gmail/Outlook addresses are rejected for
academic licences — the email domain is how they verify you are academic.

### 4b. Apply for the academic licence

Fill in the form at
**[foldxsuite.crg.eu/academic-license-info](https://foldxsuite.crg.eu/academic-license-info)**

Approval usually takes one to two days and arrives by email.

### 4c. Log in and download

Log in at **[foldxsuite.crg.eu/user/login](https://foldxsuite.crg.eu/user/login)**, then use the
**`Download → Academic License`** item in the top navigation bar.

That menu is under **Download**, *not* under "Licensing and Services" — the latter is only the
application form. The file list stays invisible until you are logged in with an approved
licence.

Pick the build for your operating system and take the newest **FoldX 5**. The filename encodes
a release date, e.g. `foldx_20251231.zip`.

### 4d. Unzip, keeping everything together

```
C:\FoldX\
├── foldx_20251231.exe      the executable (name varies by release)
├── molecules\              REQUIRED - do not move, rename, or delete
└── rotabase.txt            FoldX 4 only; FoldX 5 ignores it
```

**The `molecules/` folder is what people get wrong.** FoldX fails silently or produces nothing
if it is not sitting beside the executable.

### 4e. Verify it runs

```bash
"C:\FoldX\foldx_20251231.exe" --version
```

### 4f. Tell Salamander where it is

> **Run ONE of the next three blocks — whichever matches your shell. Not all three.**

**Either — Windows, Command Prompt** — permanent, reopen the terminal afterwards:

```bash
setx FOLDX "C:\FoldX\foldx_20251231.exe"
```

**Or — Windows, PowerShell** — current session only:

```powershell
$env:FOLDX = "C:/FoldX/foldx_20251231.exe"
```

**Or — macOS / Linux**:

```bash
export FOLDX=~/foldx/foldx
```

Or pass `--foldx <full path to executable>` on each command. On Windows, if you set neither,
Salamander picks the newest `.exe` in `C:/FoldX/` automatically.

> **Expired-licence error?** FoldX 4 licences ended 2025-12-31; you need FoldX 5. If you have
> both, **rename** `rotabase.txt` to `rotabase.txt.bak` — FoldX 5 rejects the FoldX 4 rotabase
> file. Rename rather than delete, in case you go back to FoldX 4.

---

## 5. Get a structure for your protein

Steps 3 and 4 need a 3D structure whose **residue numbering matches your sequence exactly**.
FoldX reads atomic coordinates; it cannot work from a sequence alone.

The structure does **not** have to be experimental — a prediction is fine, and is what the
published benchmark used.

### Option A — let step 0 do it (recommended)

If your protein is in UniProt, fetch the AlphaFold model by accession:

```bash
python pipeline/step0_get_structure.py --uniprot P03372 --query examples/era_lbd.fasta --out structures/era.pdb
```

This downloads the AlphaFold DB model, and — because AlphaFold covers the *whole* UniProt
entry while your query may be a domain — **locates your query inside the model, cuts it out,
and renumbers it 1..L**. It finishes by checking every residue against your FASTA and tells you
whether step 3 will accept it.

If your protein is not in UniProt, fold the sequence with ESMFold instead:

```bash
python pipeline/step0_get_structure.py --query my_protein.fasta --source esmfold --out structures/mine.pdb
```

> **ESMFold is best-effort.** It is a free public API with no account needed, roughly a
> 400-residue limit, and it is frequently overloaded — at the time of writing it returns
> HTTP 504 on new sequences. AlphaFold DB is the reliable route when your protein is in
> UniProt.

### Option B — get one yourself

Download from [rcsb.org](https://www.rcsb.org/) (experimental) or
[alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk/) (predicted, by UniProt accession) and save
it as e.g. `structures/mine.pdb`. The `structures/` folder is in `.gitignore`, so nothing you
put there gets committed.

> **Numbering must match.** If your FASTA residue 45 is `V`, the PDB must have `VAL` at residue
> 45 of the chain you pass to `--chain`. Step 3 checks this first and aborts with the exact
> mismatches listed. **This is the most common reason the pipeline stops.** Crystal structures
> often start at residue 20-something and contain gaps; AlphaFold models are numbered 1..L with
> no gaps, which is why they are usually easier.

---

## 6. Run the example end to end

The repo ships the ligand-binding domain of human estrogen receptor α at
[`examples/era_lbd.fasta`](examples/era_lbd.fasta) (250 residues).

### 6a. Install check — steps 1 and 2 only

Run this to confirm the install works. **It does not produce a usable design** — it stops
before the structure filter, so the mutation list still contains the structurally destructive
picks that step 3 exists to remove. Treat the output as a smoke test, then do 6b.

Freeze the conserved sites:

```bash
python pipeline/step1_conserved_sites.py --query examples/era_lbd.fasta --out work/era
```

Propose mutations:

```bash
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

Step 1 prints e.g. `conserved at >= 0.80 : 119 of 250 residues frozen, 131 mutable`. Step 2
prints a table of surviving mutations with probabilities.

> On ERα this stage proposes 25 mutations scoring **+1.67 kcal/mol — destabilising**. The
> structure filter cuts it to 7 mutations at **−6.65**. That gap is why 6a is not a result.

### 6b. Full run with FoldX — the real pipeline

Fetch the structure:

```bash
python pipeline/step0_get_structure.py --uniprot P03372 --query examples/era_lbd.fasta --out structures/era.pdb
```

Run all four steps:

```bash
python pipeline/run_pipeline.py --query examples/era_lbd.fasta --pdb structures/era.pdb --out work/era --repair
```

**Final result:** `work/era/final.fasta`.

---

## 7. Run it on your own protein

**Step 0 — put your sequence in a FASTA file**, one sequence only:

```
>MyProtein
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ...
```

Save it as `my_protein.fasta`.

**Get a structure:**

```bash
python pipeline/step0_get_structure.py --uniprot YOUR_ACCESSION --query my_protein.fasta --out structures/mine.pdb
```

**Step 1 — freeze the conserved sites.** Reads your FASTA; writes `conserved.json`,
`alignment.fasta`, `orthologs.fasta`:

```bash
python pipeline/step1_conserved_sites.py --query my_protein.fasta --out work/mine --cons-thresh 0.80 --n-orthologs 10
```

If OrthoDB finds nothing, or you are offline, supply your own ortholog FASTA instead:

```bash
python pipeline/step1_conserved_sites.py --query my_protein.fasta --out work/mine --source local --family my_orthologs.fasta
```

**Step 2 — propose mutations.** Reads `conserved.json`; writes `set1star.json`,
`embedding.npz`:

```bash
python pipeline/step2_design_set1star.py --work work/mine
```

**Step 3 — drop the destabilising ones. This is the step that makes the design real.** Reads
`set1star.json` + your PDB; writes `foldx_filter.json`:

```bash
python pipeline/step3_foldx_filter.py --work work/mine --pdb structures/mine.pdb --repair --jobs 6
```

Use `--repair` the first time only; the repaired structure is cached in `work/mine/foldx/`.

**Step 4 — assemble the final sequence.** Writes `final.fasta` and `final.json`:

```bash
python pipeline/step4_final_sequence.py --work work/mine --combined
```

**Or run steps 1–4 in one go:**

```bash
python pipeline/run_pipeline.py --query my_protein.fasta --pdb structures/mine.pdb --out work/mine --repair
```

---

## 8. What every file means

### Files in this repository

| Path | What it is |
|---|---|
| [`models/head_strip2.pt`](models/head_strip2.pt) | **The trained model.** 5 MB. Loaded by step 2. |
| [`pipeline/step0_get_structure.py`](pipeline/step0_get_structure.py) | Fetch or predict a structure |
| [`pipeline/step1_conserved_sites.py`](pipeline/step1_conserved_sites.py) | Ortholog search → MSA → frozen sites |
| [`pipeline/step2_design_set1star.py`](pipeline/step2_design_set1star.py) | Model inference → confident mutations |
| [`pipeline/step3_foldx_filter.py`](pipeline/step3_foldx_filter.py) | FoldX per-mutation stability filter |
| [`pipeline/step4_final_sequence.py`](pipeline/step4_final_sequence.py) | Assemble + score the whole design |
| [`pipeline/run_pipeline.py`](pipeline/run_pipeline.py) | Runs steps 1–4 in sequence |
| [`src/strip2/core.py`](src/strip2/core.py) | Embedding, checkpoint loading, positional encoding |
| [`src/strip2/orthologs.py`](src/strip2/orthologs.py) | OrthoDB and NCBI blastp retrieval |
| [`training/`](training/) | Retraining from scratch — see section 10 |
| [`training/data/combined_clustered_proteins.csv`](training/data/combined_clustered_proteins.csv) | 40 MB clustered, Tm-annotated training table |
| [`docs/method.md`](docs/method.md) | Benchmark results, the proline finding, caveats |
| [`examples/era_lbd.fasta`](examples/era_lbd.fasta) | Worked example sequence |

### Files the pipeline creates in your `--out` folder

| File | Written by | Contains |
|---|---|---|
| `orthologs.fasta` | step 1 | the ortholog sequences found |
| `alignment.fasta` | step 1 | MSA of query + orthologs |
| `conserved.json` | step 1 | per-residue conservation and the frozen list |
| `embedding.npz` | step 2 | cached ProtT5 embedding — delete to recompute |
| `set1star.json` | step 2 | proposed mutations with probabilities |
| `foldx/` | step 3 | FoldX working directories (large, safe to delete) |
| `foldx_filter.json` | step 3 | per-mutation ΔΔG, kept vs dropped |
| **`final.fasta`** | step 4 | **your redesigned sequence** |
| `final.json` | step 4 | full report: every mutation, ΔΔG, summary |

---

## 9. Every option, explained

### Step 0 — `step0_get_structure.py`

| Option | Default | What it does |
|---|---|---|
| `--out` | *required* | where to write the `.pdb` |
| `--query` | — | your FASTA; needed for `esmfold` and for the numbering check |
| `--uniprot` | — | UniProt accession, e.g. `P03372` |
| `--source` | `auto` | `auto`, `alphafold`, or `esmfold` |
| `--renumber` | off | shift numbering so the first residue is 1 |
| `--chain` | `A` | which chain to keep |

### Step 1 — `step1_conserved_sites.py`

| Option | Default | What it does |
|---|---|---|
| `--query` | *required* | FASTA with exactly one sequence |
| `--out` | *required* | output folder, created if missing |
| `--source` | `orthodb` | `orthodb`, `blast` (NCBI, slow), or `local` |
| `--family` | — | your own ortholog FASTA; required with `--source local` |
| `--n-orthologs` | `10` | how many orthologs to use |
| `--cons-thresh` | `0.80` | fraction agreeing before a column freezes. **Higher = fewer frozen = more aggressive** |
| `--min-identity` | `0.30` | reject orthologs below this identity |
| `--length-band` | `2.0` | keep orthologs within ×2 of query length |
| `--extra-freeze` | — | force-freeze residues you know matter, e.g. `"25 159 175"` (1-based) |
| `--email` | — | required by NCBI when using `--source blast` |

> With 10 orthologs, conservation is quantised to steps of 0.1 — so `--cons-thresh 0.80` and
> `0.85` give **identical** frozen sets. Raise `--n-orthologs` to 20–30 for finer control.

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

> **On `--ddg-cut 0.5`**: this default is the validated one. Across five proteins it kept ~53%
> of mutations and every design came out net stabilising — see
> [`docs/method.md`](docs/method.md). FoldX's own error is around 1 kcal/mol, so a threshold
> much below 0.5 sits inside the noise. Use `1.0` for a larger, bolder design.

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

**You do not need this to use Salamander** — `models/head_strip2.pt` is ready. This section
reproduces or retrains it. Two phases, and a GPU (CPU takes days).

### Phase A — the prokaryote base model

A1. Build training pairs, excluding yeast (~30 min, CPU):

```bash
python training/build_pairs.py --csv training/data/combined_clustered_proteins.csv --work_dir runs/prok --exclude_taxid 559292
```

A2. Cache ProtT5 embeddings (hours; GPU strongly recommended):

```bash
python training/embed_sequences.py --work_dir runs/prok --plm models/prot_t5
```

A3. Train (~2 h on a V100). Produces `runs/prok/head.pt`:

```bash
python training/train_prok.py --work_dir runs/prok --epochs 30
```

### Phase B — the yeast transfer that makes strip2

B1. Build yeast → hotter-homolog pairs:

```bash
python training/build_yeast_pairs.py --csv training/data/combined_clustered_proteins.csv --work_dir runs/yeast
```

B2. Embed the yeast sequences:

```bash
python training/embed_sequences.py --work_dir runs/yeast --plm models/prot_t5
```

B3. Fine-tune prok → strip2:

```bash
python training/train_strip2_transfer.py --work_dir runs/yeast --pretrained runs/prok/head.pt --epochs 60 --out models/head_strip2.pt
```

**Why two phases.** Yeast Tm in this dataset is *organism-level* — every *S. cerevisiae*
protein carries the same value — so yeast-vs-yeast pairs have zero ΔTm and teach nothing.
`build_yeast_pairs.py` therefore pairs each yeast protein with a **hotter homolog from another
organism** in the same cluster.

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

**First, run this. It diagnoses most of what follows:**

```bash
python check_env.py
```

**`ModuleNotFoundError: No module named 'transformers'` (or `'torch'`), or a traceback from
the ProtT5 check in section 3**

Almost always one of two things, and `check_env.py` distinguishes them:

1. **The virtual environment is not active.** Your prompt does not show `(.venv)`. Activation
   only applies to the terminal you ran it in — open a new tab and it is gone. Re-run the
   `activate` line from section 2 in *this* terminal.
2. **`python` is not Python 3.** On Windows, MGLTools, ArcGIS, Cygwin and others install their
   own `python` on your PATH, sometimes Python **2.7**, which fails on every command here.
   Check with `python --version`. If it is not 3.9+, use `py -3` in place of `python`, or
   activate the virtual environment, which fixes the name for that terminal.

**`No module named 'transformers'` with `ImportError` phrased in the old style**
That wording (rather than `ModuleNotFoundError`) is itself the tell: only Python 2 says
`ImportError: No module named x`. You are running Python 2.7. See point 2 above.

**`OSError: models/prot_t5 does not appear to have a file named config.json`**
Either the download did not finish, or you are not in the `Salamander` folder — `models/prot_t5`
is a *relative* path, so it only resolves from the repo root. `check_env.py` reports both.

**`FoldX not found - pass --foldx or set $FOLDX`**
Set the environment variable (section 4f) or pass `--foldx <path>`. Use the full path
*including* the executable filename, not just the folder.

**`RepairPDB produced nothing - check the FoldX licence and rotabase`**
Either the licence expired (FoldX 4 ended 2025-12-31 — get FoldX 5), or `molecules/` is missing
from beside the executable, or a FoldX 4 `rotabase.txt` is present and FoldX 5 is rejecting it.
Rename it to `rotabase.txt.bak`.

**`STRUCTURE MISMATCH - the PDB numbering does not match the query sequence`**
The most common failure. Step 3 lists the offending residues. Use
`pipeline/step0_get_structure.py`, which extracts and renumbers automatically. Also check
`--chain` — the default is `A`.

**`no set1* mutations to score - nothing to do`**
The filters removed everything. Lower `--wt-ratio` (try `1.10`), lower `--p-min` (try `0.05`),
or raise `--cons-thresh` in step 1 so fewer positions are frozen.

**`no orthologs found`**
OrthoDB had no match. Supply your own family with `--source local --family your_file.fasta`, or
try `--source blast --email you@example.com` (slow, several minutes).

**ESMFold returns HTTP 504**
The public ESM Atlas API is overloaded. Use `--source alphafold --uniprot <accession>` instead,
or download a structure by hand.

**Step 2 is very slow, or runs out of memory**
It is embedding on CPU. Either accept it (a few minutes for most proteins) or install a
CUDA-enabled PyTorch from [pytorch.org](https://pytorch.org/get-started/locally/).

**Everything is frozen — `250 of 250 residues frozen`**
Your "orthologs" are nearly identical to the query, so every column looks conserved. Use more
distant orthologs, raise `--min-identity` to reject near-duplicates, or lower `--cons-thresh`.

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
kcal/mol** on ERα, **Gln→Pro at +4.60** on hDnmt3a. Step 3 separates those cases.

Full benchmark results, method comparisons and limitations are in
**[`docs/method.md`](docs/method.md)**.

### Why the structure filter is the point

Every published sequence-based thermostability method stops at the model's output. Salamander
does not, and this is the measurable consequence.

**Before and after the FoldX filter**, five human proteins, strip2 designs at the OrthoDB ≥80%
freeze, FoldX combined ΔΔG in kcal/mol (negative = stabilising):

| protein | before | after | mutations kept | PROSS | beats PROSS |
|---|---|---|---|---|---|
| ERα-LBD | +1.67 | **−6.65** | 7 of 25 | −0.60 | **yes, 11×** |
| TPH1 | +13.32 | **−2.51** | 9 of 22 | −1.69 | **yes** |
| hAChE | +52.69 | **−14.44** | 29 of 58 | −13.34 | **yes**, with 43% fewer mutations |
| hSIRT6 | +9.86 | **−5.39** | 14 of 22 | −1.05 | **yes, 5×** |
| hDnmt3a | +9.93 | **−0.36** | 6 of 11 | −2.49 | no |

Every unfiltered design is destabilising. Every filtered one is stabilising. **Per-mutation
cost — which is count-normalised, so it cannot be improved just by making fewer changes — goes
negative in all five.**

**Why the model cannot do this alone.** strip2 learns from evolution, so it correctly learns
that proline substitutions are a thermostabilising move: real thermophiles carry them. What
sequence data cannot teach is *where the backbone tolerates one*. The identical substitution
type sat at both extremes of our single-point measurements:

```
Gly -> Pro    ERa      G138P    -2.85 kcal/mol     the best mutation in the whole benchmark
Gln -> Pro    hDnmt3a  Q52P     +4.60 kcal/mol     one of the worst
```

No sequence model can separate those two. A structure can, in seconds. That is the niche.

**One honest caveat**, also stated in [`docs/method.md`](docs/method.md): the mutations are
selected using FoldX and then measured with FoldX, so part of that improvement is guaranteed by
construction. An independent energy function should be used to confirm any specific design —
Rosetta `ref2015_cart` agreed with FoldX at Spearman 0.942 on the unfiltered benchmark, so it is
the natural check.

### Limitations worth knowing

- **Positions are chosen independently**, so the model cannot design cooperativity between
  mutations the way a combinatorial Rosetta protocol can. On ERα, PROSS's 24 mutations sum to
  +5.28 kcal/mol as isolated singles yet score −0.60 assembled — about 5.9 kcal/mol gained
  purely from combination. Salamander has no mechanism for that.
- **ΔΔG is a proxy for the wrong observable.** The filters optimise folding stability; the model
  predicts thermostability. They correlate but are not the same, and nothing here verifies that
  predicted ΔTm survives filtering.
- **Multi-mutant ΔΔG degrades with mutation count.** FoldX relaxes from the wild-type backbone,
  so a design changing more than ~10% of the sequence increasingly measures "how badly do these
  residues fit the original backbone". Treat large designs as directional.
- **N-terminal bias.** The positional encoding cannot tell residue 1 of a 250-mer from residue 1
  of a 600-mer, and initiator methionine dominates that position in training. Check the first
  few residues of any design.
- **Nothing here is experimentally validated.** Every number in this repository is
  computational.

---

## Licence

MIT — see [`LICENSE`](LICENSE). ProtT5 and FoldX carry their own licences and are not included
in this repository.

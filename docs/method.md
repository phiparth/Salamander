# Method and benchmark results

Background for anyone deciding how much to trust an output of this pipeline.

## Why three filters and not one

strip2 alone, applied at every mutable position, is far too aggressive. On hAChE it proposes
123 mutations; on Munc18a, 56. Those designs are badly destabilising by every energy function
we tried. Each filter removes a different failure mode:

| filter | removes | mechanism |
|---|---|---|
| conservation | active sites, fold nucleation | family MSA, modal-residue fraction |
| set1\* | low-confidence guesses | model must prefer the mutant to wild type by a margin |
| FoldX ΔΔG | structurally destructive picks | isolated point-mutation energy |

## What each filter is worth

Measured on five human proteins, design = strip2 at the OrthoDB ≥80% freeze.

**Confidence filtering (set1\*, 1.20× ratio)** cuts mutation counts roughly in half and
improves total ΔΔG in every case. But per mutation it only improves 7 of 16 cases — most of
the total gain comes from making *fewer* changes, not better ones. Do not overclaim it.

**Stability filtering** is where the large gains are. Dropping mutations above +0.5 kcal/mol:

```
protein        n   cut  keep    ddG all  ddG pruned   /mut all  /mut pruned
ERa-LBD       10     3     7      -2.47       -6.65      -0.25        -0.95
TPH1          22    13     9      13.32       -2.51       0.61        -0.28
hAChE         58    29    29      52.69      -14.44       0.91        -0.50
hSIRT6        22     8    14       9.86       -5.39       0.45        -0.39
hDnmt3a       11     5     6       9.93       -0.36       0.90        -0.06
```

All five flip from destabilising to stabilising, and per-mutation cost — which is
count-normalised and so cannot be improved just by making fewer changes — goes negative in
every case.

**Caveat**: those numbers select mutations using FoldX and then measure the result with
FoldX. Part of the improvement is guaranteed by construction. An independent method should
be used to confirm any specific design.

## Why the mutation-count confound matters

Across the whole benchmark, total ΔΔG correlates with mutation count at Spearman **+0.94**
*within* each protein. So a ranking of designs by total ΔΔG is largely a ranking by how many
mutations they make. Always report **ΔΔG per mutation** alongside the total.

The count-independent result that does survive: PROSS's per-mutation cost is negative on 7 of
8 proteins, while unfiltered strip2 is positive on nearly all. PROSS's individual
substitutions pay for themselves; strip2's, before filtering, do not.

## Method agreement

The same 48 designs scored by two independent energy functions:

| pair | Spearman |
|---|---|
| FoldX vs Rosetta `ref2015_cart` | **+0.942** |
| DynaMut2 vs FoldX, ≤20 mutations | +0.900 |
| DynaMut2 vs FoldX, >20 mutations | 0.000 |
| DDMut | unusable — hard cap at 33 mutations |

FoldX and Rosetta agree closely. DynaMut2 saturates above ~20 mutations and returns ~1
kcal/mol regardless of design quality, so it is only a useful check on small designs. Rosetta
carries its own bias here: PROSS *is* a Rosetta protocol, so scoring PROSS designs with
Rosetta flatters them — its advantage widened from −0.60 to −19.15 REU on ERα.

## The proline signature

The clearest mechanistic finding, and the reason step 3 recovers so much. Proline
substitutions appear at **both extremes** of the single-point distribution:

```
best     ERa      G138P   -2.85     (Gly→Pro: most flexible residue → most constrained)
         hAChE    S310P   -2.24
         hAChE    S68P    -1.80
worst    hDnmt3a  Q52P    +4.60
         hDnmt3a  Q30P    +3.77
         TPH1     E102P   +1.95
```

Proline rigidifies the unfolded state and is a textbook thermostabilising move — real
thermophiles carry them, so the training data teaches it. What the sequence data cannot teach
is *where the backbone tolerates one*. Forced into a helix or strand, proline is severely
destabilising. The single-point screen is what tells the two cases apart.

## Where the training signal comes from

An example from the real training set (cluster 369): a 20 °C protein and a 73 °C homolog,
ΔTm = 53. Their shared MSA is 188 columns; the cool protein occupies columns 117–137.

- **13 columns → "mutate"** — the hot protein has a different residue
- **4 columns → "keep"** — both already agree
- **4 columns → no label** — 100% conserved across all seven members, frozen

The frozen columns anchor the alignment register. At column 126 every member has `G`
(`GGGGGGG`), which is why the `N→R` assignment at column 127 is unambiguous — there are no
gaps in either sequence across that window, so nothing can slide.

Two limits of this scheme worth knowing: alignment errors are **correlated within a cluster**
(one bad MSA corrupts every pair drawn from it, which is why the split is by cluster), and
gap-placement ambiguity in low-identity regions is a genuine source of label noise that
nothing in the pipeline detects.

# Native-output reconstruction of the BHLHE40 Hi-TrAC boundary-deletion locus

**Status:** complete  
**Date:** 11 August 2026  
**Reference figure:** Hi-TrAC NAR Figure 3D  
**Locus:** hg38 `chr3:4,803,502-5,144,387` (1-based inclusive)

## Executive conclusion

Seven official de novo chromatin-contact predictors were executed successfully
on WT and inferred BHLHE40-boundary-deletion DNA. Every output was retained at
the highest resolution and in the native geometry and scientific scale of the
released checkpoint. No model output was converted into absolute Hi-TrAC PET
counts, and no single cross-model PCC was manufactured from incompatible target
normalizations.

All seven models produced a nonzero deletion response. That establishes that the
inferred sequence perturbation reached every inference path. The benchmark is
now explicitly framed as **cross-cell-type transfer to K562**: K562 sequence and,
where supported, K562 DNase are supplied to released checkpoints without
retraining. It does **not** rank their biological accuracy: only DeepC has a
K562-matched contact checkpoint; the remaining heads are HFF/H1/other-cell-type
or generic-human releases.

## 1. Experimental reference and perturbation

The published Hi-TrAC study pooled two replicates per condition and downsampled
control and deletion libraries to 33 million cis unique PETs for Figure 3D. The
local experimental panels use the registered pooled WT/deletion 1-kb matrices.

The paper reports two guide sequences but not the sequence-resolved junction of
the assayed clone. A new primary-source audit covered the complete 15-page
supplement, all nine supplemental workbooks, the authors' supplemental GitHub
repository and GEO GSE208085. None reports a clone junction. Independent exact
mapping of both guides to hg38 found one unique reverse-strand hit for each,
valid NGG PAMs, and expected SpCas9 cuts at 0-based coordinates 4,976,067 and
4,976,790. The preregistered interval is therefore exactly the implied clean
cut-to-cut deletion:

```text
chr3:4,976,067-4,976,790  (0-based, half-open)
length: 723 bp
```

It is therefore labelled **inferred clean SpCas9 cut-to-cut deletion; exact
clone junction unavailable**. This corrects the earlier, less precise phrase
"sgRNA-bounded deletion": the guides span 743 bp, while the inferred deletion
between their cleavage sites is 723 bp.

The published K562 active sub-TAD containing BHLHE40 is
`chr3:4,977,047-5,056,047` in Supplemental Table 2. Supplemental Table 4 records
the boundary-deletion read libraries but not the molecular breakpoint. Full
guide, PAM, cut-site and source details are recorded in
`results/DELETION_BREAKPOINT_REAUDIT.md`.

For fixed-length model inputs, the 723 bp were deleted, the flanks joined, and
723 bp of downstream reference sequence appended. Positions downstream of the
junction are derivative coordinates; their corresponding hg38 position is
offset by +723 bp. WT/deletion difference panels are consequently aligned by
model-output index on the derivative chromosome, the natural coordinate system
for sequence-perturbation inference.

## 2. Input contract

### Primary strict de novo arm

- DNA-only models receive WT or derivative DNA.
- Accessibility-aware models receive control K562 DNase in both conditions.
- In the deletion condition, the control signal moves with its attached DNA
  across the derivative junction; it is not replaced with deletion Hi-TrAC.
- For ChromaFold, the official AH104727 motif intervals are transformed onto
  the derivative chromosome before exact 50-bp rasterization. Four motifs
  overlapping deleted bases are removed and 1,463 intact downstream motifs in
  the registered input domain shift by exactly 723 bp. An explicit scan of all
  sequence windows crossing the new junction with FIMO 5.5.9 against the same
  three JASPAR 2022 CTCF PWMs found no hit at `p <= 1e-4`.

The native input signature of every released model is respected. DNase is not
injected into AkitaV2, DeepC, Orca, AlphaGenome or Chimaera because those
inference interfaces do not accept it. Their non-K562/generic output heads are
therefore sequence-only checkpoint-transfer tests. EPCOT and ChromaFold accept
cell-state accessibility information and are evaluated with K562 DNase; they
are accessibility-conditioned cross-cell-transfer tests. DeepC is the sole
K562-checkpoint reference rather than a transfer case.

### Secondary target-assisted sensitivity arm

For EPCOT and ChromaFold only, the experimental deletion Hi-TrAC 1D endpoint
profile was rank/quantile matched to the empirical control-DNase distribution.
The resulting 1-kb surrogate had essentially the same total and mean as the
control-DNase track and Spearman correlation 1.0 with the deletion endpoint
ranking.

This arm uses target-derived information. It is **not de novo**, is not a native
input contract for either model, and is excluded from model ranking. It answers
only whether those architectures respond when supplied with an accessibility-
shaped encoding of the experimental 1D perturbation.

## 3. Model-specific execution

| Model | Released checkpoint/channel used | Native geometry and scale | Native K562 input | Transfer interpretation |
|---|---|---|---|---|
| AkitaV2 | fold-2 human, HFF channel | 512-bin, 2,048-bp processed log(O/E); diagonal offsets 0-1 masked | DNA | Sequence-only HFF-head transfer |
| DeepC | K562 5-kb checkpoint | 83 center poles x 201 partner offsets; normalized skeleton profile | DNA | **K562-matched reference** |
| Orca | 1-Mb HFF checkpoint | 250 x 250 at 4 kb; normalized contact enrichment | DNA | Sequence-only HFF-head transfer |
| EPCOT | HFF Micro-C 1-kb head | 500 x 500 upper-triangle map; model-native predicted O/E | DNA + K562 DNase | Accessibility-conditioned transfer |
| ChromaFold motif | motif/no-coaccess checkpoint | 36 centers x 400 partners at 10 kb; HiC-DC+ Z-score | K562 DNase proxy + hg38 motif | Accessibility-conditioned transfer; bulk DNase substitutes for scATAC |
| AlphaGenome | HFFc6 Micro-C track | 512 x 512 at 2,048 bp; relative contact output | DNA | Sequence-only HFFc6-head transfer |
| Chimaera | official human release | four 32 x 128 rotated distance-coordinate images | DNA | Sequence-only generic-human transfer |

For multi-output releases, the overview uses HFF-like Micro-C channels to avoid
silently averaging incompatible cell types. These are deliberately retained as
unmatched transfer heads, not described as K562 predictions. The complete
AlphaGenome metadata table records all 28 returned contact tracks.

## 4. Native perturbation response

These statistics compare WT and deletion **within the same model and native
scale** over the displayed locus. `RMS delta / WT SD` is a scale-normalized
measure of perturbation magnitude; it is not accuracy against Hi-TrAC.

| Output | WT-vs-deletion PCC | RMS delta / WT SD |
|---|---:|---:|
| Experimental Hi-TrAC, log1p PET | 0.6648 | 0.820 |
| AkitaV2 HFF | 0.9837 | 0.257 |
| DeepC K562 | 0.9748 | 0.253 |
| Orca HFF | 0.9716 | 0.277 |
| EPCOT HFF + K562 DNase | 0.9138 | 0.479 |
| ChromaFold motif + K562 DNase proxy | 0.7556 | 0.805 |
| AlphaGenome HFFc6 Micro-C | 0.9845 | 0.177 |
| Chimaera human | 0.9229 | 0.386 |

The predicted deletions generally preserve more of the WT field than the
experimental deletion does. ChromaFold exhibits the largest normalized response,
but this cannot be interpreted as best accuracy because it uses an adapted bulk-
DNase proxy, a cross-cell checkpoint and a HiC-DC+ Z-score target. A fair ranking
would require converting experimental truth independently into each model's
training target and evaluating genome-wide held-out loci, not fitting a transform
on this single case study.

## 5. What the figures show

1. `FIGURE3D_LOCUS_AND_INPUT_TRACKS.png` registers genes, the inferred deletion,
   control DNase, derivative-shifted control DNase and the non-de-novo surrogate.
2. `FIGURE3D_DENSE_NATIVE_OUTPUTS.png` shows experimental Hi-TrAC followed by
   AkitaV2, Orca, EPCOT and AlphaGenome square outputs.
3. `FIGURE3D_NON_SQUARE_NATIVE_OUTPUTS.png` preserves DeepC poles, ChromaFold
   V-stripes and Chimaera rotated images rather than pretending they are the same
   square-matrix object.
4. `FIGURE3D_ALL_NATIVE_OUTPUTS.png` is the single overview panel.
5. `FIGURE3D_TARGET_ASSISTED_NOT_DENOVO.png` isolates the target-assisted EPCOT
   and ChromaFold sensitivity analysis.
6. `GVIZ_K562_HIC_BHLHE40_FIGURE3D.png` adds a Gviz hg38 chromosome-3 ideogram,
   Gviz locus annotations, a one-dimensional contact marginal and the exact
   paper locus extracted from the user's historical `K562.hic` at native 1-kb
   resolution using Juicer `observed NONE`.
7. `FIGURE3D_GEOMETRY_STANDARDIZED_SQUARES.png` represents all seven predictors
   as square endpoint-by-endpoint matrices over the same locus. It complements,
   rather than replaces, the native-geometry figures.

WT and deletion panels share a color scale within each row. Difference panels
use a row-specific symmetric scale. Color intensity must not be compared across
rows because every model predicts a different transformed quantity.

### 5.1 Geometry-standardized square comparison

The square comparison follows C.Origami Supplementary Figure 8 and the
associated 2026 bulk Hi-C benchmark, but deliberately standardizes **geometry
only**. It does not reproduce C.Origami's distance-stratified target
normalization or 128-bin image resizing because either would alter the model's
native scientific values.

- DeepC's 201-value center poles are mirrored and placed at their genomic
  endpoint pairs, which is coordinate-equivalent to the supplement's
  mirror/45-degree-rotation/crop workflow. Duplicate left/right estimates are
  averaged. The resulting displayed matrix is complete at 5 kb.
- ChromaFold's 400-value V-stripes are expanded into genomic pairs and duplicate
  estimates are averaged exactly as in `CBIGR/bulk_hic_benchmark`'s
  `bedpe_norm.py`. The intentionally unpredicted self-diagonal remains missing.
- Chimaera's midpoint-distance pixels are inverted to endpoint pairs. Values
  from overlapping fragments are averaged, while separations outside the
  released 8-135-kb output band remain missing rather than being imputed.
- AkitaV2 and EPCOT upper triangles are mirrored; Orca and AlphaGenome are
  already dense square matrices.

All text uses Arial. WT and deletion share one exact native-scale colorbar in
each row; the difference has a symmetric row-specific scale. Predictor color
limits are the exact finite minimum and maximum after coordinate assembly, not
quantiles. The experimental row alone uses the paper/cLoops2 display transform
`log10(PET + 1)` with fixed 0.1-0.8 clipping, white at the minimum and red at the
maximum. The raw experimental PET matrices are also retained in the archive.

The reconstructed square shapes are: experimental 342 x 342 at 1 kb, AkitaV2
166 x 166 at 2,048 bp, DeepC 68 x 68 at 5 kb, Orca 86 x 86 at 4 kb, EPCOT
340 x 340 at 1 kb, ChromaFold 34 x 34 at 10 kb, AlphaGenome 166 x 166 at
2,048 bp and Chimaera 166 x 166 at 2,048 bp. These differing shapes are a
faithful consequence of preserving native resolution, not a failure to align
the genomic locus.

## 6. Execution corrections and failures recorded

- Corrected Orca's output label from 1 kb to the observed native 4 kb.
- Corrected ChromaFold to a 4.01-Mb input and 400-partner V-stripe.
- Corrected Chimaera's shape convention to 32 distance rows x 128 position bins.
- Replaced an invalid fixed-coordinate deletion with a true derivative sequence:
  delete, join, shift attached tracks, then append downstream tail.
- Akita required isolated `tf-keras==2.16.0` compatibility; checkpoint/source
  weights were unchanged.
- DeepC required the environment's packaged CUDA libraries on its runtime path;
  model weights and inference were unchanged.
- EPCOT's released checkpoint contains 51 obsolete decoder/pretraining keys. The
  official cross-cell loader also filters to current model keys; the audited run
  reproduced that behavior and required zero missing current keys.
- ChromaFold initially failed because a SciPy sparse slice was passed where a
  dense array was required; conversion with `.toarray()` fixed data handling
  without altering the model.
- A publication audit then found that the first deletion adapter sampled the
  released 50-bp motif track at shifted bin centers. Because 723 is not divisible
  by 50, this displaced downstream motif bins and could not represent a new
  junction motif. That output is preserved under
  `results/native/chromafold/superseded_center_sampled/` but is ineligible.
  The corrected adapter starts from official AH104727 motif intervals, removes
  disrupted motifs, shifts intact downstream motifs at base-pair resolution and
  then repeats ChromaFold's 50-bp rasterization. It reproduces the released WT
  motif track over 87,200 input bins exactly (maximum difference 0), changes
  1,045 deletion-track bins relative to the superseded shortcut, and gives
  bitwise-identical predictions on an independent repeated inference run.
- All seven consolidated outputs passed shape, finiteness and SHA-256 audits.
- Re-audited every supplemental workbook and the complete supplemental PDF;
  corrected the deletion terminology from guide-bounded to SpCas9 cut-to-cut.
- Re-extracted the exact Figure 3D locus from the historical `K562.hic` without
  KR/VC balancing or O/E transformation. The 342 x 342 enclosing 1-kb grid has
  55,196 upper-triangle contacts and 8,714 nonzero upper-triangle cells.

## 7. Scientific limitations

1. The inferred allele may differ from the clone's actual junction or local
   repair sequence.
2. Six of seven primary output heads are not K562 contact-map checkpoints.
3. Bulk DNase is not native pseudobulk scATAC for ChromaFold, and no coaccessibility
   input was used because the selected released checkpoint excludes it.
4. AlphaGenome and Akita/Orca channels are selected explicitly, not averaged.
5. These normalized outputs cannot recover absolute Hi-TrAC PET depth.
6. This is a mechanistic locus case study, not a genome-wide held-out benchmark.
7. The historical `K562.hic` locus panel is not the pooled/downsampled 33-million
   cis-PET control used in the published Figure 3D; it is a separately prepared
   K562 contact map shown at the user's request.

## 8. Reproducibility and integrity

- Common input bundle SHA-256: `6040cf45f0b75005a6a18c4e1f5e19da493ccbad5a1e476bbd87e0c14ea865bb`
- Control DNase SHA-256: `7ba6ba20ed8c25acc282d03253360ac8edcdce6b9ab1a90832a980cfeaad12c6`
- Experimental truth SHA-256: `941ec6a9803e8769b2bb6dc7fb690f57b71769ebd63c58a7fe532cb8e22e6aab`
- Corrected ChromaFold native output SHA-256: `206ed50282965c16c89ce2edeb2ed39fdc6f37c9b3e8eb9c709c07e1d92ec3af`
- Official AH104727 motif resource SHA-256: `302251750af6fc23e23a783626e3b5f6ec111704fe4bcdc381dc71c6eb9bb8df`

Per-model output and checkpoint hashes are in the machine-readable audits under
`results/native/` and `results/MODEL_INPUT_OUTPUT_AUDIT.csv`.

No colleague process was stopped, restarted, reprioritized or modified. All
inference jobs completed; no training occurred.

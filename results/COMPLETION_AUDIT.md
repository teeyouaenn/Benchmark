# Completion audit

**Verdict: PASS**

- [x] Correct hg38 Figure 3D locus registered.
- [x] Two published sgRNAs independently mapped to unique hg38 reverse-strand protospacers with valid NGG PAMs.
- [x] Inferred 723-bp SpCas9 cut-to-cut deletion explicitly distinguished from an unknown clone junction.
- [x] Complete 15-page supplement and all nine supplemental workbooks searched for clone/breakpoint evidence.
- [x] Authors' supplemental repository and GEO GSE208085 checked; no clone-resolved junction found.
- [x] Fixed-length derivative chromosome constructed by delete/join/tail extension.
- [x] Control K562 DNase used for the primary deletion arm and shifted with attached DNA.
- [x] Cross-cell-transfer contract explicit: K562 DNase used only by architectures that natively accept accessibility; DNA-only models were not altered.
- [x] DeepC K562 checkpoint separated as the K562-matched reference rather than labelled transfer.
- [x] Target-derived Hi-TrAC 1D arm isolated and labelled non-de-novo.
- [x] AkitaV2 official fold-2 checkpoint complete.
- [x] DeepC official K562 5-kb checkpoint complete.
- [x] Orca official HFF and H1 1-Mb checkpoints complete; 4-kb native output verified.
- [x] EPCOT official HFF Micro-C 1-kb head complete.
- [x] ChromaFold official motif/no-coaccess checkpoint complete using motif-level derivative reconstruction.
- [x] Official AH104727 motif intervals reproduce the released WT 50-bp track exactly over all 87,200 inference bins.
- [x] Exact FIMO 5.5.9 scan of the new derivative junction finds no CTCF motif at p <= 1e-4.
- [x] ChromaFold's superseded center-sampled deletion adapter preserved and excluded; corrected inference reproduced bitwise on repeat.
- [x] AlphaGenome official all-fold local weights complete; 28 contact tracks returned.
- [x] Chimaera official human release complete.
- [x] All expected output shapes verified.
- [x] Every consolidated prediction finite.
- [x] Twenty-one active native artifacts are registered in the regenerated SHA-256 manifest; two superseded ChromaFold artifacts are preserved separately.
- [x] Seven-model machine-readable audit reports PASS.
- [x] Native dense, non-square, combined and target-assisted figures rendered.
- [x] DeepC center poles reconstructed as a square through the C.Origami Supplementary Figure 8 mirror/rotate/crop-equivalent coordinate mapping.
- [x] ChromaFold V-stripes reconstructed as a square using the associated benchmark's pair expansion and duplicate averaging.
- [x] Chimaera midpoint-distance outputs inverted to square endpoint geometry; unsupported distances remain missing.
- [x] No predictor values normalized, calibrated, quantile-clipped or fitted during square reconstruction.
- [x] Experimental Hi-TrAC rendered with `log10(PET + 1)`, fixed 0.1-0.8 paper scale, white-to-red palette; raw PET matrices retained.
- [x] Arial applied to every label in the geometry-standardized figure.
- [x] Exact per-model native and difference color scales displayed and recorded.
- [x] Geometry-standardized square archive, scale table, JSON audit and SHA-256 manifest generated.
- [x] Gviz chromosome-3 ideogram and locus tracks rendered at the exact Figure 3D interval.
- [x] Historical K562.hic extracted with Juicer observed NONE at 1 kb; 342 x 342 grid and 55,196 upper-triangle contacts independently verified.
- [x] Scientific report records failures, modifications, cell-type mismatches and target-scale caveats.
- [x] No absolute-count or cross-native-scale performance claim made.
- [x] No colleague/system process stopped or modified.

The benchmark is complete as a native-output perturbation visualization. It is
not a matched genome-wide performance ranking.

# BHLHE40 deletion-breakpoint re-audit

**Date:** 11 August 2026  
**Verdict:** the registered 723-bp interval is correct as an **inferred clean
SpCas9 cut-to-cut deletion**. It is not proven to be the exact clone junction.

## Sources exhaustively checked

1. NAR main article and Methods, including Figure 3D.
2. The complete 15-page supplemental PDF.
3. Every worksheet in all nine supplemental Excel files:
   - ENCODE metadata;
   - GM12878/K562 Hi-TrAC sub-TADs;
   - GM12878 super-enhancer components;
   - K562 boundary-deletion and CTCF-AID read summary;
   - cell-specific Hi-TrAC sub-TADs;
   - CTCF/RAD21 knock-down lost sub-TADs;
   - mouse Th17 read summary;
   - mouse Th17 peaks/segmentations/super-enhancers;
   - mouse Th17 active and Mll4-dependent sub-TADs.
4. The authors' supplemental GitHub repository at commit
   `80ee5589062e207d4ca097af3525890a14ddedd4`.
5. GEO series GSE208085 metadata and file descriptions.

The official supplemental archive was obtained through the EBI BioStudies copy
of the PMC supplement. Its SHA-256 is
`86841c4a4b6f69561770305ce490578fb614221e2d29de2fb5b3b21315a57b41`.

No checked source publishes the clone-resolved repair junction. The Methods say
that clones were genotyped by PCR and sequencing, but the sequence result is not
included. Supplemental Table 4 reports boundary-deletion library read counts,
not genomic breakpoints.

## Independent hg38 mapping

Both published guides were searched exactly on both strands of the GRCh38
primary assembly. Each has one exact local match and a valid reverse-strand PAM.

| Guide (published 5' to 3') | hg38 reference protospacer | Strand | Protospacer, 1-based inclusive | Guide-orientation PAM | Expected cut, 0-based |
|---|---|---:|---:|---:|---:|
| `TACTATCTATAGTAACTCCC` | `GGGAGTTACTATAGATAGTA` | - | chr3:4,976,065-4,976,084 | GGG | 4,976,067 |
| `TACCAGACTTCCACCGTATC` | `GATACGGTGGAAGTCTGGTA` | - | chr3:4,976,788-4,976,807 | AGG | 4,976,790 |

The inferred clean deletion is therefore:

```text
0-based half-open:  chr3:[4,976,067, 4,976,790)
1-based inclusive: chr3:4,976,068-4,976,790
length:             723 bp
```

The previous coordinates were numerically correct. The terminology was made
more precise: **cut-to-cut**, not "sgRNA-bounded." The molecular allele may
contain indels or repair sequence at the junction and remains unavailable.

## Locus and contact-map cross-check

- Published Figure 3D locus: `chr3:4,803,502-5,144,387` (340,886 bp,
  1-based inclusive).
- Published K562 BHLHE40 active sub-TAD from Supplemental Table 2:
  `chr3:4,977,047-5,056,047`.
- Historical contact source requested by the user:
  `/data/vinhtb/activhitracMIX_stage/data/user_k562_hitrac_mix/K562.hic`.
- Source SHA-256:
  `8ae9a61bf8f420d33f05f8c32747e520be57087650f1e6d0b499c2893673f9b3`.
- Extraction: Juicer `dump observed NONE`, `BP 1000`, enclosing grid
  `chr3:4,803,000-5,145,000`.
- Matrix: 342 x 342; 55,196 upper-triangle contacts; 8,714 nonzero
  upper-triangle cells.

This historical `.hic` is not the pooled and depth-matched 33-million-cis-PET
control displayed in the published Figure 3D. The two lineages must remain
separate.

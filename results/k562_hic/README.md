# Historical K562.hic locus extraction

This directory contains the exact 1-kb `observed NONE` contact records used in
the Gviz/contact-map composite.

## Source

```text
/data/vinhtb/activhitracMIX_stage/data/user_k562_hitrac_mix/K562.hic
SHA-256: 8ae9a61bf8f420d33f05f8c32747e520be57087650f1e6d0b499c2893673f9b3
```

## Extraction

The bin-aligned interval enclosing the published Figure 3D locus was dumped
directly with Juicer Tools 2.20.00:

```bash
java -jar juicer_tools.2.20.00.jar dump observed NONE \
  K562.hic \
  3:4803000:5145000 \
  3:4803000:5145000 \
  BP 1000
```

No KR, VC or VC_SQRT balancing and no observed/expected transformation were
applied. The stored table includes only the three numeric Juicer fields with a
header added by the local renderer.

## Verified result

- matrix geometry: 342 x 342;
- upper-triangle contact mass: 55,196;
- diagonal mass: 14,969;
- nonzero upper-triangle cells: 8,714;
- display only: `log10(observed count + 1)`, fixed white-red range 0-2.

This historical `.hic` is a different data lineage from the pooled and
33-million-cis-PET-downsampled Figure 3D experimental control.

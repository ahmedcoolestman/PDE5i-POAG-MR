# Genetically proxied PDE5 inhibition and primary open-angle glaucoma

Code for the drug-target Mendelian randomization study:

> El Ashry A, Higgins G, Naraine L, Elhusseiny AM. *Genetically Proxied PDE5 Inhibition and Risk of Primary Open-Angle Glaucoma: A Drug-Target Mendelian Randomization Study.* (under review)

## What's here

| Script | What it does | Paper |
|---|---|---|
| `R/01_main_analysis.R` | Main two-sample MR: SBP-lowering allele of rs2389873 at *PDE5A* against six outcomes (IGGC POAG, FinnGen POAG / primary open-angle / all glaucoma / NTG, IOP). Palindrome assessment, strand sensitivity, Steiger, BH-FDR and Bonferroni. | Table 2, Figure 2, Supp. Tables S2–S3 |
| `R/02_pqtl_sensitivity.R` | PDE5A plasma pQTL (Suhre 2017) as an alternative exposure, plus the PheWAS of rs2389873. | Supp. Tables S4–S5 |

Run them in that order. Estimates are per 1 mmHg lower SBP (script 01) and per 1 SD lower PDE5A protein (script 02).

## Data

All data are publicly available GWAS summary statistics. No individual-level data were used.

| Trait | Source | ID |
|---|---|---|
| Systolic blood pressure (exposure) | Evangelou et al. 2018, *Nat Genet* | OpenGWAS `ieu-b-38` |
| POAG (primary outcome) | IGGC, Gharahkhani et al. 2021, *Nat Commun* (UK Biobank excluded) | GWAS Catalog `GCST90011767` |
| POAG, primary open-angle glaucoma, all glaucoma, NTG | FinnGen release 5 | OpenGWAS `finn-b-H7_GLAUCOMA_POAG`, `finn-b-H7_GLAUCPRIMOPEN`, `finn-b-H7_GLAUCOMA`, `finn-b-H7_GLAUCOMA_NTG` |
| Intraocular pressure | Springelkamp et al. 2017, *Hum Mol Genet* | OpenGWAS `ebi-a-GCST004074` |
| PDE5A plasma protein (sensitivity) | Suhre et al. 2017, *Nat Commun* (KORA F4) | OpenGWAS `prot-c-5256_86_3` |

The IGGC file isn't on OpenGWAS, so download it first from the GWAS Catalog:

```
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90011001-GCST90012000/GCST90011767/
```

Save `GCST90011767_buildGRCh37.tsv.gz` in the folder you run the scripts from. Everything else is pulled from OpenGWAS through the API.

## Setup

1. R ≥ 4.3 (analysis run on R 4.6.1).
2. Install packages:
   ```r
   install.packages(c("remotes", "data.table", "dplyr", "ggplot2"))
   remotes::install_github("MRCIEU/TwoSampleMR")   # v0.7.9 used
   remotes::install_github("mrcieu/ieugwasr")      # v1.1.0.9000 used
   ```
3. Get a free OpenGWAS token at <https://api.opengwas.io> and add it to `~/.Renviron`:
   ```
   OPENGWAS_JWT=your_token_here
   ```
   Restart R after saving. Tokens expire after about two weeks.

## Running

```
Rscript R/01_main_analysis.R     # -> results_v4/
Rscript R/02_pqtl_sensitivity.R  # -> results_pqtl/
```

`sessionInfo.txt` in `results_v4/` records the exact package versions.

## Licence

MIT. See `LICENSE`.

## Contact

Ahmed El Ashry, University of Bristol

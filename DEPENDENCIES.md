# Dependencies

The exact versions used are recorded in the four manifest files below. They were generated
from the environment in which the analyses ran, not reconstructed afterwards.

| File | What it covers |
| :--- | :--- |
| `requirements.txt` | base Python 3.8.10 environment — most of the analysis |
| `requirements-scmalig.txt` | isolated Python environment — malignancy probe and Waddington-OT |
| `requirements-scib.txt` | isolated Python environment — integration benchmark |
| `R_packages.txt` | R 4.2.2 packages, with the library each was loaded from |

## Why there is more than one Python environment

Two of the tools used here cannot coexist with the rest of the stack. Installing
`scMalignantFinder` or Waddington-OT into the shared environment downgrades `torch` and
breaks `scvi-tools`, so each was given its own environment. This is a real constraint of
this analysis rather than a preference, and it is recorded because a reader trying to
reproduce the malignancy probe or the optimal-transport arm will need the second
environment, while everything else runs in the first.

The three environments used were:

```
base                 Python 3.8.10   everything except the two arms below
~/venvs/scmalig      Python 3.8.10   scMalignantFinder 1.2.0, Waddington-OT 1.0.8.post2, scvi-tools 0.15.5
~/venvs/scib         Python 3.8.10   scib 1.1.5
```

## Why the R libraries are split

R packages load from several libraries, and four of them are deliberately isolated:

```
/usr/local/lib/R/site-library          shared; most packages
~/R/x86_64-pc-linux-gnu-library/4.2    user library
~/Rlibs/fastcnv                        fastCNV, Banksy, Seurat 5.0.1, leidenAlg, sctransform
~/Rlibs/infercnv                       infercnv 1.23.0
~/Rlibs/sigsearch                      signatureSearch 1.12.0
~/Rlibs/nichenet                       nichenetr 2.2.1.1
```

Scripts must be run with the full library path on `R_LIBS`, because a plain
`Rscript --vanilla` drops `/usr/local/lib/R/site-library` and the scripts will fail on
`data.table`, which lives there. The manifest records the library for each package so the
required path can be reconstructed.

## Non-package dependencies

These are not installable from a package index and are needed by specific arms:

- **JAGS 4.3.2**, built locally (used by `infercnv`). Not the system JAGS.
- **SCMG**, the gene-perturbation library and its `CausalGenePredictor`. The source is
  cloned rather than installed, and the scripts add it to `sys.path` directly. The
  accompanying data library is `xingjiepan/SCMG_data` on HuggingFace (MIT).
- **PLIP model weights**, downloaded into `models/plip/`. The weights themselves are not in
  this repository; the configuration and tokenizer files are, so that the model version is
  pinned.
- **NicheNet reference matrices** (ligand–receptor network and ligand–target matrix), from
  the project's Zenodo release.
- **A local LINCS database.** The compound-perturbation arm queries a
  `signatureSearch`-format HDF5 file built from LINCS L1000 Level 5 (GSE70138). The file is
  5.9 GB and is not in this repository, but it is fully regenerable: `10_niche/70_write_db_fast.py`
  rebuilds it in about 26 seconds from the downloaded matrix.

## Reproducing exactly

The manifests pin versions, not hashes. Two further records in the repository pin the parts
of the analysis where a version difference would change a result:

- `exports/luad_v2_export_20260924/MANIFEST.sha256` — checksums for frozen intermediate
  artifacts.
- `docs/PARAMETERS_AND_SOURCES.md` — every analysis parameter, labelled by whether it comes
  from a source paper, a software default, or from us.

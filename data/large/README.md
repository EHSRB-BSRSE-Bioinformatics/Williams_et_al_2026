# Large Monte Carlo objects

The publication figures can be rebuilt without these files by using the
compact summaries in `data/derived/`. The files below are needed only to repeat
the expensive summary-building stage.

The Zenodo deposit is pending. After downloading it, extract or symlink the
objects into this layout:

```text
data/large/
├── DEGs_simulations/
│   ├── spheroids.RDS
│   ├── tk6.RDS
│   ├── heparg.RDS
│   └── mcf7.RDS
├── GSEA_simulations/
│   ├── spheroids_gsea.RDS
│   ├── tk6_gsea.RDS
│   ├── heparg_gsea.RDS
│   └── mcf7_gsea.RDS
├── ORA_simulations/
│   ├── spheroids_ora.RDS
│   ├── tk6_ora.RDS
│   ├── heparg_ora.RDS
│   └── mcf7_ora.RDS
└── TGx-DDI_simulations/
    └── *.RDS
```

Then run:

```bash
pixi run rebuild-summaries
```

Set `FDR_LARGE_ROOT` to use the same directory layout elsewhere. Set
`FDR_GSEA_DIR` to override only the GSEA directory. Large objects under this
directory are intentionally ignored by Git.

When the Zenodo upload is complete, add the DOI, direct archive URL, archive
checksums, and version/date here.


The local `provenance/` subdirectory retains the original study DEG script for
internal provenance. It is ignored by Git. Public readers use the cleaned
implementation and source checksum described in `docs/deseq2_methods.md`.
Selected archived DEG outputs are bundled separately in `data/reference/`;
`deg_example/` and `deg_validation/` hold disposable example and check outputs.

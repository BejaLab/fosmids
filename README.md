BejaLab fosmid assembly and annotation pipeline
==============================================

This repository contains a Snakemake workflow for fosmid processing from reads or pre-assembled contigs.

Input data
----------

The workflow reads `config/samples.xlsx` with columns:

- `id`: sample identifier (used in output paths).
- `illumina1`, `illumina2`: paired-end read files relative to `raw/`.
- `nanopore`: long-read file relative to `raw/`.
- `pre_assembled`: assembly FASTA relative to `pre_assembled/`.
- `prefix`: required; used as `locus_tag_prefix` and prepended to LOCUS in final GenBank.
- `ref_organism`: organism for PGAP `submol.yaml` (`genus_species`); defaults to `Escherichia coli` when empty.
- `organism`: organism string written into the final transformed `.gbk`.

Assembly source selection for each sample is:
1. `illumina1`/`illumina2` -> SPAdes
2. otherwise `nanopore` -> Flye
3. otherwise `pre_assembled`

Workflow overview
-----------------

Depending on available inputs, the workflow performs:

- read trimming for Illumina data
- vector filtering on reads using a middle region of the vector
- assembly (`spades` or `flye`)
- vector trimming from assembled contigs
- PGAP annotation
- taxonomy classification with Metabuli (`metabuli classify`)

Default final targets:

- `output/{id}/{id}.gbk`
- `output/{id}/report.xlsx` with sheets: `Assembly`, `Annotation`, `Classification`, `Genes`, `Vector`

How to run
----------

1. Fill `config/samples.xlsx`.
2. Place input files in `raw/` and/or `pre_assembled/` according to metadata.
3. Run Snakemake (profile enables Conda):

```bash
snakemake
```

```bash
snakemake --use-conda --cores 10
```

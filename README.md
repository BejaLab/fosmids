BejaLab's fosmid assembly and annotation pipeline
=================================================

The pipeline takes the list of fosmids in `fosmids.txt` that contains three columns: fosmid name and left and right reads.

The pipeline will run:

* quality trimming with `trim_galore`
* assembly with `spades`
* vector trimming with `cutadapt` (the vector flanks are specified in `workflow/metadata/vector_{fw,rv}.fna`)
* annotation with `pgap` (using default metadata for _Escherichia coli_)
* conversion to `.tbl` with gbf2tbl.pl (from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/toolbox/ncbi_tools/converters/scripts/gbf2tbl.pl))

The main output files are:

* `analysis/fosmids/{id}.fna` -- trimmed fosmid sequence
* `analysis/fosmids/{id}.info` -- vector trimming report
* `analysis/pgap/{id}/annot.tbl` -- annotations in `.tbl` format
* `analysis/pgap/{id}/annot.gbf` -- annotations in `.tbl` format

How to run:

* add fosmids to `fosmids.txt`
* check that the raw data are there
* run: `snakemake -c{threads} --use-conda`

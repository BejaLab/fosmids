import gzip

MIN_LEN = snakemake.params.min_len

def write_fastq_record(handle, name: str, seq: str, qual: str) -> None:
    if len(seq) != len(qual):
        raise ValueError(f"Sequence/quality length mismatch for read '{name}'")
    handle.write(f"@{name}\n{seq}\n+\n{qual}\n")

input_info ,= snakemake.input
output_trimmed = snakemake.output.trimmed
output_passed = snakemake.output.passed

with open(input_info) as inp, gzip.open(output_trimmed, "wt") as trimmed, gzip.open(output_passed, "wt") as passed:
    for raw in inp:
        line = raw.rstrip("\n")
        if not line:
            continue

        cols = line.split("\t")
        name = cols[0].replace(' ', '_')

        # No vector hit: keep the read as-is.
        if len(cols) == 4:
            seq, qual = cols[2], cols[3]
            write_fastq_record(passed, name, seq, qual)
            continue

        # Vector hit: emit one surviving chunk as-is, or two as _L/_R.
        left_seq, left_qual = cols[4], cols[8]
        right_seq, right_qual = cols[6], cols[10]

        has_left = len(left_seq) >= MIN_LEN
        has_right = len(right_seq) >= MIN_LEN

        if has_left and has_right:
            write_fastq_record(trimmed, f"{name}_L", left_seq, left_qual)
            write_fastq_record(trimmed, f"{name}_R", right_seq, right_qual)
        elif has_left:
            write_fastq_record(trimmed, name, left_seq, left_qual)
        elif has_right:
            write_fastq_record(trimmed, name, right_seq, right_qual)

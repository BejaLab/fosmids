import re
import warnings
from pathlib import Path

import pandas as pd
from openpyxl.styles import Font

warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=r"The global interpreter lock \(GIL\) has been enabled to load module 'Bio\.Align\._aligncore'.*",
)

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


GENE_LIKE_TYPES = {"gene", "CDS", "tRNA", "rRNA", "ncRNA", "tmRNA"}
GENE_CLASSIFICATION_COLUMNS = ["score", "lineage", "taxID:match_count"]


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def annotation_value_text(value) -> str:
    if isinstance(value, list):
        return "; ".join(str(v) for v in value)
    if isinstance(value, dict):
        return "; ".join(f"{k}={v}" for k, v in value.items())
    return str(value)


def format_coordinates(location) -> str:
    return str(location).replace("[", "").replace("]", "")


def with_prefix(name: str, prefix: str) -> str:
    value = str(name).strip()
    pref = f"{prefix}_"
    return value if value.startswith(pref) else f"{pref}{value}"


def update_locus_line(record: SeqRecord, prefix: str) -> str:
    old_locus = str(record.name or record.id).strip()
    if not old_locus:
        raise ValueError("Missing locus identifier in GenBank record")
    return with_prefix(old_locus, prefix)


def transform_genbank(records: list[SeqRecord], prefix: str) -> list[SeqRecord]:
    transformed: list[SeqRecord] = []
    for record in records:
        new_locus = update_locus_line(record, prefix)
        record.id = new_locus
        record.name = new_locus
        record.description = new_locus

        record.annotations["references"] = []
        transformed.append(record)

    return transformed


def mean_depth(depth_path: str | None) -> float | None:
    if not depth_path:
        return None
    path = Path(depth_path)
    if not path.exists():
        return None

    total = 0
    count = 0
    with open(path) as inp:
        for line in inp:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            total += int(cols[2])
            count += 1

    return (total / count) if count else 0.0


def vector_iteration_stats(sample_id: str, read_mode: str, flye_inputs: list[str]) -> dict[str, str]:
    default = "NA"
    if read_mode != "nanopore" or not flye_inputs:
        return default
    print(flye_inputs)
    return str(len(flye_inputs) - 1)


def optional_input_path(value) -> str | None:
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        return str(value[0])
    text = str(value)
    return text if text else None


def normalize_seq(seq: str) -> str:
    return str(seq).strip().upper()


def is_end_match(observed: str, expected_end: str, max_error_rate: float = 0.05) -> bool:
    if not observed or not expected_end:
        return False
    if len(observed) > len(expected_end):
        return False

    max_errors = int(len(observed) * max_error_rate)
    mismatches = sum(a != b for a, b in zip(observed, expected_end))
    return mismatches <= max_errors


def parse_trim_info(path: str) -> dict[str, dict[str, object]]:
    parsed: dict[str, dict[str, object]] = {}
    with open(path) as inp:
        for raw in inp:
            line = raw.strip()
            if not line:
                continue

            cols = line.split("\t")
            if len(cols) < 7:
                continue

            contig = cols[0]
            if contig.endswith(" rc"):
                contig = contig[:-3]
            matched_seq = normalize_seq(cols[5])
            if not matched_seq:
                continue
            try:
                position = int(cols[2])
            except ValueError:
                continue

            row = parsed.setdefault(contig, {"positions": [], "matched_sequences": []})
            row["positions"].append(position)
            row["matched_sequences"].append(matched_seq)
    return parsed


def build_assembly_sheet(
    trim_info_path: str, pretrim_fasta_path: str, trimmed_fasta_path: str, vector_fasta_path: str
) -> pd.DataFrame:
    trim_info = parse_trim_info(trim_info_path)
    vector_record = next(SeqIO.parse(vector_fasta_path, "fasta"), None)
    if vector_record is None:
        raise ValueError(f"No vector sequence found in FASTA: {vector_fasta_path}")
    vector_seq = normalize_seq(vector_record.seq)
    before_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(pretrim_fasta_path, "fasta")}
    after_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(trimmed_fasta_path, "fasta")}

    contigs = set(before_lengths) | set(after_lengths) | set(trim_info)

    def sort_length(contig: str) -> int:
        after_len = after_lengths.get(contig)
        if after_len is not None:
            return after_len
        before_len = before_lengths.get(contig)
        if before_len is not None:
            return before_len
        return 0

    rows = []
    for contig in sorted(contigs, key=lambda c: (-sort_length(c), c)):
        info = trim_info.get(contig, {"positions": [], "matched_sequences": []})
        matched_sequences = info["matched_sequences"]
        left_end_observed = any(
            is_end_match(seq, vector_seq[-len(seq) :]) for seq in matched_sequences
        )
        right_end_observed = any(
            is_end_match(seq, vector_seq[: len(seq)]) for seq in matched_sequences
        )
        before_len = before_lengths.get(contig)
        if before_len is None:
            before_len = after_lengths.get(contig)

        rows.append(
            {
                "Contig": contig,
                "Length before vector trimming": before_len,
                "Length after vector trimming": after_lengths.get(contig),
                "Vector found": yes_no(bool(matched_sequences)),
                "Vector on both sides of insert": yes_no(
                    left_end_observed and right_end_observed
                ),
                "Vector 3' end (left side) observed": yes_no(left_end_observed),
                "Vector 5' end (right side) observed": yes_no(right_end_observed),
            }
        )

    return pd.DataFrame(rows)


def build_annotation_sheet(records: list, prefix: str) -> pd.DataFrame:
    rows = []
    for record in records:
        structured_comment = record.annotations.get("structured_comment", {})
        genome_annotation_data = {}
        if isinstance(structured_comment, dict):
            genome_annotation_data = structured_comment.get("Genome-Annotation-Data", {})

        rows.append({"Contig": with_prefix(record.name, prefix), "Key": "", "Value": ""})
        rows.append(
            {
                "Contig": "",
                "Key": "Organism",
                "Value": annotation_value_text(record.annotations.get("organism", "")),
            }
        )
        if isinstance(genome_annotation_data, dict):
            for key, value in genome_annotation_data.items():
                if str(key) == "Annotation Provider":
                    continue
                rows.append(
                    {
                        "Contig": "",
                        "Key": str(key),
                        "Value": annotation_value_text(value),
                    }
                )
        elif genome_annotation_data:
            rows.append(
                {
                    "Contig": "",
                    "Key": "Genome-Annotation-Data",
                    "Value": annotation_value_text(genome_annotation_data),
                }
            )
        rows.append({"Contig": "", "Key": "", "Value": ""})
    return pd.DataFrame(rows)


def normalize_id(value: str) -> set[str]:
    raw = str(value).strip()
    if not raw:
        return set()
    raw = raw.lstrip(">")
    variants = {raw}
    variants.add(raw.split()[0])
    if "|" in raw:
        parts = [p for p in raw.split("|") if p]
        if parts:
            variants.add(parts[-1])
    cds_match = re.search(r"_cds_(.+)_\d+$", raw)
    if cds_match:
        variants.add(cds_match.group(1))
    tail_match = re.search(r"_cds_(.+)_\d+$", raw.split("|")[-1])
    if tail_match:
        variants.add(tail_match.group(1))
    return {v for v in variants if v}


def read_relaxed_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        engine="python",
        usecols=lambda c: not str(c).startswith("Unnamed"),
    )


def load_gene_classification_map(path: str) -> dict[str, dict[str, object]]:
    df = read_relaxed_tsv(path)
    if df.empty:
        return {}

    key_col = "name" if "name" in df.columns else df.columns[0]
    class_map: dict[str, dict[str, object]] = {}
    for _, row in df.iterrows():
        payload = {c: row[c] if c in df.columns else pd.NA for c in GENE_CLASSIFICATION_COLUMNS}
        for key in normalize_id(row[key_col]):
            if key not in class_map:
                class_map[key] = payload
    return class_map


def build_genes_sheet(records: list, class_map: dict[str, dict[str, object]]) -> pd.DataFrame:
    merged: dict[str, dict[str, object]] = {}

    for record in records:
        for feature in record.features:
            if feature.type not in GENE_LIKE_TYPES:
                continue

            locus_tags = feature.qualifiers.get("locus_tag", [])
            if not locus_tags:
                continue
            locus_tag = locus_tags[0]

            row = merged.setdefault(
                locus_tag,
                {
                    "locus_tag": locus_tag,
                    "coordinates": format_coordinates(feature.location),
                    "product": "",
                    "protein_sequence": "",
                    "nucleotide_sequence": "",
                    "score": pd.NA,
                    "lineage": pd.NA,
                    "taxID:match_count": pd.NA,
                },
            )

            product = feature.qualifiers.get("product", [""])[0]
            if product and not row["product"]:
                row["product"] = product

            if not row["nucleotide_sequence"]:
                row["nucleotide_sequence"] = str(feature.extract(record.seq))

            if feature.type == "CDS":
                translation = feature.qualifiers.get("translation", [""])[0]
                if translation:
                    row["protein_sequence"] = translation

            class_data = class_map.get(locus_tag)
            if class_data:
                for col in GENE_CLASSIFICATION_COLUMNS:
                    row[col] = class_data.get(col, pd.NA)

    rows = [merged[k] for k in sorted(merged)]
    return pd.DataFrame(rows)


def build_vector_sheet(
    insert_depth_path: str | None,
    vector_depth_path: str | None,
    read_mode: str,
    sample_id: str,
    flye_inputs: list[str],
) -> pd.DataFrame:
    insert_mean = mean_depth(insert_depth_path)
    vector_mean = mean_depth(vector_depth_path)
    ratio = "NA"
    if insert_mean is not None and vector_mean is not None and insert_mean != 0:
        ratio = f"{(vector_mean / insert_mean):.6f}"

    vector_iter = vector_iteration_stats(sample_id, read_mode, flye_inputs)
    return pd.DataFrame(
        [
            {
                "Read mode": read_mode,
                "Insert mean depth": "NA" if insert_mean is None else f"{insert_mean:.6f}",
                "Vector mean depth": "NA" if vector_mean is None else f"{vector_mean:.6f}",
                "Vector/insert depth ratio": ratio,
                "Maximum number of vector copies": vector_iter
            }
        ]
    )


def make_header_rows_bold(writer: pd.ExcelWriter, sheet_names: list[str]) -> None:
    bold_font = Font(bold=True)
    for sheet_name in sheet_names:
        ws = writer.book[sheet_name]
        for cell in ws[1]:
            cell.font = bold_font


input_dir = snakemake.input.pgap
trim_info = snakemake.input.trim
trimmed_fasta = snakemake.input.fasta
pretrim_fasta = snakemake.input.pretrim_fasta
vector_fasta = snakemake.input.vector
metabuli_classifications = snakemake.input.metabuli
metabuli_report = snakemake.input.metabuli_report
insert_depth = optional_input_path(getattr(snakemake.input, "insert_depth", None))
vector_depth = optional_input_path(getattr(snakemake.input, "vector_depth", None))
flye_inputs = [str(p) for p in getattr(snakemake.input, "flye_input", [])]
sample_id = str(snakemake.wildcards.id)

output_gbk = snakemake.output.gbk
output_report = snakemake.output.report
prefix = snakemake.params.prefix
read_mode = str(getattr(snakemake.params, "read_mode", "none"))

in_gbk = Path(input_dir) / "annot.gbk"
original_records = list(SeqIO.parse(in_gbk, "genbank"))
transformed_records = transform_genbank(list(SeqIO.parse(in_gbk, "genbank")), prefix=prefix)
SeqIO.write(transformed_records, output_gbk, "genbank")

assembly_df = build_assembly_sheet(trim_info, pretrim_fasta, trimmed_fasta, vector_fasta)
annotation_df = build_annotation_sheet(original_records, prefix)
classification_df = read_relaxed_tsv(metabuli_report)
gene_class_map = load_gene_classification_map(metabuli_classifications)
genes_df = build_genes_sheet(transformed_records, gene_class_map)
vector_df = build_vector_sheet(insert_depth, vector_depth, read_mode, sample_id, flye_inputs)

with pd.ExcelWriter(output_report, engine="openpyxl") as writer:
    assembly_df.to_excel(writer, index=False, sheet_name="Assembly")
    vector_df.to_excel(writer, index=False, sheet_name="Vector")
    annotation_df.to_excel(writer, index=False, sheet_name="Annotation")
    classification_df.to_excel(writer, index=False, sheet_name="Classification")
    genes_df.to_excel(writer, index=False, sheet_name="Genes")
    make_header_rows_bold(
        writer,
        ["Assembly", "Vector", "Annotation", "Classification", "Genes"],
    )

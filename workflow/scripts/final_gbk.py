import re
from pathlib import Path


def find_input_genbank(input_dir: str) -> Path:
    base = Path(input_dir)
    preferred = ("annot.gbf", "annot.gbk", "annot.gb")
    for name in preferred:
        candidate = base / name
        if candidate.exists():
            return candidate

    matches = []
    for pattern in ("*.gbf", "*.gbk", "*.gb"):
        matches.extend(sorted(base.glob(pattern)))

    if not matches:
        raise FileNotFoundError(f"No GenBank file found in '{input_dir}'")
    return matches[0]


def update_locus_line(line: str, prefix: str) -> str:
    match = re.match(r"^(LOCUS\s+)(\S+)(.*)$", line.rstrip("\n"))
    if not match:
        raise ValueError(f"Unexpected LOCUS line format: {line.rstrip()}")
    head, old_locus, tail = match.groups()
    new_locus = f"{prefix}_{old_locus}"
    return f"{head}{new_locus}{tail}\n", new_locus


def transform_genbank(lines: list[str], organism: str, prefix: str) -> list[str]:
    out: list[str] = []
    current_locus = None
    i = 0

    while i < len(lines):
        line = lines[i]

        if line.startswith("LOCUS"):
            updated_line, current_locus = update_locus_line(line, prefix)
            out.append(updated_line)
            i += 1
            continue

        if line.startswith("DEFINITION"):
            if current_locus is None:
                raise ValueError("DEFINITION encountered before LOCUS")
            out.append(f"DEFINITION  {current_locus}\n")
            i += 1
            while i < len(lines) and lines[i].startswith(" " * 12):
                i += 1
            continue

        if line.startswith("SOURCE"):
            source_prefix = re.match(r"^(SOURCE\s+)", line).group(1)
            out.append(f"{source_prefix}{organism}\n")
            i += 1
            continue

        if line.startswith("  ORGANISM"):
            out.append(f"  ORGANISM  {organism}\n")
            i += 1
            while i < len(lines) and lines[i].startswith(" " * 12):
                i += 1
            continue

        if line.startswith("REFERENCE"):
            i += 1
            while i < len(lines):
                nxt = lines[i]
                if nxt and not nxt.startswith(" "):
                    break
                i += 1
            continue

        if re.match(r'^\s*/db_?xref="taxon:[^"]*"\s*$', line):
            i += 1
            continue

        org_match = re.match(r'^(\s*/organism=)".*"\s*$', line)
        if org_match:
            out.append(f'{org_match.group(1)}"{organism}"\n')
            i += 1
            continue

        src_match = re.match(r'^(\s*/source=)".*"\s*$', line)
        if src_match:
            out.append(f'{src_match.group(1)}"{organism}"\n')
            i += 1
            continue

        out.append(line)
        i += 1

    return out


input_dir, = snakemake.input
output_gbk, = snakemake.output
organism = snakemake.params.organism
prefix = snakemake.params.prefix

in_gbk = find_input_genbank(input_dir)
with open(in_gbk) as inp:
    transformed = transform_genbank(inp.readlines(), organism=organism, prefix=prefix)

Path(output_gbk).parent.mkdir(parents=True, exist_ok=True)
with open(output_gbk, "w") as out:
    out.writelines(transformed)

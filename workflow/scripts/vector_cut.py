import urllib.request
import io
import warnings

warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=r"The global interpreter lock \(GIL\) has been enabled to load module 'Bio\.Align\._aligncore'.*",
)

from Bio import SeqIO
from Bio import BiopythonParserWarning
import Bio.Restriction as Restriction

warnings.filterwarnings(
    "ignore",
    category=BiopythonParserWarning,
    message=r"Attempting to parse malformed locus line:.*",
)
# Fetch parameters from Snakemake
url = snakemake.params.url
enzyme_name = snakemake.wildcards.enzyme
output_fasta ,= snakemake.output

# 1. Download the GenBank file directly into memory
try:
    req = urllib.request.Request(
        url, 
        headers={'User-Agent': 'Mozilla/5.0 (Snakemake Pipeline)'}
    )
    with urllib.request.urlopen(req) as response:
        gbk_data = response.read().decode('utf-8')
except Exception as e:
    raise RuntimeError(f"Failed to download GenBank file from URL: {url}\nError: {e}")

# 2. Parse the text data into a BioPython SeqRecord
record = SeqIO.read(io.StringIO(gbk_data), "genbank")

# 3. Look up the enzyme in BioPython's restriction database
enzyme = getattr(Restriction, enzyme_name, None)
if enzyme is None:
    raise ValueError(f"Enzyme '{enzyme_name}' not recognized by BioPython.")

# 4. Search for the restriction site (ensuring circular topology)
cut_sites = enzyme.search(record.seq, linear=False)

if len(cut_sites) == 0:
    raise ValueError(f"Enzyme {enzyme_name} site not found in the plasmid sequence.")
elif len(cut_sites) > 1:
    raise ValueError(f"Enzyme {enzyme_name} cuts {len(cut_sites)} times. It must cut exactly once.")

top_cut = cut_sites[0]
ovhg = enzyme.ovhg

# 5. Shift/linearize the sequence around the cut site
linear_seq = record.seq[top_cut - 1:] + record.seq[:top_cut - 1]

# 6. Dynamically trim single-stranded overhangs to yield blunt ends
if ovhg > 0:
    # 5' overhang: Trim bases from the start
    linear_seq = linear_seq[ovhg:]
elif ovhg < 0:
    # 3' overhang: Trim bases from the end
    linear_seq = linear_seq[:ovhg]

# 7. Update header and write out to the FASTA file
record.seq = linear_seq

SeqIO.write(record, output_fasta, "fasta")

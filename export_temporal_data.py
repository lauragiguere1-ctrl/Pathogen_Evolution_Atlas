import json
from collections import Counter

# --- Load mutation calls and collection dates ---
with open("mutation_calls.json") as f:
    genome_mutations = json.load(f)

dates = {}
with open("data/raw/genomes.ndjson") as f:
    for line in f:
        record = json.loads(line)
        dates[record["accession"]] = record["collection_date"]

def to_quarter(date_str):
    year = date_str[:4]
    month = int(date_str[5:7])
    quarter = (month -1) // 3 + 1
    return f"{year}-Q{quarter}"

# --- The three trustworthy quarters ---
GOOD_QUARTERS = ["2025-Q3", "2025-Q4", "2026-Q1"]
genomes_by_quarter = {q: [] for q in GOOD_QUARTERS}
for accession in genome_mutations:
    q = to_quarter(dates[accession])
    if q in genomes_by_quarter:
        genomes_by_quarter[q].append(accession)

quarter_n = [len(genomes_by_quarter[q]) for q in GOOD_QUARTERS]

# --- Lineage family per genome (from Nextclade) ---
lineage = {}
with open("nextclade_output/nextclade.tsv") as f:
    header = f.readline().split("\t")
    name_col = header.index("seqName")
    pango_col = header.index("Nextclade_pango")
    for line in f:
        cols = line.split("\t")
        family = cols[pango_col].strip().split(".")[0]
        lineage[cols[name_col].strip()] = family

# --- Lineage share by quarter ---
lineage_share = {}
for fam in ["XFG", "NB", "PQ"]:
    shares = []
    for q in GOOD_QUARTERS:
        genomes_this_q = genomes_by_quarter[q]
        carriers = sum(1 for acc in genomes_this_q if lineage.get(acc) == fam)
        shares.append(round(100 * carriers / len(genomes_this_q), 1))
    lineage_share[fam] = shares

# --- Mutation frequency by quarter ---
def frequency_in_quarter(mutation, accessions):
    carriers = sum(1 for acc in accessions if list(mutation) in genome_mutations[acc])
    return round(100 * carriers / len(accessions), 1)

watch_list = [
    [22995, "C", "A"], [23277, "C", "T"], [23021, "A", "G"],
    [20178, "C", "T"], [13608, "C", "T"], [10615, "C", "T"], [7113, "C", "T"],
]

mutation_freq = {}
for mut in watch_list:
    pos, ref, alt = mut
    name = f"{ref}{pos}{alt}"
    mutation_freq[name] = [frequency_in_quarter(mut, genomes_by_quarter[q]) for q in GOOD_QUARTERS]

# --- Assemble and write ---
export = {
    "quarters": GOOD_QUARTERS,
    "quarter_n": quarter_n,
    "lineage_share": lineage_share,
    "mutation_freq": mutation_freq,
}

with open("data/site/temporal.json", "w") as out:
    json.dump(export, out, indent=2)

print("Wrote data/site/temporal.json")
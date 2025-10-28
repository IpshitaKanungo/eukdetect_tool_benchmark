"""
This script extracts per-tool FASTA files for contigs where:
  ground_truth == False  AND  tool == False
(i.e., "true negatives" for each tool).

Input:
  --tsv   : Path to ground truth TSV (columns: contig_id, ground_truth, eukrep, tiara, kraken2)
  --fasta : Path to unfiltered assembly FASTA (can be .gz)

Output:
  --out-dir/{eukrep_false.fasta, tiara_false.fasta, kraken2_false.fasta}
"""
#!/usr/bin/env python3
import argparse, sys, gzip, io, os, csv
from collections import defaultdict

# To open file, also handles .gz
def open_in(path, mode='rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode=mode)
    return open(path, mode=mode, encoding=None if 'b' in mode else 'utf-8')

"""
Extract per-tool FASTA files for contigs where ground truth == false and tool column == false
Parsing command line arguments
"""
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True, help="Path to ground truth TSV")
    ap.add_argument("--fasta", required=True, help="Path to assembly FASTA")
    ap.add_argument("--id-col", default="contig_id", help="Column with contig IDs (default: contig_id).")
    ap.add_argument("--gt-col", default="ground_truth", help="Ground truth column (default: ground_truth")
    ap.add_argument("--eukrep-col", default="eukrep", help="EukRep column (default: eukrep)")
    ap.add_argument("--tiara-col", default="tiara", help="Tiara column (default: tiara)")
    ap.add_argument("--kraken2-col", default="kraken2", help="kraken column (default: kraken2)")
    ap.add_argument("--out-dir", default=".", help="Output directory (default: current dir)")
    ap.add_argument("--id-whitespace-split", action="store_true", help="Match FASTA IDs by splitting header at first whitespace")
    return ap.parse_args()

# Convert string to True/False or None
def normalize_bool_str(s):
    if s is None:
        return None
    s = s.strip().lower()
    if s in ("true","t","1","yes","y"):
        return True
    if s in ("false","f","0","no","n"):
        return False
    # return None for ambiguous/other
    return None

"""
load tsv and collect contigs where ground truth == false and tool column == false
returns a dictionary of set of contig ids/tool
"""
def load_id_sets(tsv_path, id_col, gt_col, eukrep_col, tiara_col, kraken2_col):
    tsv_ids = {
        "eukrep": set(),
        "tiara": set(),
        "kraken2": set()
    }
    n_rows = 0
    n_selected = defaultdict(int)
    with open(tsv_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in (id_col, gt_col, eukrep_col, tiara_col, kraken2_col) if c not in reader.fieldnames]
        if missing:
            sys.exit(f"ERROR: Missing columns in TSV: {missing}. Found columns: {reader.fieldnames}")
        for row in reader:
            n_rows += 1
            cid = row[id_col]
            gt = normalize_bool_str(row.get(gt_col))
            euk = normalize_bool_str(row.get(eukrep_col))
            tia = normalize_bool_str(row.get(tiara_col))
            kra = normalize_bool_str(row.get(kraken2_col))
            # Keep contig ids where ground truth is False AND tool is False
            if gt is False and euk is False:
                tsv_ids["eukrep"].add(cid)
            if gt is False and tia is False:
                tsv_ids["tiara"].add(cid)
            if gt is False and kra is False:
                tsv_ids["kraken2"].add(cid)
    return tsv_ids, n_rows

# Retrieves header and sequence from fasta file
def fasta_iter(handle, split_header=True):
    header = None
    seq_chunks = []
    for line in handle:
        if not line:
            continue
        if line[0] == '>':
            if header is not None:
                yield header, ''.join(seq_chunks)
            header = line[1:].rstrip('\n\r')
            if split_header:
                header = header.split()[0]
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())
    if header is not None:
        yield header, ''.join(seq_chunks)

# Scans the FASTA, write sequences into per tool FASTA, if the ID is selected
def write_selected_fasta(fasta_path, out_dir, id_sets, split_header=True):
    os.makedirs(out_dir, exist_ok=True)
    out_paths = {
        "eukrep": os.path.join(out_dir, "eukrep_false.fasta"),
        "tiara": os.path.join(out_dir, "tiara_false.fasta"),
        "kraken2": os.path.join(out_dir, "kraken2_false.fasta"),
    }
    outs = {k: open(v, "w", encoding="utf-8") for k, v in out_paths.items()}
    counts_found = {k: 0 for k in id_sets}
    seen = set()
    missing = {k: set() for k in id_sets}
    # Pre-fill missing with all; we'll remove those we find
    for k in id_sets:
        missing[k] = set(id_sets[k])
    # Stream FASTA and route sequences to corresponding tools FASTA
    with open_in(fasta_path, "rt") as fh:
        for hid, seq in fasta_iter(fh, split_header=split_header):
            for tool, idset in id_sets.items():
                if hid in idset:
                    outs[tool].write(f">{hid}\n")
                    # Wrap sequence at 80 chars/line for readability
                    for i in range(0, len(seq), 80):
                        outs[tool].write(seq[i:i+80] + "\n")
                    counts_found[tool] += 1
                    if hid in missing[tool]:
                        missing[tool].remove(hid)

    for f in outs.values():
        f.close()

    return out_paths, counts_found, missing

def main():
    args = parse_args()
    id_sets, n_rows = load_id_sets(args.tsv, args.id_col, args.gt_col, args.eukrep_col, args.tiara_col, args.kraken2_col)
    print(f"Loaded {n_rows} rows from {args.tsv}", file=sys.stderr)
    for tool in ("eukrep","tiara","kraken2"):
        print(f"{tool}: {len(id_sets[tool])} contig IDs selected (ground_truth=false AND {tool}=false)", file=sys.stderr)
    out_paths, counts_found, missing = write_selected_fasta(
        args.fasta, args.out_dir, id_sets, split_header=args.id_whitespace_split
    )
    for tool in ("eukrep","tiara","kraken2"):
        print(f"{tool}: wrote {counts_found[tool]} sequences to {out_paths[tool]}", file=sys.stderr)
        if missing[tool]:
            print(f"{tool}: WARNING {len(missing[tool])} IDs not found in FASTA (first 10 shown): {sorted(list(missing[tool]))[:10]}", file=sys.stderr)

if __name__ == "__main__":
    main()

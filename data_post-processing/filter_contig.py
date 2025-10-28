"""
This script merges two TSVs on the column 'contig_id':
  1) A stats table that includes at least: contig_id, length, aligned_bases, percent_aligned
  2) A labels table that includes per-tool binary calls for: eukrep, tiara, kraken2

It then adds a new column 'ground_truth' derived from percent_aligned:
  - percent_aligned >= 90  -> "True"   (eukaryotic)
  - percent_aligned <= 10  -> "False"  (non-eukaryotic)
  - percent_aligned 11<p<89  -> "ambiguous"

Inputs
- TSV with per-contig alignment statistics (must contain 'contig_id' and 'percent_aligned';
  typically also 'length' and 'aligned_bases').
- TSV with per-contig tool labels (must contain 'contig_id' and columns for tools such as
  'eukrep', 'tiara', and 'kraken2', with boolean-like values e.g. True/False).

Output
- A merged TSV with columns:
    contig_id, length, aligned_bases, percent_aligned, eukrep, tiara, kraken2, ground_truth
"""

#!/usr/bin/env python3
import argparse, gzip, sys, math
from pathlib import Path
from collections import OrderedDict

def open_in(p):
    # Open files and also handles .gz
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "rt")

def norm(s):
    # Normalize header keys for matching (case/space insensitive)
    return s.strip().lower().replace(" ", "_")

def pick(hmap, *cands):
    # Return the original header name for the first matching candidate
    for c in cands:
        k = norm(c)
        if k in hmap:
            return hmap[k]
    return None

def read_labels(path):
    """
    Read labels TSV -> (labels_by_contig, tool_cols) generated from euk_non-euk.py
      Expects a 'contig_id' column.
      Ignores a 'length' column if present in labels.
      All remaining columns are treated as tool columns (e.g., eukrep, tiara, kraken2).
    """
    labels, tool_cols = {}, []
    with open_in(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if not header or len(header) < 2:
            sys.exit("ERROR: labels TSV must have at least contig_id and one tool column.")

        hmap = {norm(h): h for h in header}
        contig_col = pick(hmap, "contig_id", "contig", "id", "name")
        if contig_col is None:
            sys.exit("ERROR: labels TSV must have a contig_id/contig column.")

        length_col = pick(hmap, "length")  # Not considered in this case
        tool_cols = [h for h in header if h not in (contig_col, length_col)]

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            row = dict(zip(header, parts))
            cid = row.get(contig_col, "").split()[0]
            if not cid:
                continue
            labels[cid] = {tool: row.get(tool, "-") for tool in tool_cols}

    return labels, tool_cols

def read_stats(path):
    """
    Read percent_aligned.tsv generated from alignment_percent.py
    Returns (rows list preserving order, header list, contig_col, length_col)
    """
    rows = []
    with open_in(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if not header or len(header) < 2:
            sys.exit("ERROR: stats TSV must have at least contig_id and length columns.")

        hmap = {norm(h): h for h in header}
        contig_col = pick(hmap, "contig_id", "contig", "id", "name")
        length_col = pick(hmap, "length", "len", "contig_length", "size")
        if contig_col is None or length_col is None:
            sys.exit(f"ERROR: cannot find contig_id and length columns in stats header: {header}")

        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            row = dict(zip(header, parts))
            rows.append(row)

    return rows, header, contig_col, length_col

def parse_float_maybe_percent(x):
    """parse string to float; return NaN if not parseable."""
    if x is None:
        return float("nan")
    s = str(x).strip()
    if s.endswith("%"):
        s = s[:-1].strip()
    try:
        return float(s)
    except ValueError:
        return float("nan")

"""Parsing command-line argument"""
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-l", "--labels", required=True) # Labels TSV (from eukrep/tiara/kraken2)
    ap.add_argument("-s", "--stats",  required=True) # Stats TSV (has contig_id, length, percent_aligned, etc.)
    ap.add_argument("-o", "--out",    required=True) # Output TSV
    ap.add_argument("--min-len", type=int, default=1500) # Minimum contig length to keep (default: 1500)
    args = ap.parse_args()

    labels, tool_cols = read_labels(args.labels)
    stats_rows, stats_header, stats_cid, stats_len = read_stats(args.stats)

    # Keep only desired stats columns (exclude contig_id/length and drops)
    other_stats_cols = [h for h in stats_header if h not in (stats_cid, stats_len)]

    # Build header and de-duplicate by normalized name
    out_header = ["contig_id", "length"] + other_stats_cols
    seen = {norm(h) for h in out_header}
    tool_cols_nd = [c for c in tool_cols if norm(c) not in seen]
    out_header += tool_cols_nd + ["ground_truth"]

    # Find the percent_aligned column to compute ground_truth from stats
    hmap_stats = {norm(h): h for h in stats_header}
    percent_col = pick(hmap_stats, "percent_aligned", "percent", "pct_aligned", "percent_euk_alignment")
    if percent_col is None:
        sys.exit("ERROR: cannot find a percent_aligned column in stats to compute ground_truth.")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    kept = drop_len = 0

    with open(args.out, "w") as out:
        out.write("\t".join(out_header) + "\n")

        for row in stats_rows:
            cid = row.get(stats_cid, "").split()[0]
            if not cid:
                continue

            # Parse length from stats
            Lraw = row.get(stats_len, "0")
            try:
                L = int(float(Lraw))
            except ValueError:
                continue
                # Apply minimun length filter
                if L < args.min_len:
                drop_len += 1
                continue

            # Compute ground_truth from percent_aligned
            p = parse_float_maybe_percent(row.get(percent_col, ""))
            if math.isnan(p):
                gt = "ambiguous"
            elif p <= 10.0:
                gt = "False"      # non-euk
            elif p >= 90.0:
                gt = "True"       # euk
            else:
                gt = "ambiguous"  # 11<p<89

            # LEFT JOIN behavior: keep stats row; fill missing tool labels with '-'
            lab = labels.get(cid, {})
            drop_targets = {"kraken2", "eukrep", "tiara"}  # Normalized names
            should_drop = any(
                norm(tool) in drop_targets and lab.get(tool, "-").strip() == "-"
                for tool in tool_cols_nd
            )
            if should_drop:
                continue
            out_row = OrderedDict()
            out_row["contig_id"] = cid
            out_row["length"]    = str(L)
            for col in other_stats_cols:
                out_row[col] = row.get(col, "")
            for tool in tool_cols_nd:
                out_row[tool] = lab.get(tool, "-")
            out_row["ground_truth"] = gt

            out.write("\t".join(out_row.get(h, "") for h in out_header) + "\n")
            kept += 1

    print(f"{kept} rows to {args.out} (dropped <{args.min_len} bp: {drop_len})", file=sys.stderr)

if __name__ == "__main__":
    main()



#!/usr/bin/env python3
import argparse, csv, gzip, os, sys
from pathlib import Path

# ---------- IO helpers ----------
def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode, encoding=None if "b" in mode else "utf-8")

def fasta_len(path):
    """Return total number of bp across all sequences in a FASTA file."""
    total = 0
    with open_maybe_gz(path, "rt") as fh:
        for line in fh:
            if not line or line.startswith(">"):
                continue
            total += len(line.strip())
    return total

def find_bins_dir(ds_dir, user_bins_dir=None):
    """Pick a bins directory. Tries user override, then common names."""
    candidates = []
    if user_bins_dir:
        candidates.append(Path(user_bins_dir))
    for name in ("metawrap_50_10_bins", "metawrap_50_10_bins.stats", "bins", "bin_fastas"):
        candidates.append(Path(ds_dir) / name)
    for c in candidates:
        if c.exists() and c.is_dir():
            return c
    return None

def find_bin_fasta(bins_dir, bin_name):
    """Find the FASTA file for a bin by trying common extensions."""
    exts = (".fa", ".fasta", ".fa.gz", ".fasta.gz", ".fna", ".fna.gz")
    for ext in exts:
        p = Path(bins_dir) / f"{bin_name}{ext}"
        if p.exists():
            return p
    # also try with any prefix/suffix pattern (fallback, slower)
    for p in Path(bins_dir).glob(f"*{bin_name}*"):
        if p.suffix in {".fa", ".fasta", ".fna"} or str(p).endswith((".fa.gz", ".fasta.gz", ".fna.gz")):
            return p
    return None

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(
        description="Collect dataset, bin, completeness, contamination, total_bp from MetaWRAP outputs."
    )
    ap.add_argument(
        "--dataset",
        action="append",
        metavar="NAME=DIR",
        required=True,
        help="Dataset name and directory, e.g. metagenome=/path/to/metagenome (repeat for each dataset).",
    )
    ap.add_argument(
        "--stats-name",
        default="metawrap_50_10_bins.stats",
        help="Filename of the stats table inside each dataset dir (default: metawrap_50_10_bins.stats).",
    )
    ap.add_argument(
        "--bins-dir-name",
        default=None,
        help="If bins are in a specific subdir name (same for all datasets), provide it; otherwise the script will guess.",
    )
    ap.add_argument(
        "-o", "--out",
        default="bins_merged.tsv",
        help="Output TSV path (default: bins_merged.tsv).",
    )
    args = ap.parse_args()

    # parse dataset mappings
    datasets = {}
    for item in args.dataset:
        if "=" not in item:
            sys.exit(f"ERROR: --dataset expects NAME=DIR, got: {item}")
        name, ddir = item.split("=", 1)
        ddir = ddir.strip()
        if not name or not ddir:
            sys.exit(f"ERROR: malformed dataset spec: {item}")
        if not Path(ddir).exists():
            sys.exit(f"ERROR: dataset dir not found: {ddir}")
        datasets[name] = ddir

    rows_out = []
    missing_fastas = []

    for name, ddir in datasets.items():
        ddir = Path(ddir)
        stats_path = ddir / args.stats_name
        if not stats_path.exists():
            sys.exit(f"[{name}] ERROR: stats file not found: {stats_path}")

        bins_dir = find_bins_dir(ddir, args.bins_dir_name)
        if bins_dir is None:
            sys.exit(f"[{name}] ERROR: could not locate bins directory under {ddir}. "
                     f"Pass --bins-dir-name if needed.")

        # read stats table
        with open(stats_path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            required_cols = {"bin", "completeness", "contamination"}
            if not required_cols.issubset(reader.fieldnames or []):
                sys.exit(f"[{name}] ERROR: stats file missing columns. "
                         f"Need {required_cols}, found {reader.fieldnames}")
            for row in reader:
                bin_name = str(row["bin"]).strip()
                comp = row["completeness"]
                cont = row["contamination"]

                fasta_path = find_bin_fasta(bins_dir, bin_name)
                if fasta_path is None:
                    missing_fastas.append((name, bin_name, str(bins_dir)))
                    total_bp = 0
                else:
                    total_bp = fasta_len(fasta_path)

                rows_out.append({
                    "dataset": name,
                    "bin": bin_name,
                    "completeness": comp,
                    "contamination": cont,
                    "total_bp": total_bp
                })

    # write output
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with open(outp, "w", encoding="utf-8", newline="") as outfh:
        w = csv.DictWriter(outfh, fieldnames=["dataset","bin","completeness","contamination","total_bp"], delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)

    print(f"Wrote {len(rows_out)} rows to {outp}")
    if missing_fastas:
        print("\nWARNING: FASTA not found for the following bins (total {}):".format(len(missing_fastas)))
        for ds, b, bd in missing_fastas[:20]:
            print(f"  [{ds}] bin={b} (searched in {bd})")
        if len(missing_fastas) > 20:
            print("  ... (more omitted)")

if __name__ == "__main__":
    main()


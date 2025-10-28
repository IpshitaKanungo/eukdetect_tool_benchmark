"""
This script splits assembly FASTA into EUK and NONEUK contigs based on Tiara predictions
- EUK classes:     {"eukarya", "mitochondrion", "plastid"}
- NON-EUK classes: {"bacteria", "archaea", "prokarya", "unknown"}
Inputs:
  --assembly : assembly FASTA (plain or .gz)
  --tiara    : Tiara output table (whitespace-separated; has sequence_id and class columns)
Outputs:
  --euk-out  : output FASTA path for EUK contigs (plain or .gz)
  --noneuk-out : output FASTA path for NON-EUK contigs (plain or .gz)
Options
  --prefer-second : prefer Tiara's second-stage class if present (and not N/A/"-")
  --accept-fasta  : only keep contigs whose IDs appear in this FASTA
"""

#!/usr/bin/env python3
import argparse, gzip
from pathlib import Path
from Bio import SeqIO

# EUK & NONEUK based on the classes
EUK  = {"eukarya", "mitochondrion", "plastid"}
NONEUK = {"bacteria", "archaea", "prokarya", "unknown"}

# Open file for reading
def open_in(p):  return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "rt")
# Open file for writing
def open_out(p): return gzip.open(p, "wt") if str(p).endswith(".gz") else open(p, "wt")


"""
Build dict: contig_token -> chosen_class (lowercased).
- Expects whitespace-separated lines with at least 3 fields.
- Uses FIRST token as the contig token (e.g., 'k141_198329').
- Chooses class from the last two columns:
    * If --prefer-second: use the final column when it's not N/A/'-'/''; else use the previous one.
    * Otherwise: use the previous (first-stage) column.
"""
def load_tiara(tbl, prefer_seccond
    m = {}
    with open(tbl, "rt") as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"): continue
            low = s.lower()
            if low.startswith("sequence_id") and "class" in low:
                continue
            parts = s.split()
            if len(parts) < 3: continue
            token = parts[0]
            fst, snd = parts[-2].lower(), parts[-1].lower() # first-stage, second-stage
            chosen = snd if (prefer_second and snd not in {"n/a","na","-",""}) else fst
            m[token] = chosen
    return m


# Return set of IDs present in a FASTA
def load_id_whitelist(fa):
    ids = set()
    with open_in(fa) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            ids.add(rec.id)
    return ids

# Split assembly into EUK/NONEUK using Tiara table
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-a", "--assembly",   required=True, help="Assembly FASTA (.fa/.fasta/.gz).")
    ap.add_argument("-t", "--tiara",      required=True, help="Tiara table (whitespace-separated).")
    ap.add_argument("-e", "--euk-out",    required=True, help="Output EUK FASTA (.fa/.gz).")
    ap.add_argument("-n", "--noneuk-out", required=True, help="Output NON-EUK FASTA (.fa/.gz).")
    ap.add_argument("--prefer-second", action="store_true", help="Prefer second-stage class if present.")
    ap.add_argument("--accept-fasta",
                    help="Only keep contigs whose IDs appear in this FASTA; others are skipped.")
    args = ap.parse_args()

    # ensure output folders exist
    Path(args.euk_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.noneuk_out).parent.mkdir(parents=True, exist_ok=True)

    # Load Tiara class map and optional whitelist IDs
    tiara = load_tiara(args.tiara, args.prefer_second)
    accept_ids = load_id_whitelist(args.accept_fasta) if args.accept_fasta else None

    # Counters for the summary
    euk_n = noneuk_n = 0
    fasta_only = 0          # in assembly but no Tiara class
    not_accepted = 0        # in assembly but not in accept_ids (if provided)
    other_class_skipped = 0 # class not in EUK/NONEUK

    # Scans the assembly and directs contigs to intended output
    with open_in(args.assembly) as inh, open_out(args.euk_out) as eukh, open_out(args.noneuk_out) as noneukh:
        for rec in SeqIO.parse(inh, "fasta"):
            token = rec.id
            if accept_ids is not None and token not in accept_ids:
                # Mark as N/A/'-' conceptually (we don't write it to either FASTA)
                not_accepted += 1
                continue
            cls = tiara.get(token)
            if not cls:
                fasta_only += 1
                continue
            if cls in EUK:
                SeqIO.write(rec, eukh, "fasta"); euk_n += 1
            elif cls in NONEUK:  # includes 'unknown'
                SeqIO.write(rec, noneukh, "fasta"); noneuk_n += 1
            else:
                other_class_skipped += 1

    # Summary
    print("prefer_second_stage\t", str(args.prefer_second).lower(), sep="")
    print("euk_contigs_written\t", euk_n, sep="")
    print("noneuk_contigs_written\t", noneuk_n, sep="")
    print("fasta_contigs_without_tiara_class\t", fasta_only, sep="")
    if accept_ids is not None:
        print("assembly_contigs_not_in_accept_fasta\t", not_accepted, sep="")
    print("other_class_skipped\t", other_class_skipped, sep="")

if __name__ == "__main__":
    main()


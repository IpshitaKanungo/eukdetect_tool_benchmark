"""
Required input: assembly fasta, Quast tsv: all_alignments-contigs.tsv 
output : tsv file with columns :contig_id, length, aligned_bases and percent_aligned
"""
#!/usr/bin/env python3
import argparse, gzip, csv, sys
from pathlib import Path
from collections import defaultdict

#open files and also handles .gz
def open_in(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "rt")

#Calculate sequence length from FASTA
def compute_lengths_from_fasta(fasta_path):
    lengths = {}
    with open_in(fasta_path) as fh:
        cid, clen = None, 0
        for line in fh:
            if line.startswith(">"):
                if cid is not None:
                    lengths[cid] = clen
                cid  = line[1:].strip().split()[0]
                clen = 0
            else:
                clen += len(line.strip())
        if cid is not None:
            lengths[cid] = clen
    return lengths

#merge overlapping contig locations and returns contig length
def merge_len(intervals):
    if not intervals:
        return 0
    iv = sorted((min(a,b), max(a,b)) for a,b in intervals)
    total = 0
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce + 1:
            if e > ce:
                ce = e
        else:
            total += (ce - cs + 1)
            cs, ce = s, e
    total += (ce - cs + 1)
    return total

def norm_header(xs):
    return [x.lstrip("\ufeff").strip().lower() for x in xs]

#Search index of header column
def find_idx_any(header, names):
    low = norm_header(header)
    for n in names:
        nlow = n.lower()
        if nlow in low:
            return low.index(nlow)
    return -1

#Parsing QUAST all alignment tsv file anf find  contig intervals (Start S/ End E)
def load_quast_intervals(path):
    intervals = defaultdict(list)
    with open_in(path) as f:
        rdr = csv.reader(f, delimiter="\t")
        try:
            header = next(rdr)
        except StopIteration:
            return intervals
        ci = find_idx_any(header, ["Contig", "contig", "ctg", "name", "contig id", "ctg_name"])
        if ci < 0:
            raise KeyError(f"Cannot find contig column in header: {header}")
        si = find_idx_any(header, ["S2","s2"])
        ei = find_idx_any(header, ["E2","e2"])
        if si < 0 or ei < 0:
            si = find_idx_any(header, ["ctg_start","contig start","contig_start","start on ctg","start"])
            ei = find_idx_any(header, ["ctg_end","contig end","contig_end","end on ctg","end"])
            if si < 0 or ei < 0:
                raise KeyError(f"Cannot find contig start/end columns in header: {header}")
        maxi = max(ci, si, ei)
        for row in rdr:
            if not row or len(row) <= maxi:
                continue
            cfield = row[ci].strip()
            if not cfield or cfield.startswith("#"):
                continue
            contig = cfield.split()[0]
            sraw, eraw = row[si].strip(), row[ei].strip()
            if not sraw or not eraw:
                continue
            try:
                s = int(sraw); e = int(eraw)
            except Exception:
                continue
            intervals[contig].append((s, e))
    return intervals

#Parsing arguments
def main():
    ap = argparse.ArgumentParser()
    #assembled fasta file
    ap.add_argument("-i","--fasta", required=True)
    #all alignement tsv file from quast
    ap.add_argument("-q","--quast-alignments", required=True)
    #output file
    ap.add_argument("-o","--out-tsv", required=True)
    args = ap.parse_args()

    """compute contig length from FASTA"""
    lengths = compute_lengths_from_fasta(args.fasta)

    #error message if no contigs found
    if not lengths:
        print(f"No contig lengths loaded.", file=sys.stderr)

    #Parse QUAST alignment and find merged coverage per contig
    intervals = load_quast_intervals(args.quast_alignments)
    aligned_bp = {c: merge_len(iv) for c, iv in intervals.items()}

    #To check existing directory
    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w") as out:
        out.write("contig_id\tlength\taligned_bases\tpercent_aligned\n") #column names
        for contig, L in lengths.items():
            cov = aligned_bp.get(contig, 0)
            if cov < 0:
                cov = 0
            if cov > L:
                cov = L
            pct_contig = (100.0 * cov / L) if L > 0 else 0.0
            out.write(f"{contig}\t{L}\t{cov}\t{pct_contig:.2f}\n")

if __name__ == "__main__":
    main()


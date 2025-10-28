"""



"""
import argparse, re, gzip
from pathlib import Path
from Bio import SeqIO

# Taxonomic roots
# EUK_ROOTS: Eukaryota ; NONEUK_ROOTS:Unclassified, Bacteria, Archaea, Viruses
EUK_ROOTS  = {"2759"}
NONEUK_ROOTS = {"0", "2", "2157", "10239"}

# Regex to extract Taxonomic ids
TX_RE = re.compile(r'(\d+):\d+')

# Open file for reading text
def oin(p):  return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "rt")
# Open file for writing text
def oout(p): return gzip.open(p, "wt") if str(p).endswith(".gz") else open(p, "wt")

"""
Checks Kraken2 taxonomic report and collects taxids under EUK and NONEUK roots
"""
def tax_sets_from_report(report, euk_roots=EUK_ROOTS, noneuk_roots=NONEUK_ROOTS):
    euk, noneuk, at = set(), set(), {}
    with oin(report) as f:
        for ln in f:
            if not ln.strip(): continue
            parts = ln.rstrip("\n").split(None, 5)
            if len(parts) < 6: continue
            taxid, ind = parts[4], parts[5]

            # Counting indentation level; 2 spaces = one level
            lvl = (len(ind) - len(ind.lstrip(" "))) // 2

            # Receives parent flag
            p_e = at.get(lvl-1, (None, False, False))[1] if lvl>0 else False
            p_n = at.get(lvl-1, (None, False, False))[2] if lvl>0 else False

            # Node is in EUK/NONEUK if the parent was or if the node is root
            in_e = p_e or (taxid in euk_roots)
            in_n = p_n or (taxid in noneuk_roots)
            if in_e: euk.add(taxid)
            if in_n: noneuk.add(taxid)
            at[lvl] = (taxid, in_e, in_n)
    return euk, noneuk


#Parse Kraken2 classification O/P and decide if contig is EUK or NON-EUK
def contig_sets_from_krak2(krak2, euk_tax, noneuk_tax, noneuk_priority=False):
    euk_ids, noneuk_ids = set(), set()
    with oin(krak2) as f:
        for s in f:
            if not s.strip(): continue
            cols = s.rstrip("\n").split("\t")
            if len(cols) < 2: continue
            cid = cols[1]
            taxids = set(TX_RE.findall(s))
            if len(cols) >= 3 and cols[2].isdigit(): taxids.add(cols[2])
            hit_e, hit_n = bool(taxids & euk_tax), bool(taxids & noneuk_tax)
            if hit_e and hit_n:
                (noneuk_ids if noneuk_priority else euk_ids).add(cid)
            elif hit_e: euk_ids.add(cid)
            elif hit_n: noneuk_ids.add(cid)
    # To check if contig appears in both
    if noneuk_priority: euk_ids -= noneuk_ids
    else:             noneuk_ids -= euk_ids
    return euk_ids, noneuk_ids

# Splits FASTA in EUK and NONEUK output
def split_fasta(assembly, euk_ids, noneuk_ids, euk_out, noneuk_out):
    Path(euk_out).parent.mkdir(parents=True, exist_ok=True)
    Path(noneuk_out).parent.mkdir(parents=True, exist_ok=True)
    euk_n = noneuk_n = skip = 0
    with oin(assembly) as inp, oout(euk_out) as eo, oout(noneuk_out) as no:
        for rec in SeqIO.parse(inp, "fasta"):
            rid = rec.id
            if rid in euk_ids:  SeqIO.write(rec, eo, "fasta");  euk_n += 1
            elif rid in noneuk_ids: SeqIO.write(rec, no, "fasta"); noneuk_n += 1
            else: skip += 1
    return euk_n, noneuk_n, skip

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tax-report", required=True, help="Kraken2 taxonomic report")
    ap.add_argument("--kraken2",    required=True, help="Kraken2 classification output")
    ap.add_argument("--assembly",   required=True, help="Assembly contigs FASTA")
    ap.add_argument("--euk-out",    required=True, help="Output euk FASTA")
    ap.add_argument("--noneuk-out",   required=True, help="Output noneuk FASTA")
    ap.add_argument("--noneuk-priority", action="store_true", help="If mixed hits, put contig in NONEUK.")
    args = ap.parse_args()

    euk_tax, noneuk_tax = tax_sets_from_report(args.tax_report)
    euk_ids, noneuk_ids = contig_sets_from_krak2(args.kraken2, euk_tax, noneuk_tax, args.noneuk_priority)
    e_n, n_n, s_n = split_fasta(args.assembly, euk_ids, noneuk_ids, args.euk_out, args.noneuk_out)

    print(f"EUK_taxids={len(euk_tax)} NONEUK_taxids={len(noneuk_tax)} | EUK_contigs={e_n} NONEUK_contigs={n_n} skipped={s_n}")

if __name__ == "__main__":
    main()

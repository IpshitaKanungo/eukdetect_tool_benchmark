"""
This script labels eukaryotic and non-eukaryotic contigs with a binary value for each tool
Required inputs : assembly fasta and eukaryotic and non eukaryotic contigs file for eukrep, kraken2 and tiara
output : tsv file with columns : contig_id, length, eukrep, kraken2, tiara
"""
#!/usr/bin/env python3
import argparse, gzip, sys
from pathlib import Path

#open files and also handles .gz
def open_in(p):
    return gzip.open(p, "rt") if p and str(p).endswith(".gz") else open(p, "rt")

def read_ids_and_lengths(fa_path):
    """
    Parsing FASTA ,returns contig id and length
    """
    ids, lens = set(), {}
    if not fa_path:
        return ids, lens
    with open_in(fa_path) as fh:
        cid, clen = None, 0
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                # removes previous contig when a new header appears
                if cid is not None:
                    ids.add(cid)
                    if cid in lens and lens[cid] != clen:
                        print(f"Length mismatch for {cid} in {fa_path}: "
                              f"{lens[cid]} vs {clen} (keeping max)", file=sys.stderr)
                        lens[cid] = max(lens[cid], clen)
                    else:
                        lens[cid] = clen
                cid  = line[1:].strip().split()[0] #new contig
                clen = 0
            else:
                clen += len(line.strip()) #calculate sequence length
        # remove the last contig
        if cid is not None:
            ids.add(cid)
            if cid in lens and lens[cid] != clen:
                print(f"Length mismatch for {cid} in {fa_path}: "
                      f"{lens[cid]} vs {clen} (keeping max)", file=sys.stderr)
                lens[cid] = max(lens[cid], clen)
            else:
                lens[cid] = clen
    return ids, lens

#retrieves only ids from FASTA
def read_ids_only(fa_path):
    ids, _ = read_ids_and_lengths(fa_path)
    return ids

#merge contig lengths, for multiple contig ids keep the maximum
def merge_lengths(into, new, label):
    for cid, L in new.items():
        if cid in into and into[cid] != L:
            print(f"Length mismatch for {cid} from {label}: "
                  f"{into[cid]} vs {L} (keeping max)", file=sys.stderr)
            into[cid] = max(into[cid], L)
        else:
            into[cid] = L
#assigns label for contig (euk/prok)
def val_for(contig, euk_set, prok_set):
    in_e = contig in euk_set
    in_p = contig in prok_set
    if in_e and in_p: return "conflict"
    if in_e:          return "True"
    if in_p:          return "False"
    return "-"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly") #Assembly FASTA; retrieves lengths and contig ids
    """
    eukaryotic or prokaryotic fasta for all the tools
    """
    # EukRep
    ap.add_argument("--eukrep-euk")
    ap.add_argument("--eukrep-prok")
    # Tiara
    ap.add_argument("--tiara-euk")
    ap.add_argument("--tiara-prok")
    # Kraken2
    ap.add_argument("--kraken2-euk")
    ap.add_argument("--kraken2-prok")

    ap.add_argument("-o", "--out-tsv", required=True)
    args = ap.parse_args()

    #determines which tools are present
    tools = []
    if args.eukrep_euk or args.eukrep_prok: tools.append("eukrep")
    if args.tiara_euk  or args.tiara_prok:  tools.append("tiara")
    if args.kraken2_euk or args.kraken2_prok: tools.append("kraken2")
    if not tools and not (args.assembly or args.lengths_tsv):
        sys.exit("ERROR: Provide at least one tool FASTA, or --assembly/--lengths-tsv to list contigs.")

    # collects euk/prok sets and lengths
    lengths = {}
    sets = {}
    #loading per tool euk/prok fasta and length merge
    if "eukrep" in tools:
        e_ids, e_lens = read_ids_and_lengths(args.eukrep_euk)
        p_ids, p_lens = read_ids_and_lengths(args.eukrep_prok)
        sets["eukrep"] = (e_ids, p_ids)
        merge_lengths(lengths, e_lens, "eukrep-euk")
        merge_lengths(lengths, p_lens, "eukrep-prok")

    if "tiara" in tools:
        e_ids, e_lens = read_ids_and_lengths(args.tiara_euk)
        p_ids, p_lens = read_ids_and_lengths(args.tiara_prok)
        sets["tiara"] = (e_ids, p_ids)
        merge_lengths(lengths, e_lens, "tiara-euk")
        merge_lengths(lengths, p_lens, "tiara-prok")

    if "kraken2" in tools:
        e_ids, e_lens = read_ids_and_lengths(args.kraken2_euk)
        p_ids, p_lens = read_ids_and_lengths(args.kraken2_prok)
        sets["kraken2"] = (e_ids, p_ids)
        merge_lengths(lengths, e_lens, "kraken2-euk")
        merge_lengths(lengths, p_lens, "kraken2-prok")

    #add or override lengths from assembly
    if args.assembly:
        _, ass_lens = read_ids_and_lengths(args.assembly)
        merge_lengths(lengths, ass_lens, "assembly")

    # Row IDs: union of all tool IDs; if assembly provided without tools, use assembly IDs
    all_ids = set()
    for pair in sets.values():
        eu, pr = pair
        all_ids |= eu
        all_ids |= pr
    if not all_ids and args.assembly:
        all_ids |= read_ids_only(args.assembly)

    # check output directory exists and write the tsv
    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w") as out:
        header = ["contig_id", "length"] + tools
        out.write("\t".join(header) + "\n")
        for cid in sorted(all_ids):
            row = [cid, str(lengths.get(cid, "-"))]
            for tool in tools:
                euk_set, prok_set = sets[tool]
                row.append(val_for(cid, euk_set, prok_set))
            out.write("\t".join(row) + "\n")

    print(f" Wrote {len(all_ids)} contigs to {args.out_tsv} "
          f"(columns: length, {', '.join(tools)})", file=sys.stderr)

if __name__ == "__main__":
    main()


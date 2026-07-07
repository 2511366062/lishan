#!/usr/bin/env python3
import gzip
import re
import sys
from pathlib import Path


gtf = Path(sys.argv[1])
out_tss = Path(sys.argv[2])
out_genes = Path(sys.argv[3])


def parse_attrs(text):
    attrs = {}
    for key, value in re.findall(r'([A-Za-z0-9_]+) "([^"]+)"', text):
        attrs[key] = value
    return attrs


opener = gzip.open if str(gtf).endswith(".gz") else open
seen = set()
with opener(gtf, "rt") as fin, out_tss.open("w") as tss, out_genes.open("w") as genes:
    for line in fin:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        chrom, _, _, start, end, _, strand, _, attrs_text = fields
        attrs = parse_attrs(attrs_text)
        gid = attrs.get("gene_id", "")
        gname = attrs.get("gene_name", gid)
        if not gid or (chrom, gid) in seen:
            continue
        seen.add((chrom, gid))
        start_i = int(start) - 1
        end_i = int(end)
        genes.write(f"{chrom}\t{start_i}\t{end_i}\t{gname}\t{gid}\t{strand}\n")
        pos = start_i if strand != "-" else end_i - 1
        tss.write(f"{chrom}\t{pos}\t{pos + 1}\t{gname}\t{gid}\t{strand}\n")

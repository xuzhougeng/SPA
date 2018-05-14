#!/usr/bin/env python3

import re
import sys
import fileinput
from os.path import join


if len(sys.argv) < 2:
    sys.exit()

bed = sys.argv[1]
gffpaths = []

params = [ tuple(line.strip().split("\t"))  for line in open(bed) ]
for param in params:
    dna     = param[0]
    d_start = param[1] if int(param[1]) - 2000 < 0 else int(param[1]) - 2000
    d_end   = param[2] if int(param[1]) - 2000 < 0 else int(param[2]) + 2000
    gffpath = join("prediction", "homology", "genewise", dna + "_" + str(d_start) + "_" + str(d_end), "genewise.gff")
    gffpaths.append(gffpath)

for line in fileinput.input(gffpaths):
    if line.startswith("//"):
        continue
    col = line.strip().split("\t")
    chrom, start, end = re.match(r"(.+?):(\d+)-(\d+)", col[0]).groups()
    pstart = int(start) + int(col[3]) - 1 # processed start
    pend   = int(start) + int(col[4]) - 1 # processed end
    if col[6] == "+":
        gff_out = [chrom, col[1],col[2], str(pstart), str(pend), col[5], col[6], col[7],col[8]]
    else:
       gff_out = [chrom, col[1],col[2], str(pend), str(pstart), col[5], col[6], col[7],col[8]]
    sys.stdout.write("\t".join(gff_out) + "\n")


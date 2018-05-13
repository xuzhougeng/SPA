#!/usr/bin/env python3

import re
import sys
import glob
import fileinput

files = glob.glob("prediction/homology/genewise/*/*.gff")

for line in fileinput.input(files):
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

    
 

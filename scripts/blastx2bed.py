#!/usr/bin/env python3
import sys
import os

if len(sys.argv) < 2:
    sys.exit(1)

blastx  = sys.argv[1]

dna     = ""
prt     = ""
d_start = 0
d_end   = 0
state   = False # record the pre-loop state 
first   = True  # omit the first line

file = open(blastx, 'r')
for line in file:
    items = line.strip().split("\t")
    unchange = True if dna == items[0] and prt == items[1] else False # first line will be false and omitted

    # if the pre-state is unchange and now change, which means next alignment
    if (state and not unchange) or (not state and not unchange) and not first:
        bed_list = [dna, str(d_start), str(d_end), prt, "0", ordered]
        print("\t".join(bed_list))
        # sequence can be substract with bedtools getfasta , seqkits subseq.
    if not unchange: # changed
        dna     = items[0]
        prt     = items[1]
        d_start = min( int(items[6]), int(items[7]) )
        d_end   = max( int(items[6]), int(items[7]) )
        ordered  = "+"  if int(items[6]) < int(items[7])  else "-"
    else: # unchanged, same as the line before
        d_start = min( d_start, int(items[6]), int(items[7]) )
        d_end   = max( d_end  , int(items[6]), int(items[7]) )

    state = unchange
    first = False

file.close()

#!/usr/bin/env python3

# version 0.1: it will run genewise parallel with Python multiprocessing library 
# version 0.2: run before testing the existed file

import os
import sys
import subprocess as sp
from multiprocessing import Pool
from os.path import exists 
from os.path import join

def sub_seq_from_bed(params):
    dna     = params[0]
    d_start = params[1] if int(params[1]) - 2000 < 0 else int(params[1]) - 2000
    d_end   = params[2] if int(params[1]) - 2000 < 0 else int(params[2]) + 2000
    prt     = params[3]
    strand  = params[5]
    genome  = params[6]
    protein = params[7]

    prefix = join("prediction", "homology", "genewise", dna + "_" + str(d_start) + "_" +str(d_end))
    if not exists(prefix):
        os.makedirs(prefix, exist_ok=True)

        # substract protein sequence
        prtpath  = join(prefix, "prt.fa")
        prtshell = "seqkit faidx {} {} -o {}".format(protein, prt, prtpath)
        sp.run(prtshell, shell=True)

    # substract dna sequence if reverse then reverse complement
    dnapath  = join(prefix, "dna.fa")
    if not exists(dnapath): 
        region   = str(d_start) + "-" + str(d_end)
        dnashell = "seqkit faidx {} {}:{} -o {}".format(genome, dna, region, dnapath )
        sp.run(dnashell, shell=True)
    
    # run genewise
    gffpath  = join(prefix, "genewise.gff")
    if not exists(gffpath): 
        logpath   = join(prefix, "genewise.err")
        if strand == "+":
            genewiseshell = "genewise {} {} -quiet -gff {} > {} 2> {}".format(prtpath, dnapath, "-tfor", gffpath, logpath)
        else:
            genewiseshell = "genewise {} {} -quiet -gff {} > {} 2> {}".format(prtpath, dnapath, "-trev", gffpath, logpath)
        #print(genewiseshell)
        sp.run(genewiseshell, shell=True)

if __name__ == "__main__":
    genome  = sys.argv[1]
    protein = sys.argv[2]
    bed     = sys.argv[3]
    params = [ tuple(line.strip().split("\t") + [genome, protein]) for line in open(bed) ]
    processes = 20
    with Pool(processes) as pool:
        pool.map(sub_seq_from_bed, params)

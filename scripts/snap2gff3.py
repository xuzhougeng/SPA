import os
import sys
import re

def multiple_exon_to_gff3(lines):
    # parse the input
    seqid  = lines.split("\t")[0] 
    source = lines.split("\t")[1]
    starts = []
    ends   = []
    for matches in re.findall(r"\t(\d+)\t(\d+)\t", lines):
        start, end = matches
        starts.append(int(start))
        ends.append(int(end))
    scores = "."
    strand = lines.split("\t")[6]
    phase  = "."
    ID     = lines.split("\t")[8].split("\n")[0]
    # contruct the output
    gene = [seqid, source, "gene", str(min(starts)), str(max(ends)), scores, strand, phase, "ID={}".format(ID)]
    mRNA = [seqid, source, "mRNA", str(min(starts)), str(max(ends)), scores, strand, phase, "ID={}-m1;Parent={}".format(ID,ID)]
    print('\t'.join(gene))
    print('\t'.join(mRNA))
    for num, (start, end) in enumerate(zip(starts, ends)):
        CDS  = [seqid, source, "CDS",  str(start), str(end), scores, strand, phase, "ID={}-cds{};Parent={}-m1".format(ID, num+1, ID)]
        exon = [seqid, source, "exon", str(start), str(end), scores, strand, phase, "ID={}-exon{};Parent={}-m1".format(ID, num+1,ID)]
        print('\t'.join(CDS))
        print('\t'.join(exon))
 

def single_exon_to_gff3(line):
    # parse the input
    seqid  = line.split("\t")[0]
    source = line.split("\t")[1]
    start  = line.split("\t")[3]
    end    = line.split("\t")[4]
    scores = "."
    strand = line.split("\t")[6]
    phase  = "."
    ID     = line.split("\t")[8].split("\n")[0]
    # contruct the output
    gene = [seqid, source, "gene", str(start), str(end), scores, strand, phase, "ID={}".format(ID)]
    mRNA = [seqid, source, "mRNA", str(start), str(end), scores, strand, phase, "ID={}-m1;Parent={}".format(ID,ID)]
    CDS  = [seqid, source, "CDS",  str(start), str(end), scores, strand, phase, "ID={}-cds{};Parent={}-m1".format(ID, 1, ID)]
    exon = [seqid, source, "exon", str(start), str(end), scores, strand, phase, "ID={}-exon{};Parent={}-m1".format(ID,1, ID)]
    print('\t'.join(gene))
    print('\t'.join(mRNA))
    print('\t'.join(CDS))
    print('\t'.join(exon))

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit()
    filepath = sys.argv[1]
    input = open(filepath)
    state = False
    for line in input:
        lines = []
        if re.search(r"\tEinit\t",  line):
            lines = []
            lines.append(line)
        elif re.search(r"\tExon\t", line):
            lines.append(line)
        elif re.search(r"\tEterm", line):
            lines.append(line)
            multiple_exon_to_gff3(''.join(lines))
            lines = []
        elif re.search(r"\tEsngl\t", line):
            single_exon_to_gff3(''.join(line))
        else:
            print(line)
    input.close()
    

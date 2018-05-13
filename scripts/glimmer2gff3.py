import io
import sys
import re

def glimmer_to_gff3(line):
    if re.search(r"\tmRNA\t", line): 
        gene = re.sub(r"\tmRNA\t",r"\tgene\t",line)
        mRNA = re.sub(r"\t(ID=(.*));", r"\t\1;Parent=\2;", line)
        mRNA = re.sub(r"gene(?=\d{1,3};P)", r"mrna",mRNA)
        sys.stdout.write( gene + mRNA )
    elif re.search(r"\tCDS\t", line):
        cds  = re.sub(r"gene(?=\d{1,3};N)", r"mrna", line)
        exon = re.sub(r"\tCDS\t",r"\texon\t", line)
        exon = re.sub(r"cds(?=\d)", "exon", exon)
        exon = re.sub(r"gene(?=\d{1,3};N)", r"mrna", exon)
        sys.stdout.write( cds + exon )
    else:
        print("warning: {}".format(line))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit()
    filepath = sys.argv[1]
    for line in io.open(filepath, 'r'):
        if not line.startswith("#"):
            glimmer_to_gff3(line)

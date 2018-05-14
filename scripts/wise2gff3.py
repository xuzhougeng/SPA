import io
import sys
import re

if len(sys.argv) <2:
    sys.exit()

filepath = sys.argv[1]
fp = open(filepath, 'r')

exon_id = 0
cds_id  = 0
choose  = 1
for line in fp:
    #column = line.split("\t")
    #col_id = column[0] + column[3] + column[4] 
    if line.startswith("#"): # or col_id in match_list:
        continue
    elif re.search(r"\tmatch\t", line):
        # construct gene feature
        gene_List    = line.strip().split("\t")
        gene_List[2] = "gene"
        gene_attr    = '_'.join([gene_List[0], gene_List[3], gene_List[4]]) + "_gene_" + str(choose)
        gene_List[8] = "ID={}".format(gene_attr)
        print("\t".join(gene_List))
        # construct mRNA feature
        mRNA_List    = line.strip().split("\t")
        mRNA_List[2] = "mRNA"
        mRNA_attr    = '_'.join([mRNA_List[0], mRNA_List[3], mRNA_List[4]]) + "_mRNA_" +  str(choose)
        mRNA_List[8] = "ID={};Parent={}".format(mRNA_attr, gene_attr)
        print("\t".join(mRNA_List))
        exon_id = 0
        mRNA_id = 0
        choose  = choose + 1
    elif re.search(r"\tcds\t", line):
        choose = choose + 1
        cds_id = cds_id + 1
        CDS_List = line.strip().split("\t")
        CDS_List[2] = "CDS"
        CDS_List[8] = "ID={}_cds{}.{};Parent={}".format(CDS_List[0], cds_id, choose, mRNA_attr)
        print("\t".join(CDS_List))
        exon_id = exon_id + 1
        exon_List = line.strip().split("\t")
        exon_List[2] = "exon"
        exon_List[8] = "ID={}_exon{}.{};Parent={}".format(exon_List[0], exon_id, choose, mRNA_attr)
        print("\t".join(exon_List))
    else:
        continue
fp.close()

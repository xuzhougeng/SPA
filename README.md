# SPA: A (s)nakemake based (p)lant genome (a)nnotation pipline

> 这个流程目前是非常的难用，估计会用的人就我一个，为了搞这个流程，我自己还写了好几个Python脚本来进行格式间的转换。

This pipeline is designed for plant genome annotation. 

If you want to use this tool, you may need to read [this](https://github.com/xuzhougeng/Notebook/blob/master/Notes/Pipeline/How-to-annotate-plant-genome.md). This documents or tutorial record the learning process when I explore genome annotaiton. 

Becauses this pipeline is writed by Snakemake, you need install this excellent program firstly

```bash
conda install snakemake
# or
pip install snakemake
```

You also need to install the following tools

For RNA-seq data process:

- [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)
- [Stringtie](https://ccb.jhu.edu/software/stringtie/)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki)

For gene prediction:

- [RepeatMasker](http://www.repeatmasker.org/) 
- [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/downloads/)
- [GeneMark](http://exon.gatech.edu/license_download.cgi)
- [GlimmerHMM](https://ccb.jhu.edu/software/glimmerhmm/)
- [BRAKER2](http://exon.gatech.edu/Braker/BRAKER2.tar.gz)
- [PASA](https://github.com/PASApipeline/PASApipeline/wiki)
- [EvidenceModeler](https://evidencemodeler.github.io/)

For intermediate file process:

- [NCBI-BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [seqkit](https://github.com/shenwei356/seqkit)
- [samtools](https://github.com/samtools/samtools)
- [bedtools](http://bedtools.readthedocs.io/en/latest/)

I haven't develop a tool to configure these tools easily, may docker or install.sh? Who knows.

Now, some tools can be installed by conda, and if not pipeline will be interrupution:

```bash
conda create -c bioconda -n align hisat2 samtools 
conda create -c bioconda -n assembly trinity
```

To run the pipline, you should follow the steps below:

- create a project directory and cp all the files under this directory to yor project direcotry
- create `input` directory under project's root directory, and cp or move the following files
    - genome.fa: genomic sequence
    - protein.fa: the homology species protein sequence
- create `RNA_seq` directory under the project's root directory, and cp or move your RNA-seq data. The sample name should be listed in the samples.txt. For example, the sample name of "AB\_1.fq, AB\_2.fq" is AB.
- some program path should be modified in the `config.yaml`, if you not configure these file, I promise you will come across with a lot of error

Before check-in, you need dry run the pipelne to check potential error:

```bash
snakemake -np
```

If no warnings errors interrupt operation, you could actually run the pipline

```bash
snakemake -j 30
# -j --jobs
```

## Scripts Usage

I write small tools for parallel running genewise and convert genewise's gff output to EvidenceModeler accepted format. Thses scripts are `blastx2bed.py`, `blastx_bed_merge.py`, `parallel_genewise.py`, `merge_parallel_genewise_out.py` and the `wise2gff3.py`.

> Make sure that your Python is 3.6+, and seqkit and genewise should be in your environments, The results will export to "prediction/homology/genewise"

```bash
blastx -query query_dna.fa -db target_protein.fa -outfmt 6 -best_hit_overhang 0.25 -best_hit_score_edge 0.25 > blastx_fmt6.out
python3 blastx2bed.py blastx_fmt6.out | bedtools sort -i - > blastx.bed.tmp
python3 blastx_bed_merge.py blastx.bed.tmp > blastx_merged.bed
python3 parallel_genewise.py query_dna.fa target_protein.fa blastx_merged.bed
python3 merge_parallel_genewise_out.py blastx_merged.bed > genewise.gff3
```

These scripts are just server for my pipeline, will not maintain, be care of these.

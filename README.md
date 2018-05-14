# SPA: A (s)nakemake based (p)lant genome (a)nnotation pipline

This pipeline is designed for plant genome annotation. 

If you want to use this tool, you may need to read <https://github.com/xuzhougeng/Notebook/blob/master/Notes/Pipeline/How-to-annotate-plant-genome.md>. This documents or tutorial record the learning process when I learning genome annotaiton. 

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

For gene prediction:

- [RepeatMasker](http://www.repeatmasker.org/) 
- [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/downloads/)
- [GeneMark](http://exon.gatech.edu/license_download.cgi)
- [SNAP](https://github.com/KorfLab/SNAP)
- [GlimmerHMM](https://ccb.jhu.edu/software/glimmerhmm/)
- [BRAKER2](http://exon.gatech.edu/Braker/BRAKER2.tar.gz)
- [PASA](https://github.com/PASApipeline/PASApipeline/wiki)
- [EvidenceModeler](https://evidencemodeler.github.io/)

For intermediate file process:

- [NCBI-BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [seqkit](https://github.com/shenwei356/seqkit)

I haven't develop a tool to configure these tools easily, may docker or install.sh? Who knows.

To run the pipline, you should follow the steps below:

- create a project directory and cp all the files under this directory to yor project direcotry
- create `input` directory under project's root directory, and cp or move the following files
    - genome.fa: genomic sequence
    - protein.fa: the homology species protein sequence
- create `RNA_seq` directory under the project's root directory, and cp or move your RNA-seq data. The sample name should be listed in the samples.txt. For example, the sample name of "AB_1.fq, AB_2.fq" is AB.
- some program path should be modified in the `config.yaml`, if you not configure these file, I promise you will come across with a lot of error

Before check-in, you need dry run the pipelne to check potential error:

```bash
snakemake -np
```

If no warnings errors interrupt operation, you could actually run the pipline

```bash
snakemake -j 30
```

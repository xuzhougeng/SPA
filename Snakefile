#!/usr/bin/env snakemake

# author: xuzhogueng
# contact mail: xuzhougeng@163.com

# version change
## version 0.1:   Finish the first version of this pipeline using braker2, SNAP, GlimmerHMM, genewise PASA and EVM
## version 0.1.1: Add cleanpasa and cleanevm for interruption of runnning PASA and EVM
## version 0.1.2: Remove SNAP as for it are not good as braker2 and GlimmerHMM when I check the results with naked eyes.
## version 0.2:   Add StringTie to assembly transcript and TransDecoder to predict Coding sequence For EvidenceModeler


# import python module
import os
from os.path import join
from os.path import abspath
from os.path import basename

# configuration
configfile: "config.yaml"
## repeat makser 
REPEATMASKER = config['REPEATMASKER']
BLAST_BIN    = config['BLAST_BIN']
## ab initio prediction
## TO DO: auto prepare species dir and file
BRAKER       = config["BRAKER"]
Glimmer      = config["Glimmer"]
GlimmerHMM   = config["GlimmerHMM"]
TRANSDECODE  = config['TRANSDECODE']
PASA_DIR     = config['PASA_DIR']
EVIDENCE_DIR = config['EVIDENCE_DIR']
## input file
REFERENCE   = abspath(join("input","genome.fa"))
SPECIES     = config['SPECIES']
SAMPLES     = [line.strip() for line in open(config['SAMPLES'])]
PROTEIN     = abspath(join("input","protein.fa"))

# directory
WORK_DIR        = os.getcwd()
INPUT_DIR       = abspath("input")
PREDICT_DIR     = abspath("prediction")
REPEAT_DIR      = abspath("repeat_masker")
RNA_SEQ_DIR     = abspath("RNA_seq")
TRANSDECODER_DIR= join(PREDICT_DIR, "transdecoder")
TEMP_DIR        = abspath("temp")
LOG_DIR         = abspath("log")

shell("mkdir -p {INPUT_DIR} {PREDICT_DIR} {RNA_SEQ_DIR} {TRANSDECODER_DIR} {REPEAT_DIR} {TEMP_DIR} {LOG_DIR}")

# output file
MASKED_GENOME = join(INPUT_DIR, "genome_hard_mask.fa")
ALL_CLEAN_R1  = [join(RNA_SEQ_DIR, "{}_clean_1.fq".format(sample)) for sample in SAMPLES]
ALL_CLEAN_R2  = [join(RNA_SEQ_DIR, "{}_clean_2.fq".format(sample)) for sample in SAMPLES]
ALL_BAM       = [join(RNA_SEQ_DIR, "{}_align.bam".format(sample)) for sample in SAMPLES]
TRANSDECODE   = join(TRANSDECODER_DIR, "transdecoder.gff3")
GLIMMER_GFF   = join(PREDICT_DIR, "ab_initio", "glimmer.gff")
AUGUSTUS_GFF  = join(PREDICT_DIR, "ab_initio", "augustus", "braker", SPECIES ,"augustus.hints.gff3")
GENEWISE_GFF  = join(PREDICT_DIR, "homology", "genewise.gff")
PASA_GFF      = join(PREDICT_DIR, "PASA", SPECIES+".pasa_assemblies.gff3")
EVM_GFF       = join(PREDICT_DIR, "EVM", "EVM.all.gff")

rule all:
    input:
        MASKED_GENOME, AUGUSTUS_GFF, GENEWISE_GFF, PASA_GFF, EVM_GFF

rule cleanpasa:
    params:
        pasa = join(PREDICT_DIR, "PASA"),
        evm  = join(PREDICT_DIR, "EVM")
    shell:"""
    rm -rf {params.pasa}
    rm -rf {params.evm}
    """

rule cleanevm:
    params:
        evm  = join(PREDICT_DIR, "EVM")
    shell:"""
    rm -rf {params.evm}
    """

# prepare the input from prediction
## genome mask with RepeatMasker
rule hard_mask:
    input: 
        ref = REFERENCE
    params:
        prefix  = basename(REFERENCE),
        engine  = "ncbi",
        species = "arabidopsis",
        outdir  = "repeat_masker"
    threads: 30
    output:
        join(REPEAT_DIR, "genome.fa.masked")
    shell:"""
    {REPEATMASKER} -e {params.engine} -species {params.species} -pa {threads} -gff -dir {params.outdir} {input}
    """
rule mask_filter:
    input:
        join(REPEAT_DIR, "genome.fa.masked")
    output:
        join(INPUT_DIR, "genome_hard_mask.fa")
    shell:"""
    python3 scripts/masked_genome_filter.py {input} > {output}
    """

## rna-seq data clean
rule rna_read_qc:
    input:
        r1 = join(RNA_SEQ_DIR, "{sample}_1.fq"),
        r2 = join(RNA_SEQ_DIR, "{sample}_2.fq")
    output:
        r1 = join(RNA_SEQ_DIR, "{sample}_clean_1.fq"),
        r2 = join(RNA_SEQ_DIR, "{sample}_clean_2.fq")
    threads: 8
    shell:"""
    fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -q 20 -l 100 -w 8
    """

## build index for hisat2
rule build_hisat2_index:
    input:
        ref = join(INPUT_DIR, "genome_hard_mask.fa")
    params:
        prefix = join(RNA_SEQ_DIR, "hisat2_idx")
    output:
        join(RNA_SEQ_DIR, "hisat2_idx.1.ht2")
    threads: 20
    shell:"""
    source activate align
    hisat2-build -p {threads} {input.ref} {params.prefix}
    """
## use hisat2 to align
rule hisat2_align:
    input:
        r1  = join(RNA_SEQ_DIR, "{sample}_clean_1.fq"),
        r2  = join(RNA_SEQ_DIR, "{sample}_clean_2.fq"),
        idx = join(RNA_SEQ_DIR, "hisat2_idx.1.ht2")
    params:
        prefix = join(RNA_SEQ_DIR, "hisat2_idx")
    output:
        join(RNA_SEQ_DIR, "{sample}_align.bam")
    log: join(LOG_DIR, "{sample}_hisat2.log")
    threads: 20
    message: "use hisat2 to align RNA-seq data"
    shell:"""
    source activate align
    hisat2 --dta --new-summary -p {threads} -x {params.prefix} \
        -1 {input.r1} -2 {input.r2} 2> {log} | samtools sort -@ 8 > {output}
    samtools index {output}
    """

## use stringtie to assmebly RNA-seq
rule StringTie:
    input:
        ALL_BAM
    params:
        tmp = join(RNA_SEQ_DIR, "merged.bam")
    output:
        join(RNA_SEQ_DIR, "merged.gtf")
    threads: 6
    message: "merge {ALL_BAM} into single bam and using StringTie to assembly"
    shell:"""
    source activate assembly
    samtools -@ {threads} {params.tmp} {input}
    stringtie -o {ouput} -p {threads} {params.tmp}
    """

rule TransDecoder:
    input:
        ref = REFERENCE,
        gtf = join(RNA_SEQ_DIR, "merged.gtf")
    params:
        wkdir = TRANSDECODER_DIR
    output:
        join(TRANSDECODER_DIR, "transdecoder.gff3")
    shell:"""
    cd {params.wkdir}
    {TRANSDECODE}/util/gtf_genome_to_cdna_fasta.pl {input.gtf} {input.ref} > {params.wkdir}/transcripts.fasta
    {TRANSDECODE}/util/gtf_to_alignment_gff3.pl {input.gtf} > {params.wkdir}/transcripts.gff3
    {TRANSDECODE}/TransDecoder.LongOrfs -t {params.wkdir}/transcripts.fasta
    {TRANSDECODE}/TransDecoder.Predict -t {params.wkdir}/transcripts.fasta
    {TRANSDECODE}/util/cdna_alignment_orf_to_genome_orf.pl \
         {params.wkdir}/transcripts.fasta.transdecoder.gff3 \
         {params.wkdir}/transcripts.gff3 \
         {params.wkdir}/transcripts.fasta > {output}
    """


## use trinity to assembly RNA-seq
rule de_novo_assembly:
    input:
        r1 = ALL_CLEAN_R1,
        r2 = ALL_CLEAN_R2,
    params:
        r1 = ','.join(ALL_CLEAN_R1),
        r2 = ','.join(ALL_CLEAN_R2),
        mem_gb = "64G",
        length = "300",
        outdir = join(RNA_SEQ_DIR, "trinity")
    output:
        join(RNA_SEQ_DIR,"trinity","Trinity.fasta")
    threads: 30
    shell:"""
    source activate assembly
    Trinity --seqType fq --max_memory 64G \
        --left {params.r1} \
        --right {params.r2} \
         --CPU {threads} \
        --min_contig_length {params.length} --output {params.outdir}
    """

# Gene structure prediction using difference methods
## ab inito prediction
rule Augustus:
    input:
        ref = join(INPUT_DIR, "genome_hard_mask.fa"),
        bam = ALL_BAM,
    params:
        bams    = ','.join(ALL_BAM),
        species = SPECIES,
        outdir  = join(PREDICT_DIR,"ab_initio","augustus")
    output:
        join(PREDICT_DIR, "ab_initio", "augustus", "braker", SPECIES ,"augustus.hints.gff3")
    threads: 40
    shell:"""
    mkdir -p {params.outdir}
    cd {params.outdir}
    {BRAKER} --gff3 --cores {threads} --species={params.species} --genome={input.ref} --bam={params.bams} \
    --useexisting
    """

rule GlimmerHMM:
    input:
        ref = join(INPUT_DIR, "genome_hard_mask.fa")
    params:
        tmp = join(TEMP_DIR, "glimmer")
    output:
        join(PREDICT_DIR, "ab_initio", "glimmer.gff")
    threads: 20 
    shell:"""
    seqkit split -f -j {threads} -i {input.ref} -O {params.tmp} 
    ls {params.tmp} | while read id; do echo "{Glimmer}  {params.tmp}/${{id}} -d {GlimmerHMM}" -g; done \
        | parallel -j {threads} > {output}
    """

## Homology based prediction
## blastx
rule blastx:
    input:
        ref = REFERENCE,
        prt = PROTEIN
    params:
        prefix = join(TEMP_DIR, basename(PROTEIN))
    output:
        protein_db = join(TEMP_DIR, basename(PROTEIN) + ".psq"),
        blastx_out = join(PREDICT_DIR, "homology", "blastx_outfmt6.txt")
    threads: 40
    shell:"""
    {BLAST_BIN}/makeblastdb -in {input.prt} -out {params.prefix} -dbtype prot
    {BLAST_BIN}/blastx -query {input.ref} -db {params.prefix} \
        -outfmt 6 -num_threads {threads} \
        -best_hit_overhang 0.25 -best_hit_score_edge 0.25 > {output.blastx_out}
    """

rule blastx2bed:
    input:
        join(PREDICT_DIR, "homology", "blastx_outfmt6.txt")
    output:
        join(PREDICT_DIR, "homology", "homology.bed")
    message: "convert blastx result in format 6 to bed, and merge the nearest entry into single"
    shell:"""
    python3 scripts/blastx2bed.py {input} | bedtools sort -i - > {output}.tmp
    python3 scripts/blastx_bed_merge.py {output}.tmp > {output}
    rm -f {output}.tmp
    """

rule parallel_wise:
    input:
        ref = REFERENCE,
        prt = PROTEIN,
        bed = join(PREDICT_DIR, "homology", "homology.bed")
    output:
        join(PREDICT_DIR, "homology", "genewise.gff")
    threads: 20
    shell:"""
    python3 scripts/parallel_genewise.py {input.ref} {input.prt} {input.bed}
    python3 scripts/merge_parallel_genewise_out.py {input.bed} > {output}
    """

## PASA integration
rule pasa_integration:
    input:
        ref = join(INPUT_DIR, "genome_hard_mask.fa"),
        transcript = join(RNA_SEQ_DIR,"trinity","Trinity.fasta")
    params:
        wkdir    = join(PREDICT_DIR, "PASA"),
        species  = SPECIES,
        align    = "80",
        identity = "80"
    output:
        join(PREDICT_DIR, "PASA", SPECIES+".pasa_assemblies.gff3")
    threads: 20
    shell:"""
    mkdir -p {params.wkdir}
    cd {params.wkdir}   
    cp -f {PASA_DIR}/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config
    sed -i -e 's/<__DATABASE__>/{params.species}/' \
        -e 's/<__MIN_PERCENT_ALIGNED__>/{params.align}/' \
        -e  's/<__MIN_AVG_PER_ID__>/{params.identity}/' alignAssembly.config
    {PASA_DIR}/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R \
        -g {input.ref} -t {input.transcript} --CPU {threads} --ALIGNERS blat,gmap
    """

# using EVM to infer 
rule EvidenceModeler_Prepartion:
    input:
        ref      = REFERENCE,
        genewise = GENEWISE_GFF,
        transdecode = TRANSDECODE,
        glimmer  = GLIMMER_GFF,
        augustus = AUGUSTUS_GFF,
        pasa     = PASA_GFF
    params:
        wkdir      = join(PREDICT_DIR, "EVM"),
        genewise   = "2",
        transdecode= "7",
        glimmer    = "1",
        augustus   = "2",
        pasa       = "6"
    output:
        join(PREDICT_DIR, "EVM", "gene_prediction.gff3"),
        join(PREDICT_DIR, "EVM", "transcript_alignments.gff3"),
        join(PREDICT_DIR, "EVM", "transdecode.gff3"),
        join(PREDICT_DIR, "EVM", "weights.txt") 
    shell:"""
    mkdir -p {params.wkdir}
    python {WORK_DIR}/scripts/glimmer2gff3.py {input.glimmer} > {params.wkdir}/glimmer.gff3 
    python {WORK_DIR}/scripts/wise2gff3.py {input.genewise} > {params.wkdir}/genewise.gff3
    ln -f {input.augustus} {params.wkdir}/augustus.gff3
    ln -f {input.pasa} {params.wkdir}/transcript_alignments.gff3
    ln -f {input.transdecode} {params.wkdir}/transdecode.gff3
    cd {params.wkdir}
    cat {params.wkdir}/augustus.gff3 {params.wkdir}/glimmer.gff3  {params.wkdir}/genewise.gff3  {params.wkdir}/transdecode.gff3 > gene_prediction.gff3
    echo -e 'ABINITIO_PREDICTION\\tAUGUSTUS\\t{params.augustus}\\nABINITIO_PREDICTION\\tGlimmerHMM\\t{params.glimmer}' > weights.txt
    echo -e 'ABINITIO_PREDICTION\\tGeneWise\\t{params.genewise}' >> weights.txt
    echo -e 'TRANSCRIPT\\tassembler-{SPECIES}\\t{params.pasa}' >> weights.txt
    echo -e 'OTHER_PREDICTION\\ttransdecoder\\t{params.transdecode}' >> weights.txt
    """ 

rule EvidenceModeler_Commands:
    input:
        ref    = REFERENCE,
        gene   = join(PREDICT_DIR, "EVM", "gene_prediction.gff3"),
        trans  = join(PREDICT_DIR, "EVM", "transcript_alignments.gff3"),
        weight = join(PREDICT_DIR, "EVM", "weights.txt") 
    params:
        wkdir    = join(PREDICT_DIR, "EVM"),
        segmentSize = 1000000,
        overlapSize = 10000,
    output:
        command  = join(PREDICT_DIR, "EVM", "commands.list")
    shell:"""
    cd {params.wkdir}
    {EVIDENCE_DIR}/EvmUtils/partition_EVM_inputs.pl --genome {input.ref} --gene_predictions {input.gene} \
         --transcript_alignments {input.trans} --segmentSize {params.segmentSize} --overlapSize {params.overlapSize}\
         --partition_listing partitions_list.out 2> /dev/null
    {EVIDENCE_DIR}/EvmUtils/write_EVM_commands.pl --genome {input.ref} --weights {input.weight} \
      --gene_predictions {input.gene} \
      --transcript_alignments {input.trans} \
      --output_file_name evm.out  --partitions partitions_list.out >  commands.list
    """

rule EvidenceModler_Result:
    input:
        ref    = REFERENCE,
        command  = join(PREDICT_DIR, "EVM", "commands.list")
    params:
        wkdir    = join(PREDICT_DIR, "EVM")
    output:
        join(PREDICT_DIR, "EVM", "EVM.all.gff")
    threads: 40
    shell:"""
    cd {params.wkdir}
    parallel --jobs {threads} < {input.command}
    {EVIDENCE_DIR}/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
    {EVIDENCE_DIR}/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out --genome {input.ref}
    find . -regex ".*evm.out.gff3" -exec cat {{}} \; > {output}
    """

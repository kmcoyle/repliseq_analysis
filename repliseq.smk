#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Krysta Coyle
# Module Author:    N/A
# Contributors:     N/A


##### SETUP #####

import os
import oncopipe as op
import pandas


SAMPLES = pandas.read_csv('/projects/rmorin_scratch/krysta_temp/repliseq/samples.tsv', sep = "\t")

##### RULES #####
out_dir = "/projects/rmorin_scratch/krysta_temp/repliseq/"


## these steps are not included in below snakefile
# bedtools makewindows -g hg38.chrom.sizes -w 50000 -s 1000 > 50kb.hg38.bed

# bedtools makewindows -g hg38.chrom.sizes -w 150 > 150bp.hg38.bed


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _bwa_mem_run:
    input:
        fastq_1 = "/projects/rmorin_scratch/krysta_temp/repliseq/00-data/{sample_id}.fastq.gz",
        fasta = "/projects/rmorin_scratch/lcr-modules-references/genomes/hg38/bwa_index/bwa-0.7.17/genome.fa"
    output:
        sam = out_dir + "01-sam/{sample_id}_out.sam"
    log:
        stderr = out_dir + "01-sam/{sample_id}/bwa.stderr.log"
    params:
        bwa_mem = '-M -R "@RG\\tID:{sample_id}\\tLB:{sample_id}\\tPL:ILLUMINA\\tSM:{sample_id}"'
    conda:
        out_dir + "ref/bwa_mem.yaml"
    resources:
        mem_mb = 50000
    threads: 8
    shell:
        op.as_one_line("""
        bwa mem -t {threads} 
        {params.bwa_mem}
        {input.fasta}
        {input.fastq_1} 
        > {output.sam}
        2> {log.stderr}
        """)


rule _bwa_mem_samtools:
    input:
        sam = str(rules._bwa_mem_run.output.sam)
    output:
        bam = out_dir + "02-bam/{sample_id}_out.bam", 
        complete = out_dir + "02-bam/{sample_id}_out.bam.complete"
    log:
        stderr = out_dir + "02-bam/{sample_id}/bwa.stderr.log"
    params:
        opts = "-bhS"
    conda:
        out_dir + "ref/samtools.yaml"
    threads: 1
    shell:
        op.as_one_line("""
        samtools view {params.opts}
        {input.sam} > {output.bam} && touch {output.complete}
        2> {log.stderr}
        """)


rule _bwa_mem_symlink_bam:
    input:
        bam = str(rules._bwa_mem_samtools.output.bam)
    output:
        bam = out_dir + "03-sort_bam/{sample_id}_out.bam", 
    run:
        op.absolute_symlink(input.bam, output.bam)

rule _utils_bam_sort:
    input:
        bam = out_dir + "03-sort_bam/{sample_id}_out.bam"
    output:
        bam = out_dir + "03-sort_bam/{sample_id}.sort.bam",
        prefix = temp(directory(out_dir + "03-sort_bam/{sample_id}_tmp"))
    log:
        stdout = out_dir + "03-sort_bam/{sample_id}_sort.stdout.log",
        stderr = out_dir + "03-sort_bam/{sample_id}_sort.stderr.log"
    params:
        memory = 500
    conda:
        out_dir + "ref/samtools.yaml"
    threads: 12
    resources: 
        mem_mb = 12000
    shell:
        op.as_one_line("""
        mkdir -p {output.prefix} &&
        samtools sort 
        -@ {threads} -m $(({params.memory}))M
        -T {output.prefix} -o {output.bam} 
        {input.bam} 
        > {log.stdout}
        2> {log.stderr}
        """)

rule _bwa_mem_symlink_sorted_bam:
    input:
        bam = out_dir + "03-sort_bam/{sample_id}.sort.bam",
        bwa_mem_bam = str(rules._bwa_mem_samtools.output.bam)
    output:
        bam = out_dir + "04-markdups/{sample_id}.sort.bam",
    run:
        op.absolute_symlink(input.bam, output.bam)
        os.remove(input.bwa_mem_bam)
        shell("touch {input.bwa_mem_bam}.deleted")


# _utils_bam_markdups: Mark duplicates in a BAM file using Picard criteria
rule _utils_bam_markdups:
    input:
        bam = out_dir + "04-markdups/{sample_id}.sort.bam",
    output:
        bam = out_dir + "04-markdups/{sample_id}.mdups.bam",
    log:
        stdout = out_dir + "04-markdups/{sample_id}_mark_dups.stdout.log",
        stderr = out_dir + "04-markdups/{sample_id}_mark_dups.stderr.log"
    params:
        opts = "--overflow-list-size 600000"
    conda:
        out_dir + "ref/sambamba.yaml"
    threads: 12
    resources: 
        mem_mb = 8000
    shell:
        op.as_one_line("""
        sambamba markdup 
        {params.opts} 
        --nthreads {threads} 
        {input.bam} {output.bam} 
        > {log.stdout}
        2> {log.stderr}
        """)


# _utils_bam_index: Index a BAM file
rule _utils_bam_index:
    input:
        bam = out_dir + "04-markdups/{sample_id}.mdups.bam"
    output:
        bam = out_dir + "04-markdups/{sample_id}.mdups.bam.bai"
    log:
        stdout = out_dir + "04-markdups/{sample_id}_index.stdout.log",
        stderr = out_dir + "04-markdups/{sample_id}_index.stderr.log"
    params:
        opts = "-b"
    conda:
        out_dir + "ref/samtools.yaml"
    threads: 6
    resources: 
        mem_mb = 4000
    shell:
        op.as_one_line("""
        samtools index 
        {params.opts} 
        -@ {threads} 
        {input.bam} 
        > {log.stdout}
        2> {log.stderr}
        """)

#samtools view -c -F 260 SAMPLE.bam
rule bam_total_reads:
    input: 
        bam = out_dir + "04-markdups/{sample_id}.mdups.bam"
    output:
        txt = out_dir + "05-reads/{sample_id}.reads.txt"
    conda:
        out_dir + "ref/samtools.yaml"
    resources: 
        mem_mb = 12000
    shell:
        op.as_one_line("""
        samtools view -c -F 260 {input.bam} > {output.txt} 
        """)

rule bam_reads_factor:
    input: 
        txt = out_dir + "05-reads/{sample_id}.reads.txt"
    output:
        txt = out_dir + "05-reads/{sample_id}.factor.txt"
    shell:
        op.as_one_line("""
        var=$(cat {input.txt}) && echo "4000000 / $var" | bc -l > {output.txt}
        """)

## calculate coverage in 150 bp bins
rule bam_150bp_coverage:
    input: 
        bam = out_dir + "04-markdups/{sample_id}.mdups.bam",
        bins = out_dir + "ref/150bp.hg38.bed"
    output:
        bed = out_dir + "06-150bp_coverage/{sample_id}.150bp.hg38.bed",
        complete = out_dir + "06-150bp_coverage/{sample_id}.150bp.complete"
    log:
        stderr = out_dir + "06-150bp_coverage/{sample_id}_150bp.stderr.log"
    conda:
        out_dir + "ref/bedtools.yaml"
    resources: 
        mem_mb = 12000
    shell:
        op.as_one_line("""
        bedtools coverage -b {input.bam} -a {input.bins} -counts > {output.bed} 
        2> {log.stderr}
        && touch {output.complete}
        """)

##apply 4M correction factor
rule bam_150bp_4M_correction:
    input: 
        bed = out_dir + "06-150bp_coverage/{sample_id}.150bp.hg38.bed",
        txt = out_dir + "05-reads/{sample_id}.factor.txt"
    output:
        bed = out_dir + "06-150bp_coverage/{sample_id}.150bp_correct.hg38.bed",
        complete = out_dir + "06-150bp_coverage/{sample_id}.150bp_correct.complete"
    log:
        stderr = out_dir + "06-150bp_coverage/{sample_id}.150bp_correct.stderr.log"
    resources: 
        mem_mb = 12000
    shell:
        op.as_one_line("""
        var=$(cat {input.txt}) && 
        awk -v var="$var" 'BEGIN {{ OFS = "\t" }} ; {{print $1, $2, $3, $4*var}}' {input.bed} > {output.bed} 
        2> {log.stderr}
        && touch {output.complete}
        """)


## edit bins > 5
rule bam_150bp_edit:
    input: 
        bed = out_dir + "06-150bp_coverage/{sample_id}.150bp_correct.hg38.bed",
    output:
        bed = out_dir + "07-150bp_edited/{sample_id}.150bp_correct_edit.hg38.bed",
        complete = out_dir + "07-150bp_edited/{sample_id}.150bp_correct_edit.complete"
    log:
        stderr = out_dir + "07-150bp_edited/{sample_id}.150bp_correct_edit.stderr.log"
    resources: 
        mem_mb = 12000
    shell:
        op.as_one_line("""
        awk 'BEGIN {{ OFS = "\t" }} ; {{ if ($4 > 5) ($4 = 5); print $0}}' {input.bed} > {output.bed} 
        2> {log.stderr}
        && touch {output.complete}
        """)


## calculate coverage in 50kb bins
rule bam_50kb_coverage:
    input: 
        bed = out_dir + "07-150bp_edited/{sample_id}.150bp_correct_edit.hg38.bed",
        bins = out_dir + "ref/50kb.hg38.bed"
    output:
        bed = out_dir + "08-50kb_bins/{sample_id}.50kb.hg38.bed",
        complete = out_dir + "08-50kb_bins/{sample_id}.50kb.complete"
    log:
        stderr = out_dir + "08-50kb_bins/{sample_id}_50kb.stderr.log"
    conda:
        out_dir + "ref/bedtools.yaml"
    resources: 
        mem_mb = 12000
    shell:
        op.as_one_line("""
        bedtools map -b {input.bed} -a {input.bins} -c 4 -o sum > {output.bed} 
        2> {log.stderr}
        && touch {output.complete}
        """)


rule merge_beds:
    input:
        bed = expand(str(rules.bam_150bp_edit.output.bed),
             zip,
             sample_id=SAMPLES['sample_id'])
    output:
        bed = out_dir + "09-merge_bam/merged.bed",
        complete = out_dir + "09-merge_bam/merged.bed.complete"
    log:
        stderr = out_dir + "09-merge_bam/merged.bed.stderr.log"
    conda:
        out_dir + "ref/bedtools.yaml"
    resources: 
        mem_mb = 120000
    shell:
        op.as_one_line("""
        cat {input.bed} | bedtools sort -i stdin | bedtools merge -i stdin -d -1 -c 4 -o sum > {output.bed} 
        2> {log.stderr}
        && touch {output.complete}
        """)


rule merge_bed_50kb_coverage:
    input: 
        bins = out_dir + "ref/50kb.hg38.bed",
        bed = str(rules.merge_beds.output.bed)
    output:
        bed = out_dir + "09-merge_bam/50kb.hg38.bed",
        complete = out_dir + "09-merge_bam/50kb.hg38.complete"
    log:
        stderr = out_dir + "09-merge_bam/50kb.hg38.stderr.log"
    conda:
        out_dir + "ref/bedtools.yaml"
    resources: 
        mem_mb = 250000
    shell:
        op.as_one_line("""
        bedtools sort -i {input.bins} | bedtools map -b {input.bed} -a stdin -c 4 -o sum > {output.bed} 
        2> {log.stderr}
        && touch {output.complete}
        """)

rule all:
    input:
        str(rules.merge_bed_50kb_coverage.output.bed),
        expand(str(rules.bam_50kb_coverage.output.bed),
            zip,
            sample_id=SAMPLES['sample_id'])

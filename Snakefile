configfile: "stacks2_env.yaml"

import pandas as pd

RUNS = ["ac_finland","ac_norway","ac_sweden"]

LANES = dict()

for run in RUNS:
     PATH = "barcodes/ac_scandinavia/" + run + "_barcodes.txt"
     LANES[run] = pd.read_csv(PATH, header=None, sep="\t").iloc[:,2].tolist()

SAMPLES = dict()

for lane, samples in LANES.items():
    for sample in samples:
        SAMPLES[sample] = lane

rule all:
    input:
        "stacks_output/populations.snps.vcf"

rule fastp:
    output:
        read1 = "fastp/{indir}_R1_001.fastq.gz",
        read2 = "fastp/{indir}_R2_001.fastq.gz"
    threads: 20
    input:
        read1 = "raw_reads/ac_scandinavia/{indir}_R1_001.fastq.gz",
        read2 = "raw_reads/ac_scandinavia/{indir}_R2_001.fastq.gz"
    log:
        "logs_fastp/{indir}.log"        
    shell:
        r"""
        fastp -i {input.read1} -I {input.read2} -o {output.read1} -O {output.read2} --html {wildcards.indir}.html &> {log}
        """
 
rule process_radtags:
    output:
       demultiplexed = directory("demultiplexed_{indir}"), 
    input:
        read1 = "fastp/{indir}_R1_001.fastq.gz",
        read2 = "fastp/{indir}_R2_001.fastq.gz",
        barcodes = "barcodes/ac_scandinavia/{indir}_barcodes.txt"
    threads:20
    shell:
        r"""
        mkdir {output.demultiplexed}
        process_radtags -1 {input.read1} -2 {input.read2}  -o {output.demultiplexed} -b {input.barcodes} --threads 20  -c -q --inline_inline --renz_1 sbfI --renz_2 nlaIII
        """

rule samples:
    output:
        read1 = "samples_scandinavia/{sample}.1.fq.gz",
        read2 = "samples_scandinavia/{sample}.2.fq.gz"
    input:
        lambda x: f"demultiplexed_{SAMPLES[x.sample]}",
    shell:
        r"""
        mkdir -p samples_scandinavia
        mv {input}/{wildcards.sample}.1.fq.gz {output.read1}
        mv {input}/{wildcards.sample}.2.fq.gz {output.read2}
        """

rule alignment_sam:
    output:
        sam = temporary("sam/{name}.sam")
    input:
        read1 = "samples_scandinavia/{name}.1.fq.gz",
        read2 = "samples_scandinavia/{name}.2.fq.gz" 
    log:
        "logs_sam/{name}.log"
    threads:20    
    shell:
        r"""
        mkdir -p sam
        bowtie2 -x /home/csas0002/Ref_genomes/arctic_charr -p 16 --phred33 --very-sensitive \
             -1 {input.read1} -2 {input.read2} -S {output.sam} &> {log}
        """

rule alignment_bam:
    output:
        bam = "bam/{name}.bam"
    input:
         sam = "sam/{name}.sam"
    threads: 20     
    shell:
        r"""
        mkdir -p bam
        samtools view -b {input.sam} | samtools sort --threads 16  > {output.bam}
        """

rule bam_index:
    output:
        bam_index = "bam/{name}.bam.bai"
    input:
        bam = "bam/{name}.bam"
    threads:20
    shell:
        r"""
        samtools index -@ 4 {input.bam}
        """

rule gstakcs:
    output:
        gstacks = directory("stacks_output"),
        catalog =  "stacks_output/catalog.fa.gz",
        gstack_details = "stacks_output/gstacks.details.gz"
    input:
        map = "populations_map/ac_scandinavia/ac_scandinavia_map.txt",
        bai = expand("bam/{name}.bam.bai", name=SAMPLES.keys()),
        bam = expand("bam/{name}.bam", name=SAMPLES.keys())
    threads: 20
    shell:
        r"""
        gstacks -I bam -M {input.map} -O {output.gstacks}  -t 16 --gt-alpha 0.001 --details --min-mapq 20
        """

rule populations:
    output:
        snps = "{indir}/populations.snps.vcf"
    input: 
        "{indir}/catalog.fa.gz",
        "{indir}/gstacks.details.gz",
        map = "populations_map/ac_scandinavia/ac_scandinavia_map.txt"
    threads: 16
    shell:
        r"""
        populations -P {wildcards.indir} -M {input.map} -t 16 -r 0.8 --min-maf 0.05 --min-mac 3 \
            --max-obs-het 0.6 --write-single-snp --smooth-fstats --smooth-popstats --vcf 
        """ 






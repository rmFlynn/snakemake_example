

rule bowtie2_build_database:
    input:
       config['inputs']['bowtie2_database_raw']
    output:
       directory(config['inputs']['bowtie2_database_built'])
    threads: workflow.cores
    # log:
    #     f"logs/database{os.path.getsize(config['inputs']['bowtie2_database_raw'])}/bowtie2_build_database.log"
    params:
        larg_index='--large-index',
    benchmark:
        f"benchmarks/bowtie2_build_database/database{os.path.getsize(config['inputs']['bowtie2_database_raw'])}.tsv"
    shell:
       """
       mkdir {output}
       bowtie2-build \\
           {params.larg_index} \\
           {input} \\
           {output}/bowtie_DB \\
           --threads {threads}
       """


rule bowtie2_map:
    input:
       bowtie_database=config['inputs']['bowtie2_database_built'],
       r1="results/{sample}_bined_R1.fastq",
       r2="results/{sample}_bined_R2.fastq"
    output:
       sam_file=temp("results/bowtie2/{sample}.sam")
    threads: workflow.cores
    # log:
    #     "logs/{sample}/bowtie2_map.log"
    benchmark:
        "benchmarks/bowtie2_map/{sample}.tsv"
    shell:
       """
       bowtie2 \\
           -D 10 -R 2 -N 0 -L 22 -i S,0,2.50 \\
           -p {threads} \\
           -x {input.bowtie_database}/bowtie_DB \\
           -S {output.sam_file} \\
           -1 {input.r1} \\
           -2 {input.r2}
       """

rule samtools_sam_to_bam:
    input:
       "results/counting/{sample}.final.sam",
    output:
       "results/{sample}.bam"
    threads: workflow.cores
    # log:
    #     "logs/{sample}/samtools_sam_to_bam.log"
    benchmark:
        "benchmarks/samtools_sam_to_bam/{sample}.tsv"
    shell:
       """
      samtools view -bS -@ {threads} {input} > {output}
       """

# rule samtools_sam_to_bam:
#     input:
#        sam_file="results/bowtie2/{sample}.sam"
#     output:
#        bam_file=temp("results/bowtie2/{sample}.bam")
#     threads: workflow.cores
#     # log:
#     #     "logs/{sample}/samtools_sam_to_bam.log"
#     benchmark:
#         "benchmarks/samtools_sam_to_bam/{sample}.tsv"
#     shell:
#        """
#       samtools view -bS -@ {threads} {input} > {output}
#        """


# rule reformat_bam:
#     input:
#        bam_in="results/bowtie2/{sample}.bam"
#     output:
#        bam_out=temp("results/bowtie2/{sample}_filtered.bam")
#     threads: workflow.cores
#     # log:
#     #     "logs/{sample}/reformat_bam.log"
#     benchmark:
#         "benchmarks/reformat_bam/{sample}.tsv"
#     shell:
#        """
#        reformat.sh \\
#            in={input.bam_in}\\
#            out={output.bam_out} \\
#            minidfilter=1 \\
#            primaryonly=t pairedonly=f
#        """


# rule samtools_bam_sort:
#     input:
#        bam_in="results/bowtie2/{sample}_filtered.bam"
#     output:
#        bam_out="results/bowtie2/{sample}_sorted_filtered.bam"
#     threads: workflow.cores
#     benchmark:
#         "benchmarks/samtools_bam_sort/{sample}.tsv"
#     # log:
#     #     "logs/{sample}/samtools_bam_sort.log"
#     shell:
#        """
#        samtools sort -n -@ {threads} -o {output} {input}
#        """



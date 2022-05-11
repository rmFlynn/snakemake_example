
ruleorder: reformat_sam > samtools_sam_sort > rich_script_sam_file > counting > join_counting > samtools_sam_to_bam


rule reformat_sam:
    input:
       sam_in="results/bowtie2/{sample}.sam"
    output:
       sam_out=temp("results/counting/{sample}.reformatted.sam")
    threads: workflow.cores
    # log:
    #     "logs/{sample}/reformat_sam.log"
    benchmark:
        "benchmarks/reformat_sam/{sample}.tsv"
    shell:
       """
       reformat.sh \\
           in={input.sam_in}\\
           out={output.sam_out} \\
           minidfilter=1 \\
           primaryonly=t pairedonly=f
       """


rule samtools_sam_sort:
    input:
       sam_in="results/counting/{sample}.reformatted.sam"
    output:
       sam_out=temp("results/counting/{sample}.sorted.reformatted.sam")
    threads: workflow.cores
    benchmark:
        "benchmarks/samtools_sam_to_bam/{sample}.tsv"
    # log:
    #     "logs/{sample}/samtools_sam_sort.log"
    shell:
       """
       samtools sort -n -@ {threads} -o {output} {input}
       """


rule rich_script_sam_file:
    input:
       sam_in="results/counting/{sample}.sorted.reformatted.sam"
    output:
       sam_out=temp("results/counting/{sample}.sorted.reformatted.filtered.sam")
    threads: workflow.cores
    params:
       percent_read=config['counting']['percent_read'],
       max_mismatch=config['counting']['max_mismatch'],
       mult_align='T'
    benchmark:
        "benchmarks/rich_script_sam_file/{sample}.tsv"
    shell:
       """
       /usr/bin/python /ORG-Data/scripts/sam_file.py \\
                       -i {input.sam_in} \\
                       -o {output.sam_out} \\
                       --percent_read {params.percent_read} \\
                       --max_mismatch {params.max_mismatch} \\
                       --mult_align {params.mult_align}
       """


rule counting:
    input:
       sam="results/counting/{sample}.sorted.reformatted.filtered.sam",
       gff=config['inputs']['gff']
    output:
        temp('results/counting/{sample}.counts.tsv')
    threads: workflow.cores
    params:
       percent_read=config['counting']['percent_read'],
       featuretype=config['counting']['featuretype'],
       idattr=config['counting']['idattr'],
       stranded='reverse',
       nonunique='Rory'
    benchmark:
        "benchmarks/counting/{sample}.tsv"
    shell:
       """
       python3 workflow/scripts/flexcount.py \\
                       {input.sam} \\
                       {input.gff} \\
                       --type {params.featuretype} \\
                       --idattr  {params.idattr} \\
                       --stranded {params.stranded} \\
                       --nonunique {params.nonunique} \\
                       --counts_output {output}
       """


rule join_counting:
    input:
       expand("results/counting/{sample}.counts.tsv", sample=sample_dict.keys())
    output:
      "results/all_counts.tsv"
    threads: 
       workflow.cores
    benchmark:
        "benchmarks/join_counting.tsv"
    run:
        shell(f"cat {input[0]} > {{output}}")
        for i in input[1:]:
            shell(f"sed 1,1d {i} >> {{output}}")


rule samtools_sam_to_bam:
    input:
       "results/counting/{sample}.sorted.reformatted.filtered.sam",
    output:
       temp("results/{sample}.bam")
    threads: workflow.cores
    # log:
    #     "logs/{sample}/samtools_sam_to_bam.log"
    benchmark:
        "benchmarks/samtools_sam_to_bam/{sample}.tsv"
    shell:
       """
      samtools view -bS -@ {threads} {input} > {output}
       """

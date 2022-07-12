
# ruleorder: reformat_sam > samtools_sam_sort > rich_script_sam_file > counting > join_counting > samtools_sam_to_bam

import os
from glob import glob
SAMPLES = [os.path.basename(i)[:-11] for i in glob("/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaT_mapping/bowtie_mapping_final/*_sorted.bam")]


# rule reformat_sam:
#     input:
#        sam_in="results/bowtie2/{sample}.sam"
#     output:
#        sam_out=temp("results/counting/{sample}.reformatted.sam")
#     threads: workflow.cores
#     # log:
#     #     "logs/{sample}/reformat_sam.log"
#     benchmark:
#         "benchmarks/reformat_sam/{sample}.tsv"
#     shell:
#        """
#        reformat.sh \\
#            in={input.sam_in}\\
#            out={output.sam_out} \\
#            minidfilter=1 \\
#            primaryonly=t pairedonly=f
#        """
# 
# 
# rule samtools_sam_sort:
#     input:
#        sam_in="results/counting/{sample}.reformatted.sam"
#     output:
#        sam_out=temp("results/counting/{sample}.sorted.reformatted.sam")
#     threads: workflow.cores
#     benchmark:
#         "benchmarks/samtools_sam_to_bam/{sample}.tsv"
#     # log:
#     #     "logs/{sample}/samtools_sam_sort.log"
#     shell:
#        """
#        samtools sort -n -@ {threads} -o {output} {input}
#        """
# 
# 
# rule rich_script_sam_file:
#     input:
#        sam_in="results/counting/{sample}.sorted.reformatted.sam"
#     output:
#        sam_out=temp("results/counting/{sample}.final.sam")
#     threads: workflow.cores
#     params:
#        percent_read=config['counting']['percent_read'],
#        max_mismatch=config['counting']['max_mismatch'],
#        mult_align='T'
#     benchmark:
#         "benchmarks/rich_script_sam_file/{sample}.tsv"
#     shell:
#        """
#        /usr/bin/python /ORG-Data/scripts/sam_file.py \\
#                        -i {input.sam_in} \\
#                        -o {output.sam_out} \\
#                        --percent_read {params.percent_read} \\
#                        --max_mismatch {params.max_mismatch} \\
#                        --mult_align {params.mult_align}
#        """

# rule bam_to_sam:
#     input:
#        bam="/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaT_mapping/bowtie_mapping_final/{sample}_sorted.bam",
#     output:
#        sam=temp("results/counting/{sample}.final.sam"),
#     resources:
#         mem=1000,
#         time='1-00:00:00'
#     threads: 1
#     shell:
#             "samtools view -h {input.bam} > {output.sam}"

rule counting:
    input:
       sam="/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaT_mapping/bowtie_mapping_final/{sample}_sorted.bam",
       gff=config['inputs']['gff']
    output:
        temp('results/counting/{sample}.counts.tsv')
    resources:
        mem=4000,
        time='1-00:00:00'
    threads: 1
    # conda:
    # we can't have nice things untill pysam (depencacy)
    # gets its act together
    #         "../envs/counting.yaml"
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
       expand("results/counting/{sample}.counts.tsv", sample=SAMPLES)
    output:
      "results/all_counts.tsv"
    threads: 
       1
    benchmark:
        "benchmarks/join_counting.tsv"
    run:
        shell(f"cat {input[0]} > {output}")
        for i in input[1:]:
            shell(f"sed 1,1d {i} >> {output}")



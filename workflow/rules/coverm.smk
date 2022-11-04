"""
Here we both run coverM and combine the output using a more complex algorithum, to make the coverm_trimmed_zero data all the data with 0s inserted and the coverm_trimed_filterd data sets, the option to make one or both is in the config, in the coverm section. 
After coverm joiner joins the coverm runs for a sample the coverm_merger merges the sample into one ouput file 
"""
COVERM_ENV = '../envs/coverm.yaml'

# BAM = "/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaG_mapping/mapping_bowtie_012722/{sample}.bam"
# FASTA_DIR = '/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/DRAM_renamed_bins' 

# SAMPLES=[ os.path.basename(i)[:-4] for i in 
#                 glob("/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaG_mapping/mapping_bowtie_012722/*sorted.bam")]
# wildcard_constraints:
#     sample = '|'.join(SAMPLES)

rule coverm_trimmed_mean:
    input: 
        bam = "results/{sample}.bam",
        fasta_dir = config['coverm']['fasta_dir']
    output:
        mean = temp("results/coverm/{sample}_coverm_trimmed_mean.tsv"),
        stats = temp("results/coverm/stats_{sample}_coverm_trimmed_mean.tsv")
    benchmark:
        "benchmarks/coverm_trimmed_mean/{sample}.txt"
    conda:
        COVERM_ENV
    threads: 
        workflow.cores
    params:
        min_read_percent_identity_pair = config['coverm']['min_read_percent_identity_pair'],
        genome_fasta_extension = config['coverm']['genome_fasta_extension']
    shell:
        """
           coverm genome \\
               --proper-pairs-only \\
               --genome-fasta-extension {params.genome_fasta_extension} \\
               --genome-fasta-directory {input.fasta_dir} \\
               --bam-files {input.bam} \\
               --threads {threads} \\
               --min-read-percent-identity-pair {params.min_read_percent_identity_pair} \\
               -m trimmed_mean \\
               --output-file {output.mean} &> {output.stats} 
        """

rule coverm_reads_per_base:
    input: 
        bam = "results/{sample}.bam",
        fasta_dir = config['coverm']['fasta_dir']
    output:
        base = temp("results/coverm/{sample}_coverm_reads_per_base.tsv"),
        stats = temp("results/coverm/stats_{sample}_coverm_reads_per_base.tsv")
    conda:
        COVERM_ENV
    benchmark:
        "benchmarks/coverm_reads_per_base/{sample}.txt"
    threads: 
        workflow.cores
    params:
        min_read_percent_identity_pair= config['coverm']['min_read_percent_identity_pair'],
        min_covered_fraction=0,  # must be 0 for this config['coverm']['min_covered_fraction'],
        genome_fasta_extension= config['coverm']['genome_fasta_extension']
    shell:
        """
           coverm genome \\
               --proper-pairs-only \\
               --genome-fasta-extension {params.genome_fasta_extension} \\
               --genome-fasta-directory {input.fasta_dir} \\
               --bam-files {input.bam} \\
               --threads {threads} \\
               --min-read-percent-identity-pair {params.min_read_percent_identity_pair} \\
               --min-covered-fraction {params.min_covered_fraction} \\
               -m reads_per_base \\
               --output-file {output.base} &> {output.stats}
        """



rule coverm_covered_fraction:
    input: 
        bam = "results/{sample}.bam",
        fasta_dir = config['coverm']['fasta_dir']
    output:
        frac = temp("results/coverm/{sample}_coverm_covered_fraction.tsv"),
        stats = temp("results/coverm/stats_{sample}_coverm_covered_fraction.tsv")
    benchmark:
            "benchmarks/coverm_covered_fraction/{sample}.txt"
    threads:
        workflow.cores
    conda: 
        COVERM_ENV
    params:
        min_read_percent_identity_pair = config['coverm']['min_read_percent_identity_pair'],
        min_covered_fraction = config['coverm']['min_covered_fraction'],
        genome_fasta_extension = config['coverm']['genome_fasta_extension']
    shell:
        """
           coverm genome \\
               --proper-pairs-only \\
               --genome-fasta-extension {params.genome_fasta_extension} \\
               --genome-fasta-directory  {input.fasta_dir} \\
               --bam-files {input.bam} \\
               --threads {threads} \\
               --min-read-percent-identity-pair {params.min_read_percent_identity_pair} \\
               --min-covered-fraction {params.min_covered_fraction} \\
               -m mean \\
               --output-file {output.frac} &> {output.stats}

        """


rule coverm_join:
    input: 
        mean = "results/coverm/{sample}_coverm_trimmed_mean.tsv",
        frac = "results/coverm/{sample}_coverm_covered_fraction.tsv",
        base = "results/coverm/{sample}_coverm_reads_per_base.tsv",
        bam = "results/{sample}.bam"
    output:
        zero = "results/coverm/{sample}_coverm_trimmed_zeroed.tsv",
        filt = "results/coverm/{sample}_coverm_trimmed_filtered.tsv"
    benchmark:
        "benchmarks/coverm_join/{sample}.txt"
    params:
        min_depth = config['coverm']['min_depth'],
        sequence_lenth= config['coverm']['sequence_lenth'],
        calculate_per_seq_length= config['coverm']['calculate_per_seq_length']
    script:
        '../scripts/coverm_joiner.py'

rule coverm_merge:
    input:
        expand("results/coverm/{sample}_coverm_trimmed_zeroed.tsv", 
               sample=SAMPLE_DICT.keys())
    output:
        "results/coverm_all_output_merged.tsv"
    run:
        pd.concat([pd.read_csv(i, sep='\t', index_col='Genome') for i in input], axis=1).to_csv(output[0])


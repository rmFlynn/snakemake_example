"""
Here we both run coverM and combine the output to make the coverm_trimmed_zero data and the coverm_trimed_filterd data sets, the option to make one or both is in the config, it the checkM section. 

"""
rule coverm_trimmed_mean:
    input: 
        bam:"{sample}.bam",
        fasta_dir: '/home/projects-wrighton/NIH_Salmonella/Salmonella/Metagenomes/MAG_database/MQ_HQ_bins/bins/MQHQ_fastas/scaffolds_renamed' 
    output:
        mean = temp("results/coverm_trimmed_mean/{sample}_coverm_trimmed_mean.tsv"),
        sats = temp("results/coverm_trimmed_mean/stats_{sample}_coverm_trimmed_mean.tsv")
    conda: "some_env"
    params:
        min_read_percent_identity_pair: config['coverm']['min_read_percent_identity_pair'],
        genome-fasta-extension: config['coverm']['genome-fasta-extension']
    shell:
        """
           coverm genome \\
               --proper-pairs-only \\
               --genome-fasta-extension {params.genome-fasta-extension} \\
               --genome-fasta-directory {input.fasta_dir} \\
               --bam-files} {input.bam} \\
               --threads {threads} \\
               --min-read-percent-identity-pair {params.min_read_percent_identity_pair} \\
               -m trimmed_mean \\
               --output-file {output.mean} &> {output.stats} 
        """


rule coverm_reads_per_base:
    input: 
        bam:"{sample}.bam",
        fasta_dir: '/home/projects-wrighton/NIH_Salmonella/Salmonella/Metagenomes/MAG_database/MQ_HQ_bins/bins/MQHQ_fastas/scaffolds_renamed' 
    output:
        base = temp("results/coverm_reads_per_base/{sample}_coverm_reads_per_base.tsv"),
        sats = temp("results/coverm_reads_per_base/stats_{sample}_coverm_reads_per_base.tsv")
    conda: "some_env"
    params:
        min_read_percent_identity_pair: config['coverm']['min_read_percent_identity_pair'],
        min_covered_fraction: config['coverm']['min_covered_fraction'],
        genome-fasta-extension: config['coverm']['genome-fasta-extension']
    shell:
        """
           coverm genome \\
               --proper-pairs-only \\
               --genome-fasta-extension {params.genome-fasta-extension} \\
               --genome-fasta-directory {input.fasta_dir} \\
               --bam-files {input.bam} \\
               --threads {threads} \\
               --min-read-percent-identity-pair {params.min_read_percent_identity_pair} \\
               --min-covered-fraction {params.min_covered_fraction} \\
               -m reads_per_base \\
               --output-file {output.mean} &> {output.stats}
        """


rule coverm_covered_fraction:
    input: 
        bam:"{sample}.bam",
        fasta_dir: '/home/projects-wrighton/NIH_Salmonella/Salmonella/Metagenomes/MAG_database/MQ_HQ_bins/bins/MQHQ_fastas/scaffolds_renamed' 
    output:
        frac = temp("results/coverm_covered_fraction/{sample}_coverm_covered_fraction.tsv"),
        sats = temp("results/coverm_covered_fraction/stats_{sample}_coverm_covered_fraction.tsv")
    conda: "some_env"
    threads:
        workflow.cores
    params:
        min_read_percent_identity_pair: config['coverm']['min_read_percent_identity_pair'],
        min_covered_fraction: config['coverm']['min_covered_fraction'],
        genome-fasta-extension: config['coverm']['genome-fasta-extension']
    shell:
        """
           coverm genome \\
               --proper-pairs-only \\
               --genome-fasta-extension {params.genome-fasta-extension} \\
               --genome-fasta-directory  {input.fasta_dir} \\
               --bam-files {input.bam} \\
               --threads {threads} \\
               --min-read-percent-identity-pair {params.min_read_percent_identity_pair} \\
               --min-covered-fraction {params.min_covered_fraction} \\
               -m mean \\
               --output-file {output.mean} &> {output.stats}

        """



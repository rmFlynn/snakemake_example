import gzip
import shutil
from workflow.scripts.singel_thread_test import fastas_dup_check


def tempfunk2(w):
    print(w.sample)
    return 


rule unzip_inputfiles:
    input:
        lambda w:  SAMPLE_DICT[w.sample][w.type + '_gz'],
    wildcard_constraints:
            type = "R1,R2,inter"
    output:
        read="results/raw_files/{sample}/raw_{type}.fastq",
    resources:
        mem=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 10,
        time='1-00:00:00'
    run:
        shell("mkdir -p output.outdir")
        with gzip.open(input[0], 'rb') as f_in:
            with open(output.read, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)


rule deinterleave_fastq_inputfiles:
    input:
        il=lambda w: SAMPLE_DICT[w.sample]['inter'],
    output:
        outdir=temp(directory("results/raw_files/{sample}")),
        r1="results/raw_files/{sample}/raw_R1.fastq",
        r2="results/raw_files/{sample}/raw_R2.fastq"
    resources:
        mem=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 10,
        time='1-00:00:00'
    shell:
        """
        mkdir -p output.outdir
        /ORG-Data/scripts/deinterleave_fastq.sh < {input.il} {output.r1} {output.r1}
        """


rule check_dups_fasta:
    input:
        reads=expand("{fasta}", 
                        fasta=[i for j in SAMPLE_DICT.values() for i in [j['R1'], j['R2']]])
    output:
        "results/duplicate_names_note.txt"
    resources:
        mem=1000,
        time='1-00:00:00'
    benchmark:
        "benchmarks/duplicate/check"
    run:
        if not config['backend']['check_fasta_for_dups']:
            msg='Duplicate test waved'
        else:
            if fastas_dup_check(input.reads, '@'):
                msg='There are no duplicates'
        with open(output[0], 'w') as f:
            f.write(msg)


rule fastqc_get_quality:
    input:
        r1=lambda w: SAMPLE_DICT[w.sample]['R1'] if w.quality == 'raw' else "results/{sample}_bined_R1.fastq",
        r2=lambda w: SAMPLE_DICT[w.sample]['R2'] if w.quality == 'raw' else "results/{sample}_bined_R2.fastq",
        no_dups="results/duplicate_names_note.txt"
    output:
       folder = directory("results/read_quality/{sample}_{quality}"),
    wildcard_constraints:
            raw = "raw|trimmed_filtered"
    resources:
        mem=2570,
        time='1-00:00:00'
    params:
       folder_temp ="results/read_quality/{sample}_{quality}_temp",
    benchmark:
        "benchmarks/fastqc/{sample}_{quality}.tsv"
    shell:
       """
       mkdir -p {params.folder_temp}
       fastqc {input.r1} {input.r2} -o {params.folder_temp}
       mkdir -p {output.folder}
       mv {params.folder_temp}/*.zip {output.folder}
       rm -r {params.folder_temp}
       """

rule bbduk_trim_reads:
    input:
        r1=lambda w: SAMPLE_DICT[w.sample]['R1'],
        r2=lambda w: SAMPLE_DICT[w.sample]['R2']
    output:
       spath=temp(directory("results/binning/{sample}_bbduk")),
       out1 = "results/binning/{sample}_bbduk/trimmed_R1.fastq",
       out2 = "results/binning/{sample}_bbduk/trimmed_R2.fastq"
    threads: 
        workflow.cores
    resources:
        mem=lambda wildcards, input, attempt: (input.size//100000000) * attempt * 10,
        time='1-00:00:00'
    # log:
    #     "logs/{sample}/bbduk_trim_reads.log"
    benchmark:
        "benchmarks/bbduk_trim_reads/{sample}.tsv"
    shell:
       """
       bbduk.sh threads={threads} \\
           overwrite=t \\
           in1={input.r1} \\
           in2={input.r2} \\
           ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=20 minlength=75 \\
           ref=/opt/bbtools/bbmap/resources/adapters.fa  \\
           out1={output.out1} \\
           out2={output.out2}
       """


rule sickle_trimreads:
    input:
        r1=lambda w: SAMPLE_DICT[w.sample]['R1'],
        r2=lambda w: SAMPLE_DICT[w.sample]['R2']
    output:
        file_path_rna=temp(directory("results/binning/{sample}_sickel")),
        r1="results/binning/{sample}_sickel/trimmed_R1.fastq",
        r2="results/binning/{sample}_sickel/trimmed_R2.fastq",
        singles = "results/binning/{sample}_sickel/trimmed_singles.fastq"
    threads: 
        1
    resources:
        mem= 7,
        time='1-00:00:00'
    benchmark:
        "benchmarks/sickle_trim_reads/{sample}_bin.tsv"
    params:
        sickle_quality_type = config['binning']['sickle_quality_type']
    shell:
       """
       sickle pe \\
           -f {input.r1} \\
           -r {input.r2} \\
           -t {params.sickle_quality_type} \\
           -o {output.r1} \\
           -p {output.r2} \\
           -s {output.singles}
       """


rule rqcfilter2_filter_reads:
    input:
       file = f"results/binning/{{sample}}_{'bbduk' if config['type'] == 'meta_t' else 'sickel'}",
       in1  = f"results/binning/{{sample}}_{'bbduk' if config['type'] == 'meta_t' else 'sickel'}/trimmed_R1.fastq",
       in2  = f"results/binning/{{sample}}_{'bbduk' if config['type'] == 'meta_t' else 'sickel'}/trimmed_R2.fastq"
    output:
       file_path=temp(directory("results/filtered_trimed_reads/{sample}")),
       file_path_rna=temp("results/filtered_trimed_reads/{sample}/RNA.fq.gz"),
       interleaved=temp("results/filtered_trimed_reads/{sample}/trimmed_R1.anqrpht.fastq.gz")
    threads: workflow.cores
    resources:
       mem=309062, # Needs tweekinglambda wildcards, input, attempt: (input.size//1000000000) * attempt * 10,
       time='14-00:00:00'
    benchmark:
        "benchmarks/rqcfilter2_filter_reads/{sample}.tsv"
    shell:
       """
       mkdir -p {output.file_path}
       rqcfilter2.sh jni=t \\
           threads={threads} \\
           in1={input.in1} \\
           in2={input.in2} \\
           path={output.file_path} \\
           outribo={output.file_path_rna} \\
           rna=t trimfragadapter=t qtrim=r trimq=0 maxns=1 \\
           maq=10 minlen=51 mlf=0.33 phix=t removeribo=t \\
           removehuman=t removedog=t removecat=t removemouse=t \\
           khist=t removemicrobes=t mtst=t sketch kapa=t \\
           clumpify=t tmpdir=null barcodefilter=f trimpolyg=5 \\
           -Xmx300g  rqcfilterdata=/home/opt/RQCFilterData
       """

rule unzip_file:
    # input: "{path}.gz"
    input: 
       #  rules.rqcfilter2_filter_reads.output
       file_path="results/filtered_trimed_reads/{sample}",
       file_path_rna="results/filtered_trimed_reads/{sample}/RNA.fq.gz",
       interleaved="results/filtered_trimed_reads/{sample}/trimmed_R1.anqrpht.fastq.gz"
    output: 
        temp("results/filtered_trimed_reads/{sample}/trimmed_R1.anqrpht_unziped.fastq"),
    resources:
        mem=lambda wildcards, input, attempt: (input.size//10000000) * attempt * 10,
        time='1-00:00:00'
    # log:
    #     "logs/{sample}/unzip_filtered_file.log"
    run:
        with gzip.open(input['interleaved'], 'rb') as f_in:
            with open(output[0], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

rule deinterleave_fastq:
    input:
       file_path="results/filtered_trimed_reads/{sample}",
       file="results/filtered_trimed_reads/{sample}/trimmed_R1.anqrpht_unziped.fastq"
    output:
       r1=temp("results/{sample}_bined_R1.fastq"),
       r2=temp("results/{sample}_bined_R2.fastq")
    resources:
        mem=10, # This is invariant to the best I can tell
        time='1-00:00:00'
    # log:
    #     "logs/{sample}/deinterleave_fastq.log"
    benchmark:
        "benchmarks/deinterleave_fastq/{sample}.tsv"
    shell:
       """
       bash /ORG-Data/scripts/deinterleave_fastq.sh < {input.file} {output.r1}  {output.r2}
       """


"""
If you have separate files for forward and reverse reads:
Usage: sickle pe [options] -f <paired-end forward fastq file> -r <paired-end reverse fastq file> -t <quality type> -o <trimmed PE forward file> -p <trimmed PE reverse file> -s <trimmed singles file>

If you have one file with interleaved forward and reverse reads:
Usage: sickle pe [options] -c <interleaved input file> -t <quality type> -m <interleaved trimmed paired-end output> -s <trimmed singles file>

If you have one file with interleaved reads as input and you want ONLY one interleaved file as output:
Usage: sickle pe [options] -c <interleaved input file> -t <quality type> -M <interleaved trimmed output>
Paired-end separated reads
--------------------------
-f, --pe-file1, Input paired-end forward fastq file (Input files must have same number of records)
-r, --pe-file2, Input paired-end reverse fastq file
-o, --output-pe1, Output trimmed forward fastq file
-p, --output-pe2, Output trimmed reverse fastq file. Must use -s option.

Paired-end interleaved reads
----------------------------
-c, --pe-combo, Combined (interleaved) input paired-end fastq
-m, --output-combo, Output combined (interleaved) paired-end fastq file. Must use -s option.
-M, --output-combo-all, Output combined (interleaved) paired-end fastq file with any discarded read written to output file as a single N. Cannot be used with the -s option.

Global options
--------------
-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)
-s, --output-single, Output trimmed singles fastq file
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.
-x, --no-fiveprime, Don't do five prime trimming.
-n, --truncate-n, Truncate sequences at position of first N.
-g, --gzip-output, Output gzipped files.
--quiet, do not output trimming info
--help, display this help and exit
--version, output version information and exit
"""

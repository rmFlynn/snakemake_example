import gzip
import shutil
from workflow.scripts.singel_thread_test import fastas_dup_check

def tempfunk2(w):
    print(w.sample)
    return 


rule unzip_inputfiles:
    input:
        r1=lambda w: SAMPLE_DICT[w.sample]['forward_gz'],
        r2=lambda w: SAMPLE_DICT[w.sample]['reversed_gz']
    output:
        outdir=temp(directory("results/raw_files/{sample}")),
        r1="results/raw_files/{sample}/raw_R1.fastq",
        r2="results/raw_files/{sample}/raw_R2.fastq"
    # log:
    #     "logs/{sample}/unzip_inputfiles.log"
    resources:
        mem=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 10,
        time='1-00:00:00'
    run:
        shell("mkdir -p output.outdir")
        with gzip.open(input.r1, 'rb') as f_in:
            with open(output.r1, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        with gzip.open(input.r2, 'rb') as f_in:
            with open(output.r2, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
 

rule check_dups_fasta:
    input:
        reads=expand("{fasta}", 
                        fasta=[i for j in SAMPLE_DICT.values() for i in [j['forward'], j['reversed']]])
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
        r1=lambda w: SAMPLE_DICT[w.sample]['forward'] if w.quality == 'raw' else "results/{sample}_bined_R1.fastq",
        r2=lambda w: SAMPLE_DICT[w.sample]['reversed'] if w.quality == 'raw' else "results/{sample}_bined_R2.fastq",
       
       
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
        r1=lambda w: SAMPLE_DICT[w.sample]['forward'],
        r2=lambda w: SAMPLE_DICT[w.sample]['reversed']
    output:
       spath=temp(directory("results/trimmed_reads/{sample}"))
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
           out1={output.spath}/trimmed_R1.fastq \\
           out2={output.spath}/trimmed_R2.fastq
       """


rule rqcfilter2_filter_reads:
    input:
       "results/trimmed_reads/{sample}"
    output:
       file_path=temp(directory("results/filtered_trimed_reads/{sample}")),
       file_path_rna=temp("results/filtered_trimed_reads/{sample}/RNA.fq.gz"),
       interleaved=temp("results/filtered_trimed_reads/{sample}/trimmed_R1.anqrpht.fastq.gz")
    threads: workflow.cores
    resources:
        mem=100000, #lambda wildcards, input, attempt: (input.size//1000000000) * attempt * 10,
        time='14-00:00:00'
    # log:
    #     "logs/{sample}/rqcfilter2_filter_reads.log"
    benchmark:
        "benchmarks/rqcfilter2_filter_reads/{sample}.tsv"
    shell:
       """
       rqcfilter2.sh jni=t \\
           threads={threads} \\
           in1={input}/trimmed_R1.fastq \\
           in2={input}/trimmed_R2.fastq \\
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
        mem=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 10,
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
        mem=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 10,
        time='1-00:00:00'
    # log:
    #     "logs/{sample}/deinterleave_fastq.log"
    benchmark:
        "benchmarks/deinterleave_fastq/{sample}.tsv"
    shell:
       """
       bash /ORG-Data/scripts/deinterleave_fastq.sh < {input.file} {output.r1}  {output.r2}
       """



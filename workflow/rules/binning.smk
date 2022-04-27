
wildcard_constraints:
    sample="[A-z,0-9,_,-]+"

def tempfunk2(w):
    print(w.sample)
    return 

rule unzip_inputfiles:
    input:
        r1=lambda w: sample_dict[w.sample]['forward_gz'],
        r2=lambda w: sample_dict[w.sample]['reversed_gz']
    output:
        outdir=directory("results/raw_files/{sample}"),
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
 

rule fastqc_raw:
    input:
        r1=lambda w: sample_dict[w.sample]['forward'],
        r2=lambda w: sample_dict[w.sample]['reversed']
    output:
       folder = directory("results/raw_read_quality/{sample}"),
       report_html_r1=report("results/raw_read_quality/{sample}/raw_R1_fastqc.html", category='Sample Quality'),
       report_html_r2=report("results/raw_read_quality/{sample}/raw_R2_fastqc.html", category='Sample Quality')
    resources:
        mem=2570,
        time='1-00:00:00'
    # log:
    #     "logs/{sample}/fastqc_raw.log"
    benchmark:
        "benchmarks/fastqc_raw/{sample}.tsv"
    shell:
       """
       fastqc {input.r1} {input.r2} -o {output}
       """


rule bbduk_trim_reads:
    input:
        r1=lambda w: sample_dict[w.sample]['forward'],
        r2=lambda w: sample_dict[w.sample]['reversed']
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
       file_path=directory("results/filtered_trimed_reads/{sample}"),
       file_path_rna="results/filtered_trimed_reads/{sample}/RNA.fq.gz",
       interleaved="results/filtered_trimed_reads/{sample}/trimmed_R1.anqrpht.fastq.gz"
    threads: workflow.cores
    resources:
        mem=lambda wildcards, input, attempt: (input.size//10000000) * attempt * 10,
        time='1-00:00:00'
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

rule unzip_filtered_file:
    # input: "{path}.gz"
    input: 
        "results/filtered_trimed_reads/{sample}/trimmed_R1.anqrpht.fastq.gz"
    output: 
        temp("results/filtered_trimed_reads/{sample}/trimmed_R1.anqrpht_unziped.fastq")
    resources:
        mem=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 10,
        time='1-00:00:00'
    # log:
    #     "logs/{sample}/unzip_filtered_file.log"
    run:
        with gzip.open(input[0], 'rb') as f_in:
            with open(output[0], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

rule deinterleave_fastq:
    input:
       "results/filtered_trimed_reads/{sample}/trimmed_R1.anqrpht_unziped.fastq"
    output:
       r1=temp("results/filtered_trimed_reads_deinterleave/{sample}/R1.fastq"),
       r2=temp("results/filtered_trimed_reads_deinterleave/{sample}/R2.fastq")
    resources:
        mem=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 10,
        time='1-00:00:00'
    # log:
    #     "logs/{sample}/deinterleave_fastq.log"
    benchmark:
        "benchmarks/deinterleave_fastq/{sample}.tsv"
    shell:
       """
       bash /ORG-Data/scripts/deinterleave_fastq.sh < {input} {output.r1}  {output.r2}
       """


rule fastqc_trimm_filter:
    input:
       r1="results/filtered_trimed_reads_deinterleave/{sample}/R1.fastq",
       r2="results/filtered_trimed_reads_deinterleave/{sample}/R2.fastq"
    output:
       directory("results/trimmed_filtered_read_quality/{sample}"),
       report_html_r1=report("results/trimmed_filtered_read_quality/{sample}/R1_fastqc.html", category='Sample Quality'),
       report_html_r2=report("results/trimmed_filtered_read_quality/{sample}/R2_fastqc.html", category='Sample Quality')
    resources:
        mem=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 10,
        time='1-00:00:00'
    # log:
    #     "logs/{sample}/fastqc_trimm_filter.log"
    benchmark:
        "benchmarks/fastqc_trimm_filter/{sample}.tsv"
    shell:
       """
       mkdir -p {output[0]}
       fastqc {input.r1} {input.r2} -o {output}
       """

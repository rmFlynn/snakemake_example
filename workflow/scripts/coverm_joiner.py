import pandas as pd
import click
import subprocess
from io import StringIO


# os.system(
# COUNT_OPTION_1 = """
# echo "fasta\tlength" > bin_lengths.txt
# for fasta in /home/projects-wrighton/NIH_Salmonella/Salmonella/Metagenomes/MAG_database/MQ_HQ_bins/bins/MQHQ_fastas/scaffolds_renamed/*fa
# do
#   length=`grep -v ">" $fasta | wc | awk '{print $3-$1}'`
#   name=`basename $fasta`
#   echo "${name}\t${length}" >> bin_lengths.txt
# done
# """
# COUNT_OPTION_2 = 

# os.system("""
# echo "fasta\tlength" > bin_lengths.txt
# for fasta in /home/projects-wrighton/NIH_Salmonella/Salmonella/Metagenomes/MAG_database/MQ_HQ_bins/bins/MQHQ_fastas/scaffolds_renamed/*fa
# do
#   length=`grep -v ">" $fasta | wc | awk '{print $3-$1}'`
#   name=`basename $fasta`
#   echo "${name}\t${length}" >> bin_lengths.txt
# done
# """)
# min_depth = 3
# read_per_base_path = './data/JLP_reads_per_base.tsv'
# covered_fract_path = './data/JLP_covered_fraction.tsv'
# coverm_trimed_path = './data/JLP_trimmed_mean.tsv'

def get_bin_lengths(bam:str):
    with subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE) as sam:
        sam_lines = [i for i in sam.stdout]
    split_lines = [str(i).split('\\t') for i in  sam_lines]
    bad_lines = [i for i in  split_lines if len(i) < 10]
    if len(bad_lines) > 0:
        breakpoint()
    lenths = [(i[9]) for i in split_lines]
    names = [i[2] for i in split_lines]
    data = pd.DataFrame({'length': lenths}, index=names)
    return data

def combine_coverm(coverm_trimed_path, read_per_base_path,
                   covered_fract_path, bam, out_zero, out_filterd,
                   min_depth, sequence_lenth, calculate_per_seq_length,
                   untrimmed_output):
    read_per_base = pd.read_csv(read_per_base_path, 
                                sep='\t', index_col='Genome').fillna(0)
    covered_fract = pd.read_csv(covered_fract_path, 
                                sep='\t', index_col='Genome').fillna(0)
    coverm_trimed = pd.read_csv(coverm_trimed_path, 
                                sep='\t', index_col='Genome').fillna(0)

    read_per_base.T[read_per_base.apply(lambda x: x.isna().any(), axis=0)]
    if isinstance(bam, str):
        bam = [bam]
    if calculate_per_seq_length:
        # this should never hapen but if it dose it will be slow
        lengths = get_bin_lengths(bam[0])
        read_per_base = read_per_base.merge(lengths, left_index=True, 
                                            right_index=True, how='right')
        read_per_base = read_per_base.iloc[:,:-1].apply(
            lambda x : x * read_per_base['length'])
    else:
        read_per_base = read_per_base * sequence_lenth

    coverm_trimed[covered_fract == 0] = 0
    coverm_trimed[read_per_base >= min_depth] = 0

    coverm_trimed.to_csv(out_zero, sep='\t')
    keep = ~(coverm_trimed==0)
    coverm_trimed.loc[
        keep.any(axis=1), 
        keep.any(axis=0)].to_csv(out_filterd, sep='\t')

combine_coverm(snakemake.input['mean'], snakemake.input['base'], 
               snakemake.input['frac'], snakemake.input['bam'], 
                snakemake.output['zero'],
                snakemake.output['filt'],
                snakemake.params['min_depth'],
                snakemake.params['sequence_lenth'],
                snakemake.params['calculate_per_seq_length'],
                untrimmed_output=True)

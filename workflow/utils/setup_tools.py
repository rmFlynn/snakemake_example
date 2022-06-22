import os
from glob import glob


def ungz(path):
    if path.endswith('.gz'):
        return path.replace('.gz', '')
    return path


def clean(string, r1:str, r2:str):
    removable = [r1, r2, '.fastq', '.gz', '.fq']
    for i in removable:
         string = string.replace(i, '')
    return string
   

def find_match(name, files, r1:str, r2:str):
    for i in files:
        if clean(i, r1, r2) == name:
            return i
   

def get_reads(paths:list, r1:str, r2:str):
    # Check for dupes in fastas

    freads = [i for i in paths if r1 in i]
    breads = [i for i in paths if r2 in i]
    if len(freads) != len(breads):
        raise ValueError("Somehow there is a read that has no pair."
                        f"the forward reads are: {freads}"
                        f"the reverse reads are: {breads}"
                        )
    names = [clean(i, r1, r2) for i in freads]
    paths = {
        os.path.basename(i):{'R1_gz': find_match(i, freads, r1, r2), 
           'R2_gz': find_match(i, breads, r1, r2)}
        for i in names}
    for i in paths:
        if (gzf:=paths[i]['R1_gz']).endswith('.gz'):
            paths[i]['R1'] = f"results/raw_files/{i}/raw_R1.fastq"
        else:
            paths[i]['R1'] = gzf
        if (gzr:=paths[i]['R2_gz']).endswith('.gz'):
            paths[i]['R2'] = f"results/raw_files/{i}/raw_R2.fastq"
        else:
            paths[i]['R2'] = gzr
        paths[i]['inter'] = None
        paths[i]['inter_gz'] = None
    return  paths


def get_inter_reads(paths:list):
    # Check for dupes in fastas

    names = [clean(i, "", "") for i in freads]
    paths = { os.path.basename(i):{'inter_gz': freads}
             for i in names}
    for i in paths:
        if (gzf:=paths[i]['inter_gz']).endswith('.gz'):
            paths[i]['inter'] = f"results/raw_files/{i}/raw_inter.fastq"
        else:
            paths[i]['inter'] = gzf
        paths[i]['R1'] = f"results/raw_files/{i}/raw_R2.fastq"
        paths[i]['R2'] = f"results/raw_files/{i}/raw_R2.fastq"
        paths[i]['R1_gz'] = None
        paths[i]['R2_gz'] = None
    return  paths


def get_named_reads(named_paths):
    def make_paths(k:str, v:dict):
        out:dict = {}
        r1 = v.get('R1')
        r2 = v.get('R2')
        inter = v.get('inter')
        if r1 is not None and r1.endswith('gz'):
            out['R1_gz'] = r1
            out['R1'] = f"results/raw_files/{k}/raw_R1.fastq"
        else:
            out['R1_gz'] = None
            out['R1'] = r1
        if r2 is not None and r2.endswith('gz'):
            out['R2_gz'] = r2
            out['R2'] = f"results/raw_files/{k}/raw_R2.fastq"
        else:
            out['R2_gz'] = None
            out['R2'] = r2
        if inter is not None and inter.endswith('gz'):
            out['inter_gz'] = inter
            out['inter'] = f"results/raw_files/{k}/raw_inter.fastq"
        else:
            out['inter_gz'] = None
            out['inter'] = inter
        return out
    return {k: make_paths(k, v) for k, v in named_paths.items()}


def get_sample_dict(config:dict):
    r1 = config['binning']['forward_id']
    r2 = config['binning']['backward_id']
    sample_dict = {}
    if (inter_reads := config['inputs'].get('interleaved_reads')) is not None:
        paths = [i for p in inter_reads for i in glob(p)]
        sample_dict.update(get_inter_reads(paths, r1, r2))
    if (paired_reads := config['inputs'].get('paired_reads')) is not None:
        paths = [i for p in paired_reads for i in glob(p)]
        sample_dict.update(get_reads(paths, r1, r2))
    if (named_paths := config['inputs'].get('named_reads')) is not None:
        sample_dict.update(get_named_reads(named_paths))
    return sample_dict

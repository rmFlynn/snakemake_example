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
        os.path.basename(i):{'forward_gz': find_match(i, freads, r1, r2), 
           'reversed_gz': find_match(i, breads, r1, r2)}
        for i in names}
    for i in paths:
        if (gzf:=paths[i]['forward_gz']).endswith('.gz'):
            paths[i]['forward'] = f"results/raw_files/{i}/raw_R1.fastq"
        else:
            paths[i]['forward'] = gzf
        if (gzr:=paths[i]['reversed_gz']).endswith('.gz'):
            paths[i]['reversed'] = f"results/raw_files/{i}/raw_R2.fastq"
        else:
            paths[i]['reversed'] = gzr
    return  paths


def get_named_reads(named_paths):
    def make_paths(k:str, v:dict):
        out:dict = {}
        if v['forward'].endswith('gz'):
            out['forward_gz'] = v['forward']
            out['forward'] = f"results/raw_files/{k}/raw_R1.fastq"
        else:
            out['forward_gz'] = None
            out['forward'] = v['forward']
        if v['reversed'].endswith('gz'):
            out['reversed_gz'] = v['reversed']
            out['reversed'] = f"results/raw_files/{k}/raw_R1.fastq"
        else:
            out['reversed_gz'] = None
            out['reversed'] = v['reversed']
        return out
    return {k: make_paths(k, v) for k, v in named_paths.items()}


def get_sample_dict(config:dict):
    r1 = config['binning']['forward_id']
    r2 = config['binning']['backward_id']
    if (paired_reads := config['inputs'].get('paired_reads')) is not None:
        paths = [i for p in paired_reads for i in glob(p)]
        sample_dict = get_reads(paths, r1, r2)
    else:
        sample_dict = {}
    if (named_paths := config['inputs'].get('named_reads')) is not None:
        sample_dict.update(get_named_reads(named_paths))
    return sample_dict

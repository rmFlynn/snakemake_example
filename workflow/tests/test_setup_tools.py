from workflow.scripts.setup_tools import get_sample_dict

def test_get_sample_dict():
    assert get_sample_dict({ 
        'inputs':{
            'paired_reads':{
                "workflow/tests/data/test_finder*.fa"
            },
            'named_reads':{
                'any_name_you_want': {
                    'forward': "workflow/tests/data/test_finder_R1.fa",
                    'reversed': "workflow/tests/data/test_finder_R2.fa"
                }
            }
        },
        'binning': {
            'forward_id': '_R1',
            'backward_id': '_R2'
        }
    }) =={'test_finder.fa': {'forward_gz': 'workflow/tests/data/test_finder_R1.fa', 'reversed_gz': 'workflow/tests/data/test_finder_R2.fa', 'forward': 'workflow/tests/data/test_finder_R1.fa', 'reversed': 'workflow/tests/data/test_finder_R2.fa'}, 'any_name_you_want': {'forward_gz': None, 'forward': 'workflow/tests/data/test_finder_R1.fa', 'reversed_gz': None, 'reversed': 'workflow/tests/data/test_finder_R2.fa'}}

    assert get_sample_dict({ 
        'inputs':{
            'paired_reads':{
                "workflow/tests/data/test_finder*.fa.gz"
            },
            'named_reads':{
                'any_name_you_want': {
                    'forward': "workflow/tests/data/test_finder_R1.fa.gz",
                    'reversed': "workflow/tests/data/test_finder_R2.fa.gz"
                }
            }
        },
        'binning': {
            'forward_id': '_R1',
            'backward_id': '_R2'
        }
    }) == {'test_finder.fa': {'forward_gz': 'workflow/tests/data/test_finder_R1.fa.gz', 'reversed_gz': 'workflow/tests/data/test_finder_R2.fa.gz', 'forward': 'results/raw_files/test_finder.fa/raw_R1.fastq', 'reversed': 'results/raw_files/test_finder.fa/raw_R2.fastq'}, 'any_name_you_want': {'forward_gz': 'workflow/tests/data/test_finder_R1.fa.gz', 'forward': 'results/raw_files/any_name_you_want/raw_R1.fastq', 'reversed_gz': 'workflow/tests/data/test_finder_R2.fa.gz', 'reversed': 'results/raw_files/any_name_you_want/raw_R1.fastq'}}

    assert get_sample_dict({ 
        'inputs':{
            'named_reads':{
                'any_name_you_want': {
                    'forward': "workflow/tests/data/test_finder_r1.fa.gz",
                    'reversed': "workflow/tests/data/test_finder_r2.fa.gz"
                }
            }
        },
        'binning': {
            'forward_id': '_R1',
            'backward_id': '_R2'
        }
    }) == {'any_name_you_want': {'forward_gz': 'workflow/tests/data/test_finder_r1.fa.gz', 'forward': 'results/raw_files/any_name_you_want/raw_R1.fastq', 'reversed_gz': 'workflow/tests/data/test_finder_r2.fa.gz', 'reversed': 'results/raw_files/any_name_you_want/raw_R1.fastq'}}

    assert get_sample_dict({ 
        'inputs':{
            'paired_reads':{
                "workflow/tests/data/test_finder*.fa.gz"
            }
        },
        'binning': {
            'forward_id': '_R1',
            'backward_id': '_R2'
        }
    }) == {'test_finder.fa': {'forward_gz': 'workflow/tests/data/test_finder_R1.fa.gz', 'reversed_gz': 'workflow/tests/data/test_finder_R2.fa.gz', 'forward': 'results/raw_files/test_finder.fa/raw_R1.fastq', 'reversed': 'results/raw_files/test_finder.fa/raw_R2.fastq'}}

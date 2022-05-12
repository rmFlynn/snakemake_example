import os
import pytest
from workflow.scripts.singel_thread_test import fastas_dup_check

def test_fastas_dup_check():
    assert fastas_dup_check([os.path.join('workflow', 'tests', 'data', i) for i in ['1.fa', '2.fa', '3.fa']])
    with pytest.raises(ValueError, match=r"contains "
                       "duplicate headers, you must correct this before continuing. "
                       "The duplicate headers are: \['this'\]"):
        fastas_dup_check([os.path.join('workflow', 'tests', 'data', i) for i in ['1.fa', '2.fa', '3.fa', 'dup.fa']])

        
    with pytest.raises(ValueError, match=r"contains duplicate headers, you must correct this before continuing."
                       " The duplicate headers are \['two'\]"):
        fastas_dup_check([os.path.join('workflow', 'tests', 'data', i) for i in ['1.fa', '2.fa', '3.fa', '1dup.fa']])

from typing import Dict

from tools.resources import PURINES, PYRIMIDINES, COMPLEMENT_DNA, COMPLEMENT_RNA


def is_purine(base: str) -> bool:
    return base.upper() in PURINES


def is_pyrimidine(base: str) -> bool:
    return base.upper() in PYRIMIDINES


def is_transition(base1: str, base2: str) -> bool:
    return (base1 != base2)\
           and ((is_purine(base1) and is_purine(base2))
                or (is_pyrimidine(base1) and is_pyrimidine(base2)))


def is_transversion(base1: str, base2: str) -> bool:
    return (base1 != base2)\
           and ((is_purine(base1) and is_pyrimidine(base2))
                or (is_pyrimidine(base1) and is_purine(base2)))


def is_complement(base1: str, base2: str, type: str = 'DNA') -> bool:
    comp_dict: Dict[str, str]
    if type == 'DNA':
        comp_dict = COMPLEMENT_DNA
    elif type == 'RNA':
        comp_dict = COMPLEMENT_RNA
    else:
        raise ValueError('type must be RNA or DNA')
    return comp_dict[base1.upper()] == base2.upper()
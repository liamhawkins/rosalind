from typing import Dict, List, Tuple

from tools.Fasta import Fasta
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


def overlap(f1: Fasta, f2: Fasta, overlap: int) -> bool:
    return f1[-overlap:] == f2[:overlap]


def max_overlap(f1: Fasta, f2: Fasta) -> int:
    min_size: int = min([len(f1), len(f2)])  # Get length of shortest seq
    max_: int = 0
    for i in range(1, min_size + 1):
        if overlap(f1, f2, overlap=i):
            max_ = i
    return max_


def join_reads(f1: Fasta, f2: Fasta, min_overlap: int = 0) -> Fasta:
    max_o = max_overlap(f1, f2)
    if max_o < min_overlap:
        raise ValueError(f'Reads do not overlap by min_overlap: {min_overlap}')
    return Fasta(id_=f'{f1.id}_{f2.id}', sequence=f'{f1.sequence}{f2.sequence[max_o:]}')


def bin_search(i: int, sorted_array: List[int], abs_start_index: int = 0) -> int:
    mid_index = len(sorted_array) // 2
    mid_value = sorted_array[mid_index]
    if len(sorted_array) == 1 and i != mid_value:
        return -1
    elif i == mid_value:
        return mid_index + abs_start_index + 1
    elif i < mid_value:
        return bin_search(i, sorted_array[:mid_index], abs_start_index=abs_start_index)
    else:
        return bin_search(i, sorted_array[mid_index:], abs_start_index=abs_start_index + mid_index)


def ins_sort(array: List[int]) -> Tuple[List[int], int]:
    swap_count: int = 0
    for i in range(1, len(array)):
        k: int = i
        while k > 0 and array[k] < array[k - 1]:
            array[k - 1], array[k] = array[k], array[k - 1]  # Swap
            swap_count += 1
            k = k - 1
    return array, swap_count


def hamming(s: str, t: str) -> int:
    num: int = 0
    for x, y in zip(s, t):
        num += x != y
    return num


if __name__ == '__main__':
    sorted_array = [0,1,2,3,4,5,6,7]
    print(bin_search(8, sorted_array))
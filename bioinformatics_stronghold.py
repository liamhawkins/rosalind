from typing import List, Tuple, Dict, Set

from tools.Fasta import Fasta
from tools.functions import is_transition, is_transversion
from tools.resources import CODON_MAP


def count_bases(s: str) -> Tuple[int, int, int, int]:
    """
    Given: A DNA string s of length at most 1000 nt.
    Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G',
    and 'T' occur in s.
    """
    return s.count('A'), s.count('C'), s.count('G'), s.count('T')


def transcribe(t: str) -> str:
    """
    Given: A DNA string t having length at most 1000 nt.
    Return: The transcribed RNA string of t.
    """
    return t.replace('T', 'U')


def reverse_complement(s: str) -> str:
    """
    Given: A DNA string s of length at most 1000 bp.
    Return: The reverse complement s^c of s.
    """
    comp: Dict[str, str] = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
    }
    return ''.join([comp[base] for base in reversed(s)])


def num_point_mutations(s: str, t: str) -> int:
    """
    Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).
    Return: The Hamming distance dH(s,t).
    """
    num: int = 0
    for x, y in zip(s, t):
        num += x != y
    return num


def max_gc_content(s: str) -> None:
    """
    Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).
    Return: The ID of the string having the highest GC-content, followed by the GC-content of that string.
    """
    fastas: List[Fasta] = Fasta.from_str(s)
    max_gc: Fasta = Fasta(id_='null', sequence='ATAT')
    for f in fastas:
        if f.gc > max_gc.gc:
            max_gc = f
    print(max_gc.id)
    print(max_gc.gc)


def motif_indices(s: str, t: str) -> List[int]:
    """
    Given: Two DNA strings s and t (each of length at most 1 kbp).
    Return: All locations of t as a substring of s.
    """
    locations: List[int] = []
    for index in range(len(s)-len(t)):
        if t == s[index:index+len(t)]:
            locations.append(index + 1)
    return locations


def translate_protein(s: str) -> None:
    """
    Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).
    Return: The protein string encoded by s.
    """
    fasta: Fasta = Fasta('test', s)
    print(fasta.translate_first_orf())


def calculating_molecular_weight(s: str) -> float:
    """
    Given: A protein string P of length at most 1000 aa.
    Return: The total weight of P. Consult the monoisotopic mass table.
    """
    fasta: Fasta = Fasta('test', s)
    return fasta.molecular_weight()


def mendels_first_law(k: int, m: int, n: int) -> float:
    """
    Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: k individuals
    are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.
    Return: The probability that two randomly selected mating organisms will produce an individual possessing a
    dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.
    """
    total: int = k + m + n
    k_k: float = (k/total) * ((k-1)/(total-1)) * 1.0
    k_m: float = (k/total) * (m/(total-1)) * 1.0
    k_n: float = (k/total) * (n/(total-1)) * 1.0

    m_k: float = (m/total) * (k/(total-1)) * 1.0
    m_m: float = (m/total) * ((m-1)/(total-1)) * 0.75
    m_n: float = (m/total) * (n/(total-1)) * 0.5

    n_k: float = (n/total) * (k/(total-1)) * 1.0
    n_m: float = (n/total) * (m/(total-1)) * 0.5
    n_n: float = (n/total) * ((n-1)/(total-1)) * 0.0

    return sum([k_k, k_m, k_n, m_k, m_m, m_n, n_k, n_m, n_n])


def protein_to_RNA(s: str) -> int:
    """
    Given: A protein string of length at most 1000 aa.
    Return: The total number of different RNA strings from which the protein could have been translated,
    modulo 1,000,000. (Don't neglect the importance of the stop codon in protein translation.)
    """
    num: int = 1
    codon_list: List[str] = list(CODON_MAP.values())
    for c in s:
        num *= codon_list.count(c)
    num *= codon_list.count('Stop')
    return num % 1_000_000


def rabbits_and_recurrance(n: int, s: int) -> int:
    """
    Given: Positive integers n≤40 and k≤5.
    Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each
    generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
    """
    if n == 0:
        return 0
    elif n == 1:
        return 1
    else:
        return rabbits_and_recurrance(n-1, s) + rabbits_and_recurrance(n-2, s) * s


def open_reading_frames(s: str) -> Set[str]:
    """
    Given: A DNA string s of length at most 1 kbp in FASTA format.
    Return: Every distinct candidate protein string that can be translated from ORFs of s. Strings can be returned in
    any order.
    """
    fasta: Fasta = Fasta.from_str(s)[0]
    return {y.sequence for x in fasta.translate_all_orfs().values() for y in x}


def RNA_splicing(s: str) -> None:
    """
    Given: A DNA string s (of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings
    are given in FASTA format.
    Return: A protein string resulting from transcribing and translating the exons of s. (Note: Only one solution will
    exist for the dataset provided.)
    """
    fastas = Fasta.from_str(s)
    main = fastas.pop(0)
    for f in fastas:
        main = main.simple_splice(f)
    for f in main.translate_first_orf():
        print(f)


def finding_a_spliced_motif(inp: str) -> None:
    """
    Given: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.
    Return: One collection of indices of s in which the symbols of t appear as a subsequence of s. If multiple solutions exist, you may return any one.
    """
    s: str
    t: str
    s, t = [x.sequence for x in Fasta.from_str(inp)]

    first_index_group: List[str] = []
    last_ind: int = 0
    for x in t:
        ind: int = s.find(x, last_ind+1) + 1
        first_index_group.append(str(ind))
        last_ind = ind
    print(' '.join(first_index_group))


def transition_transversion_ratio(s: Fasta, t: Fasta) -> float:
    """
    Given: Two DNA strings s1 and s2 of equal length (at most 1 kbp).
    Return: The transition/transversion ratio R(s1,s2).
    """
    transition: int = 0
    transversion: int = 0
    for b1, b2 in zip(s.sequence, t.sequence):
        transition += is_transition(b1, b2)
        transversion += is_transversion(b1, b2)
    return transition/transversion


def calculate_expected_offspring(s):
    """
    Given: Six nonnegative integers, each of which does not exceed 20,000. The integers correspond to the number of
    couples in a population possessing each genotype pairing for a given factor. In order, the six given integers
    represent the number of couples having the following genotypes:
    AA-AA
    AA-Aa
    AA-aa
    Aa-Aa
    Aa-aa
    aa-aa
    Return: The expected number of offspring displaying the dominant phenotype in the next generation, under the
    assumption that every couple has exactly two offspring.
    """
    s = [int(x) for x in s.split()]
    return (s[0] + s[1] + s[2] + (s[3]*0.75) + (s[4]*0.5)) * 2


if __name__ == '__main__':
    print(calculate_expected_offspring("19777 17806 19543 19861 19821 19960"))

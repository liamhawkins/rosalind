from typing import List, Tuple, Dict, Set, Optional

from tools.Fasta import Fasta
from tools.Graph import Graph
from tools.functions import is_transition, is_transversion, overlap, max_overlap, join_reads, hamming, get_consensus, \
    pdistance
from tools.resources import CODON_MAP


def counting_DNA_nucleotides(s: str) -> str:
    """
    Given: A DNA string s of length at most 1000 nt.

    Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G',
    and 'T' occur in s.
    """
    t: Tuple[int, ...] = (s.count('A'), s.count('C'), s.count('G'), s.count('T'))
    return ' '.join([str(x) for x in t])


def transcribing_DNA_into_RNA(t: str) -> str:
    """
    Given: A DNA string t having length at most 1000 nt.

    Return: The transcribed RNA string of t.
    """
    return t.replace('T', 'U')


def complementing_a_strand_of_DNA(s: str) -> str:
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


def counting_point_mutations(s: str) -> int:
    """
    Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

    Return: The Hamming distance dH(s,t).
    """
    t: str
    s, t = [x.strip() for x in s.split('\n')]
    return hamming(s, t)


def computing_gc_content(s: str) -> str:
    """
    Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

    Return: The ID of the string having the highest GC-content, followed by the GC-content of that string.
    """
    fastas: List[Fasta] = Fasta.from_str(s)
    max_gc: Fasta = Fasta(id_='null', sequence='ATAT')
    for f in fastas:
        if f.gc > max_gc.gc:
            max_gc = f
    return f'{max_gc.id}\n{max_gc.gc:.6f}'


def finding_a_motif_in_DNA(s: str) -> str:
    """
    Given: Two DNA strings s and t (each of length at most 1 kbp).

    Return: All locations of t as a substring of s.
    """
    t: str
    s, t = [x.strip() for x in s.split('\n')]
    locations: List[str] = []
    for index in range(len(s) - len(t)):
        if t == s[index:index + len(t)]:
            locations.append(str(index + 1))
    return ' '.join(locations)


def translating_RNA_into_protein(s: str) -> str:
    """
    Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).

    Return: The protein string encoded by s.
    """
    fasta: Fasta = Fasta('test', s)
    proteins = fasta.translate_first_orf()
    return proteins[0].sequence  # There can be multiple proteins so return the first


def calculating_protein_mass(s: str) -> float:
    """
    Given: A protein string P of length at most 1000 aa.

    Return: The total weight of P. Consult the monoisotopic mass table.
    """
    fasta: Fasta = Fasta('test', s)
    return fasta.molecular_weight()


def mendels_first_law(s: str) -> float:
    """
    Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: k individuals
    are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.

    Return: The probability that two randomly selected mating organisms will produce an individual possessing a
    dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.
    """
    k: int
    m: int
    n: int
    k, m, n = (int(x) for x in s.split())
    total: int = k + m + n
    k_k: float = (k / total) * ((k - 1) / (total - 1)) * 1.0
    k_m: float = (k / total) * (m / (total - 1)) * 1.0
    k_n: float = (k / total) * (n / (total - 1)) * 1.0

    m_k: float = (m / total) * (k / (total - 1)) * 1.0
    m_m: float = (m / total) * ((m - 1) / (total - 1)) * 0.75
    m_n: float = (m / total) * (n / (total - 1)) * 0.5

    n_k: float = (n / total) * (k / (total - 1)) * 1.0
    n_m: float = (n / total) * (m / (total - 1)) * 0.5
    n_n: float = (n / total) * ((n - 1) / (total - 1)) * 0.0

    return sum([k_k, k_m, k_n, m_k, m_m, m_n, n_k, n_m, n_n])


def inferring_mRNA_from_protein(s: str) -> int:
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


def rabbits_and_recurrance_relations(s: str) -> int:
    """
    Given: Positive integers n≤40 and k≤5.

    Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each
    generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
    """

    def breed(n: int, k: int):
        if n == 0:
            return 0
        elif n == 1:
            return 1
        else:
            return breed(n - 1, k) + breed(n - 2, k) * k

    n: int
    k: int
    n, k = (int(x) for x in s.split())
    return breed(n, k)


def open_reading_frames(s: str) -> Set[str]:
    """
    Given: A DNA string s of length at most 1 kbp in FASTA format.

    Return: Every distinct candidate protein string that can be translated from ORFs of s. Strings can be returned in
    any order.
    """
    fasta: Fasta = Fasta.from_str(s)[0]
    return {y.sequence for x in fasta.translate_all_orfs().values() for y in x}


def RNA_splicing(s: str) -> str:
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
    return main.translate_first_orf()[0].sequence


def finding_a_spliced_motif(inp: str) -> str:
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
        ind: int = s.find(x, last_ind + 1) + 1
        first_index_group.append(str(ind))
        last_ind = ind
    return ' '.join(first_index_group)


def transition_transversion_ratio(s: str) -> float:
    """
    Given: Two DNA strings s1 and s2 of equal length (at most 1 kbp).

    Return: The transition/transversion ratio R(s1,s2).
    """
    fastas = Fasta.from_str(s)

    transition: int = 0
    transversion: int = 0
    for b1, b2 in zip(fastas[0].sequence, fastas[1].sequence):
        transition += is_transition(b1, b2)
        transversion += is_transversion(b1, b2)
    return transition / transversion


def calculate_expected_offspring(s: str) -> float:
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
    i: List[int] = [int(x) for x in s.split()]
    return (i[0] + i[1] + i[2] + (i[3] * 0.75) + (i[4] * 0.5)) * 2


def overlap_graphs(s):
    fastas: List[Fasta] = Fasta.from_str(s)
    edges: List[str] = []
    for f1 in fastas:
        for f2 in fastas:
            if f1 != f2 and overlap(f1, f2, 3):
                edges.append(f'{f1.id} {f2.id}')
    return '\n'.join(edges)


def genome_assembly_as_shortest_superstring(s: str) -> str:
    fastas: List[Fasta] = Fasta.from_str(s)
    min_overlap: int = min([len(x) for x in fastas]) // 2

    def overlap(n1, n2):
        return n1 != n2 and max_overlap(n1, n2) > min_overlap

    g = Graph(directed=True)
    g.add_nodes(fastas)
    g.join_nodes_by_func(overlap)

    if g.is_linear():
        order = [g.starting_nodes()[0]]

        while len(order) < len(g.edges) + 1:
            order.append(order[-1].outgoing_nodes[0])

        ret = order[0].obj
        for f in order[1:]:
            ret = join_reads(ret, f.obj)
        return ret.sequence
    else:
        raise ValueError('Genome cannot be assembled from non-linear graph')


def completing_a_tree(s: str) -> int:
    """
    Given: A positive integer n (n≤1000) and an adjacency list corresponding to a graph on n nodes that contains no
    cycles.

    Return: The minimum number of edges that can be added to the graph to produce a tree.
    """
    g: Graph = Graph.from_str(s)
    return len(g.disconnected_subgraphs()) - 1


def error_correction_in_reads(s: str) -> str:
    """
    Given: A collection of up to 1000 reads of equal length (at most 50 bp) in FASTA format. Some of these reads were
    generated with a single-nucleotide error. For each read s in the dataset, one of the following applies:
    s was correctly sequenced and appears in the dataset at least twice (possibly as a reverse complement);
    s is incorrect, it appears in the dataset exactly once, and its Hamming distance is 1 with respect to exactly one
    correct read in the dataset (or its reverse complement).

    Return: A list of all corrections in the form "[old read]->[new read]". (Each correction must be a single symbol
    substitution, and you may return the corrections in any order.)
    """
    fastas: List[Fasta] = Fasta.from_str(s)
    rev_comp_fastas: List[Fasta] = [f.reverse_complement() for f in fastas]
    error_reads: List[Fasta] = list()
    for f in fastas:
        if fastas.count(f) + rev_comp_fastas.count(f) == 1:
            error_reads.append(f)

    all_reads: Set[Fasta] = set(fastas + rev_comp_fastas) - set(error_reads)
    ret: List[str] = []
    for read in error_reads:
        hamming_1 = [hm for hm in all_reads if hamming(hm.sequence, read.sequence) == 1][0]
        ret.append(f'{read.sequence}->{hamming_1.sequence}')

    return '\n'.join(ret)


def concensus_and_profile(s: str) -> str:
    """
    Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.

    Return: A consensus string and profile matrix for the collection. (If several possible consensus strings exist,
    then you may return any one of them.)
    """
    fastas: List[Fasta] = Fasta.from_str(s)
    concensus, profile = get_consensus(fastas)
    profile_str: str = f'A: {" ".join([str(x) for x in profile["A"]])}\nC: {" ".join([str(x) for x in profile["C"]])}\nG: {" ".join([str(x) for x in profile["G"]])}\nT: {" ".join([str(x) for x in profile["T"]])}'
    return f'{concensus[0].sequence}\n{profile_str}'


def creating_a_distance_matrix(s: str) -> str:
    """
    Given: A collection of n (n≤10) DNA strings s1,…,sn of equal length (at most 1 kbp). Strings are given in FASTA
    format.

    Return: The matrix D corresponding to the p-distance dp on the given strings. As always, note that your answer is
    allowed an absolute error of 0.001.
    """
    fastas: List[Fasta] = Fasta.from_str(s)
    distances: List[str] = []
    for f in fastas:
        f_dists: str = " ".join(["{0:.5f}".format(pdistance(f, fasta)) for fasta in fastas])
        distances.append(f_dists)
        
    return "\n".join(distances)


def finding_a_shared_motif(s: str) -> str:
    """
    Given: A collection of k (k≤100) DNA strings of length at most 1 kbp each in FASTA format.

    Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single
    solution.)
    """

    def find_motif(fastas: List[Fasta], seq_len: int) -> Optional[str]:
        for f in fastas:
            motifs: Set[str] = {f[i:i+seq_len] for i in range(len(f)-seq_len + 1) if all([f[i:i+seq_len] in fasta for fasta in fastas])}
            if len(motifs) > 0:
                for first in motifs: return first  # Hacky way of getting one elem from set
        return

    def bin_motif_search(fastas: List[Fasta], lengths: List[int], last_found_motif: str) -> str:
        if len(lengths) == 1:
            new: str
            if new := find_motif(fastas, lengths[0]):
                return new
            else:
                return last_found_motif
        mid_length_index: int = len(lengths) // 2
        mid_length_value: int = lengths[mid_length_index]
        motif: Optional[str] = find_motif(fastas, mid_length_value)
        if not motif:
            return bin_motif_search(fastas, lengths[:mid_length_index], last_found_motif)
        else:
            return bin_motif_search(fastas, lengths[mid_length_index:], motif)

    fastas: List[Fasta] = Fasta.from_str(s)
    seq_len: int = max([len(f) for f in fastas])

    return bin_motif_search(fastas, list(range(seq_len)), '')


if __name__ == '__main__':
    pass

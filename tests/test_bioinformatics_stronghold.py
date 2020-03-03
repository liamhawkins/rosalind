import unittest
from typing import Set

import pytest

from bioinformatics_stronghold import counting_DNA_nucleotides, transcribing_DNA_into_RNA, \
    complementing_a_strand_of_DNA, counting_point_mutations, computing_gc_content, finding_a_motif_in_DNA, \
    translating_RNA_into_protein, calculating_protein_mass, mendels_first_law, inferring_mRNA_from_protein, \
    rabbits_and_recurrance_relations, open_reading_frames, RNA_splicing, finding_a_spliced_motif, \
    transition_transversion_ratio, calculate_expected_offspring, overlap_graphs, genome_assembly_as_shortest_superstring


def test_counting_DNA_nucleotides() -> None:
    in_: str = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
    out: str = '20 12 17 21'
    assert counting_DNA_nucleotides(in_) == out


def test_transcribing_DNA_into_RNA() -> None:
    in_: str = 'GATGGAACTTGACTACGTAAATT'
    out: str = 'GAUGGAACUUGACUACGUAAAUU'
    assert transcribing_DNA_into_RNA(in_) == out


def test_complementing_a_strand_of_DNA() -> None:
    in_: str = 'AAAACCCGGT'
    out: str = 'ACCGGGTTTT'
    assert complementing_a_strand_of_DNA(in_) == out


def test_counting_point_mutations() -> None:
    in_: str = """GAGCCTACTAACGGGAT\nCATCGTAATGACGGCCT"""
    out: int = 7
    assert counting_point_mutations(in_) == out


def test_computing_gc_content() -> None:
    in_: str = """>Rosalind_6404
    CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
    TCCCACTAATAATTCTGAGG
    >Rosalind_5959
    CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
    ATATCCATTTGTCAGCAGACACGC
    >Rosalind_0808
    CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
    TGGGAACCTGCGGGCAGTAGGTGGAAT"""
    out: str = """Rosalind_0808\n60.919540"""
    assert computing_gc_content(in_) == out


def test_finding_a_motif_in_DNA() -> None:
    in_: str = 'GATATATGCATATACTT\nATAT'
    out: str = '2 4 10'
    assert finding_a_motif_in_DNA(in_) == out


def test_translating_RNA_into_protein() -> None:
    in_: str = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
    out: str = 'MAMAPRTEINSTRING'
    assert translating_RNA_into_protein(in_) == out


def test_calculating_protein_mass() -> None:
    in_: str = 'SKADYEK'
    out: float = 821.392
    assert calculating_protein_mass(in_) == pytest.approx(out)


def test_mendels_first_law() -> None:
    in_: str = '2 2 2'
    out: float = 0.78333
    assert mendels_first_law(in_) == pytest.approx(out, rel=0.0001)


def test_inferring_mRNA_from_protein() -> None:
    in_: str = 'MA'
    out: int = 12
    assert inferring_mRNA_from_protein(in_) == out


def test_rabbits_and_recurrance_relations() -> None:
    in_: str = '5 3'
    out: int = 19
    assert rabbits_and_recurrance_relations(in_) == out


def test_open_reading_frames() -> None:
    in_: str = """>Rosalind_99
AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"""
    out: Set[str] = {'MLLGSFRLIPKETLIQVAGSSPCNLS', 'M', 'MGMTPRLGLESLLE', 'MTPRLGLESLLE'}
    assert open_reading_frames(in_) == out


def test_RNA_splicing() -> None:
    in_: str = """>Rosalind_10
ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG
>Rosalind_12
ATCGGTCGAA
>Rosalind_15
ATCGGTCGAGCGTGT"""
    out: str = 'MVYIADKQHVASREAYGHMFKVCA'
    assert RNA_splicing(in_) == out


def test_finding_a_spliced_motif() -> None:
    in_: str = """>Rosalind_14
ACGTACGTGACG
>Rosalind_18
GTA"""
    out: str = '3 8 10'
    assert finding_a_spliced_motif(in_) == out


def test_transition_transversion_ratio() -> None:
    in_: str = """>Rosalind_0209
GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGA
AGTACGGGCATCAACCCAGTT
>Rosalind_2200
TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGC
GGTACGAGTGTTCCTTTGGGT"""
    out: float = 1.21428571429
    assert transition_transversion_ratio(in_) == pytest.approx(out)


def test_calculate_expected_offspring() -> None:
    in_: str = '1 0 0 1 0 1'
    out: float = 3.5
    assert calculate_expected_offspring(in_) == out


def test_overlap_graphs() -> None:
    in_: str = """>Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATCCC
>Rosalind_5013
GGGTGGG"""
    out: str = """Rosalind_0498 Rosalind_2391
Rosalind_0498 Rosalind_0442
Rosalind_2391 Rosalind_2323"""
    assert overlap_graphs(in_) == out


def test_genome_assembly_as_shortest_superstring() -> None:
    in_: str = """>Rosalind_56
ATTAGACCTG
>Rosalind_57
CCTGCCGGAA
>Rosalind_58
AGACCTGCCG
>Rosalind_59
GCCGGAATAC"""
    out: str = 'ATTAGACCTGCCGGAATAC'
    assert genome_assembly_as_shortest_superstring(in_) == out


if __name__ == '__main__':
    unittest.main()

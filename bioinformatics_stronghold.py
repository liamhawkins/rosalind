def count_bases(s):
    """
    Given: A DNA string s of length at most 1000 nt.
    Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G',
    and 'T' occur in s.
    """
    return s.count('A'), s.count('C'), s.count('G'), s.count('T')


def transcribe(t):
    """
    Given: A DNA string t having length at most 1000 nt.
    Return: The transcribed RNA string of t.
    """
    return t.replace('T', 'U')


def reverse_complement(s):
    """
    Given: A DNA string s of length at most 1000 bp.
    Return: The reverse complement s^c of s.
    """
    comp = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
    }
    return ''.join([comp[base] for base in reversed(s)])


def num_point_mutations(s, t):
    """
    Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).
    Return: The Hamming distance dH(s,t).
    """
    num = 0
    for x, y in zip(s, t):
        num += x != y
    return num


if __name__ == '__main__':
    pass

def count_bases(s):
    return s.count('A'), s.count('C'), s.count('G'), s.count('T')


def translate(s):
    return s.replace('T', 'U')


def reverse_complement(s):
    comp = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
    }
    return ''.join([comp[base] for base in reversed(s)])


def num_point_mutations(s, t):
    num = 0
    for x, y in zip(s, t):
        num += x != y
    return num


if __name__ == '__main__':
    pass
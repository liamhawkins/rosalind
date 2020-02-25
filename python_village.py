from collections import Counter
from functools import reduce


def installing_python():
    import this


def variable_and_some_arithmetic(a, b):
    return a ** 2 + b ** 2


def strings_and_lists(s, a, b, c, d):
    return f'{s[a:b + 1]} {s[c:d + 1]}'


def conditions_and_loops(a, b):
    odds = [x for x in range(a, b + 1) if x % 2 != 0]
    return reduce(lambda x, y: x+y, odds)


def working_with_files(in_path, out_path):
    with open(in_path, 'r') as f:
        lines = f.readlines()

    with open(out_path, 'w') as f:
        for line in lines[1::2]:
            f.write(line)


def dictionaries(s):
    for key, value in dict(Counter(s.split())).items():
        print(f'{key} {value}')


if __name__ == '__main__':
    pass

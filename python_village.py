from collections import Counter
from functools import reduce


def installing_python():
    """
    After downloading and installing Python, type import this into the Python command line and see what happens.
    Then, click the "Download dataset" button below and copy the Zen of Python into the space provided
    """
    import this


def variable_and_some_arithmetic(a: int, b: int):
    """
    Given: Two positive integers a and b, each less than 1000.
    Return: The integer corresponding to the square of the hypotenuse of the right triangle whose legs have lengths
    a and b.
    """
    return a ** 2 + b ** 2


def strings_and_lists(s: str, a: int, b: int, c: int, d: int):
    """
    Given: A string s of length at most 200 letters and four integers a, b, c and d.
    Return: The slice of this string from indices a through b and c through d (with space in between), inclusively.
    In other words, we should include elements s[b] and s[d] in our slice.
    """
    return f'{s[a:b + 1]} {s[c:d + 1]}'


def conditions_and_loops(a: int, b: int):
    """
    Given: Two positive integers a and b (a<b<10000).
    Return: The sum of all odd integers from a through b, inclusively.
    """
    odds = [x for x in range(a, b + 1) if x % 2 != 0]
    return reduce(lambda x, y: x+y, odds)


def working_with_files(in_path: str, out_path: str):
    """
    Given: A file containing at most 1000 lines.
    Return: A file containing all the even-numbered lines from the original file. Assume 1-based numbering of lines.
    """
    with open(in_path, 'r') as f:
        lines = f.readlines()

    with open(out_path, 'w') as f:
        for line in lines[1::2]:
            f.write(line)


def dictionaries(s: str):
    """
    Given: A string s of length at most 10000 letters.
    Return: The number of occurrences of each word in s, where words are separated by spaces. Words are case-sensitive,
    and the lines in the output can be in any order.
    """
    for key, value in dict(Counter(s.split())).items():
        print(f'{key} {value}')


if __name__ == '__main__':
    pass

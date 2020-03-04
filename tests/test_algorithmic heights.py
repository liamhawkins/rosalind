import unittest

from algorithmic_heights import fibonacci_numbers, degree_array, binary_search, insertion_sort, double_degree_array


def test_fibonacci_numbers() -> None:
    in_: int = 6
    out: int = 8
    assert fibonacci_numbers(in_) == out


def test_degree_array() -> None:
    in_: str = """6 7
1 2
2 3
6 3
5 6
2 5
2 4
4 1"""
    out: str = "2 4 2 2 2 2"
    assert degree_array(in_) == out


def test_binary_search() -> None:
    in_: str = """5
6
10 20 30 40 50
40 10 35 15 40 20"""
    out: str = '4 1 -1 -1 4 2'
    assert binary_search(in_) == out


def test_insertion_sort() -> None:
    in_: str = """6
6 10 4 5 1 2"""
    out: int = 12
    assert insertion_sort(in_) == out


def test_double_degree_array() -> None:
    in_: str = """5 4
1 2
2 3
4 3
2 4"""
    out: str = '3 5 5 5 0'
    assert double_degree_array(in_) == out


if __name__ == '__main__':
    unittest.main()
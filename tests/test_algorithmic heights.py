import unittest

from algorithmic_heights import fibonacci_numbers, degree_array


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


if __name__ == '__main__':
    unittest.main()
import unittest

from algorithmic_heights import fibonacci_numbers


def test_fibonacci_numbers() -> None:
    in_: int = 6
    out: int = 8
    assert fibonacci_numbers(in_) == out


if __name__ == '__main__':
    unittest.main()
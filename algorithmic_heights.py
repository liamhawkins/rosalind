from typing import List

from tools.Graph import Graph
from tools.functions import bin_search, ins_sort


def fibonacci_numbers(n: int) -> int:
    """
    Given: A positive integer n≤25.

    Return: The value of Fn.
    """
    if n == 0:
        return 0
    elif n == 1:
        return 1
    return fibonacci_numbers(n-1) + fibonacci_numbers(n-2)


def degree_array(s: str) -> str:
    """
    Given: A simple graph with n≤103 vertices in the edge list format.

    Return: An array D[1..n] where D[i] is the degree of vertex i.
    """
    g = Graph.from_str(s)

    # Iterate over nodes, and join node degree into output string
    degrees: List[str] = []
    for i in range(1, len(g.nodes)+1):
        degrees.append(str(g.get_node_by_obj(i).degree))
    return ' '.join(degrees)


def binary_search(s: str) -> str:
    """
    Given: Two positive integers n≤105 and m≤105, a sorted array A[1..n] of integers from −105 to 105 and a list of m
    integers −105≤k1,k2,…,km≤105.

    Return: For each ki, output an index 1≤j≤n s.t. A[j]=ki or "-1" if there is no such index.
    """
    inp: List[str] = s.split('\n')
    sorted_array: List[int] = [int(x) for x in inp[2].split()]
    m: List[int] = [int(x) for x in inp[3].split()]

    ret = []
    for num in m:
        ret.append(str(bin_search(num, sorted_array)))

    return ' '.join(ret)


def insertion_sort(s: str) -> int:
    """
    Given: A positive integer n≤103 and an array A[1..n] of integers.

    Return: The number of swaps performed by insertion sort algorithm on A[1..n].
    """
    array: List[int] = [int(x) for x in s.split('\n')[1].split()]
    swap_count: int
    _, swap_count = ins_sort(array)
    return swap_count


def double_degree_array(s: str) -> str:
    """
    Given: A simple graph with n≤103 vertices in the edge list format.

    Return: An array D[1..n] where D[i] is the sum of the degrees of i's neighbors.
    """
    g = Graph.from_str(s)

    # Iterate over nodes, and join node degree into output string
    double_degrees: List[str] = []
    for i in range(1, len(g.nodes) + 1):
        node = g.get_node_by_obj(i)
        neighbour_degree_sum = sum([n.degree for n in node.neighbors])
        double_degrees.append(str(neighbour_degree_sum))

    return ' '.join(double_degrees)


def majority_element(s: str) -> str:
    """
    Given: A positive integer k≤20, a positive integer n≤104, and k arrays of size n containing positive integers not
    exceeding 105.

    Return: For each array, output an element of this array occurring strictly more than n/2 times if such element
    exists, and "-1" otherwise.
    """
    inp: List[str] = s.strip().split('\n')
    n: int = [int(x) for x in inp[0].split()][1]
    arrays: List[str] = inp[1:]
    ret: List[str] = []
    for array in arrays:
        ints: List[int] = [int(x) for x in array.split()]
        maj_elem: List[str] = [str(x) for x in set(ints) if ints.count(x) > n/2]
        if maj_elem:
            ret.append(maj_elem[0])
        else:
            ret.append('-1')
    return ' '.join(ret)


if __name__ == '__main__':
    pass

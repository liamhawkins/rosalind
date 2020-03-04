from typing import List

from tools.Graph import Graph


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
    edge_list: List[tuple] = [(int(a), int(b)) for a, b in [c.split() for c in [j for j in s.split('\n')]]]
    g = Graph()
    # Construct graph
    for e in edge_list[1:]:  # First line is num of nodes and edges
        node1, node2 = g.add_nodes(e)
        g.add_edge(node1, node2)

    # Iterate over nodes, and join node degree into output string
    degrees: List[str] = []
    for i in range(1, len(g.nodes)+1):
        degrees.append(str(g.get_node_by_obj(i).degree))
    return ' '.join(degrees)


if __name__ == '__main__':
    pass
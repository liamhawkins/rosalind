from itertools import permutations
from typing import List, Any


class Node:
    def __init__(self, obj: Any) -> None:
        self.obj: Any = obj
        self.edges: List[Edge] = []
        self.incoming_nodes: List[Node] = []
        self.outgoing_nodes: List[Node] = []

    def __str__(self) -> str:
        return str(self.obj)

    def __repr__(self) -> str:
        return f'Node({repr(self.obj)})'

    @property
    def degree(self) -> int:
        return len(self.edges)

    @property
    def indegree(self) -> int:
        return len(self.incoming_nodes)

    @property
    def outdegree(self) -> int:
        return len(self.outgoing_nodes)

    def is_connected(self) -> bool:
        return self.degree > 0

    def add_edge(self, edge: 'Edge') -> None:
        self.edges.append(edge)

    def add_outgoing(self, node: 'Node') -> None:
        self.outgoing_nodes.append(node)

    def add_incoming(self, node: 'Node') -> None:
        self.incoming_nodes.append(node)

    def is_leaf(self):
        return self.degree == 1

    def has_incoming(self):
        return len(self.incoming_nodes) > 0

    def has_outgoing(self):
        return len(self.outgoing_nodes) > 0


class Edge:
    def __init__(self, node1: Node, node2: Node, directed: bool = False):
        node1.add_edge(self)
        node2.add_edge(self)

        node1.add_outgoing(node2)
        node2.add_incoming(node1)

        if not directed:
            node1.add_incoming(node2)
            node2.add_outgoing(node1)

        self.node1: Node = node1
        self.node2: Node = node2
        self.directed: bool = directed

    def __str__(self) -> str:
        if self.directed:
            return f'{self.node1}-->{self.node2}'
        else:
            return f'{self.node1}<--->{self.node2}'

    def __repr__(self) -> str:
        return f'Edge({repr(self.node1)}, {repr(self.node2)}, {self.directed})'


class Graph:
    def __init__(self, directed: bool = False) -> None:
        self.nodes: List[Node] = []
        self.edges: List[Edge] = []
        self.directed: bool = directed

    def add_node(self, obj: Any) -> Node:
        node: Node = Node(obj)
        self.nodes.append(node)
        return node

    def add_nodes(self, objs: List[Any]) -> List[Node]:
        for o in objs:
            self.add_node(o)
        return self.nodes

    def add_edge(self, node1: Node, node2: Node) -> Edge:
        edge: Edge = Edge(node1, node2, directed=self.directed)
        self.edges.append(edge)
        return edge

    def join_nodes_by_func(self, f):
        for n1, n2 in permutations(self.nodes, r=2):
            if f(n1.obj, n2.obj):
                self.add_edge(n1, n2)

    def leaf_nodes(self) -> List[Node]:
        return [n for n in self.nodes if n.is_leaf()]

    def starting_nodes(self) -> List[Node]:
        if not self.directed:
            raise(ValueError('Starting nodes cannot exist on undirected graph'))
        return [n for n in self.nodes if not n.has_incoming() and n.is_connected()]

    def ending_nodes(self) -> List[Node]:
        if not self.directed:
            raise(ValueError('Ending nodes cannot exist on undirected graph'))
        return [n for n in self.nodes if not n.has_outgoing() and n.is_connected()]

    def is_linear(self) -> bool:
        return all([x.degree <= 2 for x in self.nodes]) and len(self.starting_nodes()) == 1 and len(self.ending_nodes()) == 1


if __name__ == '__main__':
    def is_double(a, b):
        return a * 2 == b

    g = Graph(directed=True)
    n1 = g.add_node(4)
    n2 = g.add_node(2)
    n3 = g.add_node(8)
    n4 = g.add_node(16)

    g.join_nodes_by_func(is_double)

    print(g.edges)
    if g.is_linear():
        order = [g.starting_nodes()[0]]

        while len(order) < len(g.edges) + 1:
            order.append(order[-1].outgoing_nodes[0])

        print(order)




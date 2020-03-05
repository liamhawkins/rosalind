from itertools import permutations
from typing import List, Any, Union, Tuple, Optional, Set


class NodeNotInGraphError(Exception):
    pass


class Node:
    def __init__(self, obj: Any) -> None:
        self.obj: Any = obj
        self.edges: List[Edge] = []
        self.incoming_nodes: List[Node] = []
        self.outgoing_nodes: List[Node] = []
        self.neighbors: List[Node] = []

    def __str__(self) -> str:
        return str(self.obj)

    def __repr__(self) -> str:
        return f'Node({repr(self.obj)})'

    def __hash__(self):
        return hash(self.obj)

    def __contains__(self, item):
        return item == self.obj

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
        if node not in self.outgoing_nodes:
            self.outgoing_nodes.append(node)
        if node not in self.neighbors:
            self.neighbors.append(node)

    def add_incoming(self, node: 'Node') -> None:
        if node not in self.incoming_nodes:
            self.incoming_nodes.append(node)
        if node not in self.neighbors:
            self.neighbors.append(node)

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

    def __contains__(self, item: Any) -> bool:
        if isinstance(item, Node):
            return item.obj in [n.obj for n in self.nodes]
        else:
            return item in [n.obj for n in self.nodes]

    def add_node(self, obj: Any) -> Node:
        if obj in self:
            print(f'{repr(obj)} already in graph')
            return self.get_node_by_obj(obj)
        node: Node = Node(obj)
        self.nodes.append(node)
        return node

    def add_nodes(self, objs: Union[List[Any], Tuple[Any, ...]]) -> List[Node]:
        nodes: List[Node] = []
        for o in objs:
            node = self.add_node(o)
            nodes.append(node)
        return nodes

    def add_edge(self, node1: Node, node2: Node) -> Edge:
        if node1 not in self or node2 not in self:
            raise NodeNotInGraphError("Can't create edge between one or more nodes not in the graph")
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

    def get_node_by_obj(self, obj: Any) -> Node:
        for n in self.nodes:
            if n.obj == obj:
                return n
        raise ValueError(f'{repr(obj)} not in graph')

    @classmethod
    def from_str(cls, s: str) -> 'Graph':
        """
        Construct graph from string in edge list format.
        e.g
            5 4
            1 2
            2 3
            4 3
            2 4
        """
        if not s:
            raise ValueError('Graph cannot be made from empty string')
        inp: List[str] = [j for j in s.split('\n')]
        num_nodes = int(inp[0].split()[0])
        edge_list = [(int(a), int(b)) for a, b in [c.split() for c in inp[1:]]]
        g: Graph = Graph()

        # Create nodes according to num_nodes
        g.add_nodes(list(range(1, num_nodes + 1)))

        # Create edges according to edge_list
        for e in edge_list:  # First line is num of nodes and edges
            node1 = g.get_node_by_obj(e[0])
            node2 = g.get_node_by_obj(e[1])
            g.add_edge(node1, node2)
        return g

    def disconnected_subgraphs(self) -> List[Set[Node]]:
        fills: List[Set[Node]] = []
        for node in self.nodes:
            f = self.fill_from_node(node)
            if f not in fills:
                fills.append(f)

        return fills

    def fill_from_node(self, node: Node, prev: Optional[Node] = None, filled: Optional[Set[Node]] = None) -> Set[Node]:
        if filled and node in filled:
            return filled

        if not filled:
            filled = {node}
        else:
            filled.add(node)

        next_nodes: Set[Node]
        if prev:
            next_nodes = set(node.neighbors) - {prev}
        else:
            next_nodes = set(node.neighbors)

        if len(next_nodes) == 0:
            return filled

        for neigh in next_nodes:
            filled.union(self.fill_from_node(neigh, node, filled))
        return filled


if __name__ == '__main__':
    g = Graph()
    n1 = g.add_node(1)
    n2 = g.add_node(2)
    n3 = g.add_node(3)
    n4 = g.add_node(4)

    g.add_edge(n1, n2)
    g.add_edge(n2, n3)
    g.add_edge(n2, n4)

    n5 = g.add_node(5)
    n6 = g.add_node(6)
    g.add_edge(n5, n6)

    n7 = g.add_node(7)

    print(g.disconnected_subgraphs())
    print(len(g.disconnected_subgraphs()))





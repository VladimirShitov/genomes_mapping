from Node import Node


class Graph:  # TODO: document

    def __init__(self):
        self.nodes = set()
        self.components = 0

    def add_node(self, data, type_):
        self.nodes.add(Node(data, type_))

    def __getitem__(self, item):
        for node in self.nodes:
            if node.data == item:
                return node
        return None






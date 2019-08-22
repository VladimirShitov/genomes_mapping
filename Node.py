

class Node:  # TODO: document

    def __init__(self, data, type_):
        self.data = data
        self.neighbors = set()
        self.type = type_

    def add_neighbors(self, neighbors):

        # To prevent possible iteration by string
        if isinstance(neighbors, str):
            neighbors = (neighbors, )

        for neighbor in neighbors:
            self.neighbors.add(neighbor)

    def __repr__(self):
        return self.data + '\n' + 'neighbours: ' + len(self.neighbors)



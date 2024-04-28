import numpy as np
import random
import matplotlib.pyplot as plt

from abc import ABC, abstractmethod

class RandomWalker(ABC):

    def __init__(self, node_count, grid_size):

        if grid_size > 100:
            raise ValueError("grid_size must be a number < 100")

        self.node_count = node_count
        self.grid_size = grid_size
        self.walk = self.create_walk()

    def __str__(self):
        return f"Self Avoiding Random Walk:\n\tnode_count = {self.node_count}\n\tgrid_size = {self.grid_size}"

    @abstractmethod
    def create_walk(self):
        pass

    def get_coords(self, walk):
        x = [pos[0] for pos in walk]
        y = [pos[1] for pos in walk]
        z = [pos[2] for pos in walk]
        return x, y, z

    def plot(self):
        x, y, z = self.get_coords(self.walk)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot(x, y, z)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()

    def get_grid_volume(self):
        return pow(self.grid_size * 2 + 1, 3)

    def get_ratio(self, a, b):
        return round(a / b * 100, 2)

    def get_walk_completeness(self):
        return self.get_ratio(len(self.walk), self.node_count)

    def show_info(self):
        print(f"should be filled in {self.get_ratio(self.node_count, self.get_grid_volume())}%")
        print(f"has been filled in {self.get_ratio(len(self.walk), self.get_grid_volume())}%")
        print(f"the walk is {self.get_walk_completeness()}% complete")

class SARW(RandomWalker):

    def __init__(self, node_count, grid_size):
        super().__init__(node_count, grid_size)

    def create_walk(self):

        MOVES = np.array([(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)])

        visited = set()
        pos = (0, 0, 0)
        walk = [pos]

        for i in range(self.node_count):
            move = random.choice(MOVES)
            new_pos = tuple(x + y for x, y in zip(pos, move))

            if any(abs(coord) > self.grid_size for coord in new_pos):
                continue

            if new_pos in visited:
                continue
            visited.add(new_pos)
            pos = new_pos
            walk.append(pos)

        return walk

class AdvancedSARW(RandomWalker):

    def __init__(self, node_count, grid_size):
        super().__init__(node_count, grid_size)

    def create_walk(self):
        
        shape = np.repeat(self.grid_size*2+1, 3)
        

class Grid():

    def __init__(self, size):
        self.size = size
        self.filled_fields = 0

        n = size*2+1
        self.visited = np.full((n, n, n), False, dtype=bool)

    def is_valid(self, field):
        if field.x >= abs(self.size):
            return False
        if field.y >= abs(self.size):
            return False
        if field.z >= abs(self.size):
            return False
        return True
    
    def is_empty(self, field):
        return not self.visited[field.x, field.y, field.z]

    def fill_field(self, field):
        self.visited[field.x, field.y, field.z] = True
        self.filled_fields += 1

    def empty_field(self, field):
        self.visited[field.x, field.y, field.z] = False
        self.filled_fields -= 1

    

class Field():

    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def get_neighbours_at_distance(self, distance):
        #   TODO
        #   implement
        pass
import numpy as np
import random
import matplotlib.pyplot as plt

from abc import ABC, abstractmethod

class RandomWalker(ABC):

    MOVES = np.array([(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)])

    def __init__(self, node_count, grid_size):
        self.node_count = node_count
        self.grid_size = grid_size
        self.path = self.create_path()

    def __str__(self):
        return f"Self Avoiding Random Walk:\n\tnode_count = {self.node_count}\n\tgrid_size = {self.grid_size}"

    @abstractmethod
    def create_path(self):
        pass

    def get_coords(self, path):
        x = [pos[0] for pos in path]
        y = [pos[1] for pos in path]
        z = [pos[2] for pos in path]
        return x, y, z

    def plot_path(self):
        x, y, z = self.get_coords(self.path)
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

    def path_completeness_ratio(self):
        return self.get_ratio(len(self.path), self.node_count)

    def show_info(self):
        print(f"should be filled in {self.get_ratio(self.node_count, self.get_grid_volume())}%")
        print(f"has been filled in {self.get_ratio(len(self.path), self.get_grid_volume())}%")
        print(f"the path is {self.path_completeness_ratio()}% complete")

    
class SARW(RandomWalker):

    def __init__(self, node_count, grid_size):
        super().__init__(node_count, grid_size)

    def create_path(self):
        print("creating a path for simple walker")
        visited = set()
        pos = (0, 0, 0)
        path = [pos]

        for i in range(self.node_count):
            move = random.choice(self.MOVES)
            new_pos = tuple(x + y for x, y in zip(pos, move))

            if any(abs(coord) > self.grid_size for coord in new_pos):
                continue

            if new_pos in visited:
                continue
            visited.add(new_pos)
            pos = new_pos
            path.append(pos)

        return path
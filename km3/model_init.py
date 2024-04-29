import numpy as np
import random
import matplotlib.pyplot as plt

from abc import ABC, abstractmethod

class RandomWalker(ABC):

    def __init__(self, node_count, grid_size):

        self.node_count = node_count
        self.grid = Grid(grid_size)
        self.walk = self.create_walk()

    def __str__(self):
        # for field in self.walk:
        #     print(field)
        return f"Self Avoiding Random Walk:\n\tnode_count = {self.node_count}\n\tgrid_size = {self.grid.size}"

    @abstractmethod
    def create_walk(self):
        pass

    def get_coords(self, walk):
        x = [pos[0] for pos in walk]
        y = [pos[1] for pos in walk]
        z = [pos[2] for pos in walk]
        return x, y, z

    def get_coords2(self, walk):
        x = [field.x for field in walk]
        y = [field.y for field in walk]
        z = [field.z for field in walk]
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

    def plot2(self):
        x, y, z = self.get_coords2(self.walk)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot(x, y, z)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()

    def get_ratio(self, a, b):
        return round(a / b * 100, 2)

    def get_walk_completeness(self):
        return self.get_ratio(len(self.walk), self.node_count)

    def show_info(self):
        print(f"should be filled in {self.get_ratio(self.node_count, self.grid.volume)}%")
        print(f"has been filled in {self.get_ratio(len(self.walk), self.grid.volume)}%")
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

            if any(abs(coord) > self.grid.size for coord in new_pos):
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
        
        if self.grid.volume < self.node_count:
            raise ValueError(f"Invalid parameters for walk object:\n{self.__str__()}\nGrid is too small to contain this walk. Consider changing grid_size or node_count values.")

        walk = np.empty(self.node_count, dtype=Field)
        current_field = Field()
        for i in range(self.node_count):
            walk[i] = current_field
            current_field = self.get_child(current_field)

        return walk
        
    def get_child(self, field):

        potential_children = np.empty(0, dtype=Field)
        distance = 0
        while len(potential_children) == 0:
            distance += 1
            potential_children = field.get_neighbours_at_distance(distance)
            potential_children = self.grid.validate(potential_children)
            potential_children = self.grid.check_availability(potential_children)

        chosen_child = np.random.choice(potential_children)
        return chosen_child

class Grid():

    def __init__(self, size):
        if size > 100:
            raise ValueError("grid_size must be a number < 100")
        self.size = size
        self._filled_fields = 0

        n = size*2+1
        self.visited = np.full((n, n, n), False, dtype=bool)

    def is_valid(self, field):
        if abs(field.x) >= self.size:
            return False
        if abs(field.y) >= self.size:
            return False
        if abs(field.z) >= self.size:
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

    def validate(self, fields):
        mask = np.array([self.is_valid(field) for field in fields])
        return fields[mask]

    def check_availability(self, fields):
        mask = np.array([self.is_empty(field) for field in fields])
        return fields[mask]
    
    @property
    def volume(self):
        return pow(self.size*2+1, 3)

    @property
    def filled_fields(self):
        return self._filled_fields

    @property
    def empty_fields(self):
        return self.volume - self.filled_fields

class Field():

    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z
    
    def __eq__(self, other):
        if isinstance(other, Field):
            return self.x == other.x and self.y == other.y and self.z == other.z
        return False

    def __lt__(self, other):
        if isinstance(other, Field):
            return (self.x < other.x) or (self.x == other.x and self.y < other.y) or (self.x == other.x and self.y == other.y and self.z < other.z)
        return NotImplemented

    def __le__(self, other):
        if isinstance(other, Field):
            return self == other or self < other
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, Field):
            return not (self <= other)
        return NotImplemented

    def __ge__(self, other):
        if isinstance(other, Field):
            return not (self < other)
        return NotImplemented

    def __str__(self):
        return f"Field: ({self.x}, {self.y}, {self.z})"

    def get_neighbours_at_distance(self, distance):

        neighbours = np.empty(0, dtype=Field)
        
        for dx in range(0, distance+1):
            for dy in range(0, distance+1-dx):
                dz = distance-(dx+dy)

                new_neighbours = np.array([
                    Field(self.x+dx, self.y+dy, self.z+dz),
                    Field(self.x+dx, self.y+dy, self.z-dz),
                    Field(self.x+dx, self.y-dy, self.z+dz),
                    Field(self.x+dx, self.y-dy, self.z-dz),
                    Field(self.x-dx, self.y+dy, self.z+dz),
                    Field(self.x-dx, self.y+dy, self.z-dz),
                    Field(self.x-dx, self.y-dy, self.z+dz),
                    Field(self.x-dx, self.y-dy, self.z-dz)
                ], dtype=Field)
                new_neighbours = np.unique(new_neighbours)
                neighbours = np.append(neighbours, new_neighbours)
        
        return neighbours
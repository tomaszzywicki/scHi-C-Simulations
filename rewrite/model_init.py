import numpy as np
import matplotlib.pyplot as plt

from abc import ABC, abstractmethod

class RandomWalk(ABC):

    def __init__(self, node_count, grid_size):
        self.node_count = node_count
        self.grid_size = grid_size
        self._filled_fields = 0

        self.create_walk()

    def __str__(self):
        return f"Self Avoiding Random Walk:\n\tnode_count = {self.node_count}\n\tgrid_size = {self.grid_size}"
        
    @abstractmethod
    def create_walk(self):
        if self.volume < self.node_count:
            raise ValueError(f"Invalid parameters for walk object:\n{self.__str__()}\nGrid is too small to contain this walk. Consider changing grid_size or node_count values.")

    def get_coords(self, walk=None):
        if walk is None:
            walk = self.walk
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

    def get_ratio(self, a, b):
        return round(a / b * 100, 2)

    def get_walk_completeness(self):
        return self.get_ratio(len(self.walk), self.node_count)

    def show_info(self):
        print(f"should be filled in {self.get_ratio(self.node_count, self.volume)}%")
        print(f"has been filled in {self.get_ratio(len(self.walk), self.volume)}%")
        print(f"the walk is {self.get_walk_completeness()}% complete")

    def is_valid(self, field):
        if abs(field.x) >= self.grid_size:
            return False
        if abs(field.y) >= self.grid_size:
            return False
        if abs(field.z) >= self.grid_size:
            return False
        return True

    def is_empty(self, field):
        return not self.has_field(field)

    def validate(self, fields):
        mask = np.array([self.is_valid(field) for field in fields])
        return fields[mask]

    def check_availability(self, fields):
        mask = np.array([self.is_empty(field) for field in fields])
        return fields[mask]

    def has_field(self, field):
        return field in self.walk
    
    @property
    def volume(self):
        return pow(self.grid_size*2+1, 3)

    @property
    def filled_fields(self):
        return self._filled_fields

    @property
    def empty_fields(self):
        return self.volume - self.filled_fields

    
class SARW(RandomWalk):

    def __init__(self, node_count, grid_size):
        super().__init__(node_count, grid_size)

    def create_walk(self):
        # overall time complexity: n^2
        # could be nlogn, but implementation hasn't yet proved to be worth the hassle
        super().create_walk()

        self.walk = np.empty(self.node_count, dtype=Field)
        current_field = Field()
        for i in range(self.node_count):
            self.walk[i] = current_field
            current_field = self.get_child(current_field)

    def get_child(self, field):
        potential_children = np.empty(0, dtype=Field)
        distance = 0
        while len(potential_children) == 0:
            distance += 1
            potential_children = field.get_neighbours_at_distance(distance)
            potential_children = self.validate(potential_children)
            potential_children = self.check_availability(potential_children)

        chosen_child = np.random.choice(potential_children)
        return chosen_child
    
    def get_field(self, index):
        return self.walk[index]

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

    @staticmethod
    def get_distance(field1, field2):
        distance = pow(
            pow(abs(field1.x - field2.x), 2) + 
            pow(abs(field1.y - field2.y), 2) + 
            pow(abs(field1.z - field2.z), 2),
            1/2)
        return distance

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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import copy
import pickle

import model_init

class scData():
    
    def __init__(self, filepath, sep='\t', chromosome=None):
        self.df = pd.read_csv(filepath, sep=sep)
        if chromosome:
            self.isolate_chromosome(chromosome)
        
    def prep(self, resolution=1e7):
        self.resolution = resolution
        self.chromosomes = self.get_all_chromosomes()
        self.bin_count = self.chromosomes.iloc[-1]["last_bin"] + 1
        self.contact_matrix = self.generate_contact_matrix()
        self.theta_matrix = self.generate_theta_matrix()

    def generate_contact_matrix(self):
        contact_matrix = np.zeros((self.bin_count, self.bin_count), dtype=int)
        
        for _, row in self.df.iterrows():
            bin_index1 = (self.get_bin(row["chrom1"], row["coord1"]))
            bin_index2 = (self.get_bin(row["chrom2"], row["coord2"]))
            if (bin_index1 == bin_index2):
                continue

            contact_matrix[bin_index1, bin_index2] += 1
            contact_matrix[bin_index2, bin_index1] += 1

        return contact_matrix

    def plot_matrix(self, matrix):

        plt.figure(figsize=(8, 6))
        plt.imshow(matrix, cmap='Reds', interpolation='nearest')
        plt.colorbar()
        plt.xlabel('Genomic Bins')
        plt.ylabel('Genomic Bins')

        plt.show()

    def get_all_chromosomes(self):
        chromosomes1 = self.df["chrom1"].unique()
        chromosomes2 = self.df["chrom2"].unique()
        chrom_list = np.unique(np.concatenate((chromosomes1, chromosomes2), axis=0))
        sorted_chromosomes = sorted(chrom_list, key=lambda x: int(x) if x.isdigit() else float('inf'))

        records = []
        end_bin = -1
        for chromosome in sorted_chromosomes:
            min_coord, max_coord = self.get_coordinate_range(chromosome)
            bin_count = self.get_bin_count(min_coord, max_coord)
            record = {'chromosome': chromosome, 'first_bin': end_bin+1, 'last_bin': end_bin+bin_count, 'min_coord': min_coord, 'max_coord': max_coord}
            records.append(record)
            end_bin += bin_count
        
        return pd.DataFrame(records)

    def isolate_chromosome(self, chromosome):
        self.df = self.df[(self.df["chrom1"] == chromosome) & (self.df["chrom2"] == chromosome)]

    def get_coordinate_range(self, chromosome):
        df_for_chromosome = self.df[(self.df['chrom1'] == chromosome) & (self.df["chrom2"] == chromosome)]
        coords = np.concatenate([df_for_chromosome["coord1"].values, df_for_chromosome["coord2"].values])
        min_coord = min(coords)
        max_coord = max(coords)
        return (min_coord, max_coord)
    
    def get_bin_count(self, min_coord, max_coord):
        return math.ceil((max_coord - min_coord) / self.resolution)

    def get_bin(self, chromosome, coord):
        
        record = self.chromosomes[self.chromosomes["chromosome"] == chromosome].iloc[0]
        
        range = record["max_coord"] - record["min_coord"] + 1
        num_bins = record["last_bin"] - record["first_bin"] + 1
        bin_size = range/num_bins

        bin = int((coord - record["min_coord"]) // bin_size)
        return bin + record["first_bin"]

    def get_contacts(self):
        
        contacts = []
        for i in range(self.bin_count-1):
            for j in range(i+1, self.bin_count):
                if self.contact_matrix[i][j] > 0:
                    contacts.append((i, j))
        return contacts

    def generate_theta_matrix(self):
        
        d0 = self.bin_count
        theta_matrix = np.zeros((self.bin_count, self.bin_count))

        print("preping theta matrix...")
        contacts = self.get_contacts()

        for i in range(self.bin_count-1):
            print(i, " of ", self.bin_count)
            for j in range(i+1, self.bin_count):
                value = 0

                for bin_index1, bin_index2 in contacts:
                    if abs(bin_index1 - i) > d0 or abs(bin_index2 - j) > d0:
                        continue
                    value += self.zyrafa(i, j, bin_index1, bin_index2)

                theta_matrix[i][j] = min(1, value)

        return theta_matrix

    def zyrafa(self, i, j, x, y):

        mu2 = 2
        return math.exp(-(pow(x-i, 2) / mu2 + pow(y-j, 2) / mu2))
        

class Model():

    def __init__(self, data):
        self.data = data

        # dependant on the data
        # determine bin count -> will be the node count for sarw
        # determine grid_size -> dependant on i don't know what yet
        
        self.path = model_init.AdvancedSARW(self.data.bin_count, 100)

    def evolve(self, iterations=500):

        
        for i in range(iterations):

            best_evaluation = self.evaluate(self.path.walk)
            print(f"after iteration {i}: ", round(best_evaluation, 2))
            best_candidate = self.path.walk
            candidates = self.generate_sibling_walks()
            for candidate in candidates:
                evaluation = self.evaluate(candidate)
                if evaluation < best_evaluation:
                    best_evaluation = evaluation
                    best_candidate = candidate
                
            self.path.walk = best_candidate

    def generate_sibling_walks(self, count=10, index_to_modify=None):

        if not index_to_modify:
            index_to_modify = random.randint(0, self.data.bin_count-1)

        new_walks = []
        for i in range(count):
            new_walk = copy.deepcopy(self.path.walk)
            new_walk[index_to_modify].x += random.uniform(-5,5)
            new_walk[index_to_modify].y += random.uniform(-5,5)
            new_walk[index_to_modify].z += random.uniform(-5,5)
            new_walks.append(new_walk)

        return new_walks

    def delta(self, i, j):
        return self.delta0 / pow(min(1, self.data.theta_matrix[i][j]), 1/3)

    @property
    def delta1(self):
        return self.delta0 / pow(min(1, self.theta1), 1/3)

    def d(self, walk, i, j):
        
        return model_init.Field.get_distance(walk[i], walk[j])

    def evaluate(self, walk):

        self.delta0 = 8
        self.theta1 = 0.7
        self.beta = 1
        self.tau = 1
        self.mu1 = 20
        self.rho = 1
        self.phi = 0.1
        
        result = 0
        for i in range(self.data.bin_count-1):
            for j in range(i+1, self.data.bin_count):
                if self.data.contact_matrix[i][j] >= 1 or self.data.theta_matrix[i][j] == 1:
                    result += pow(self.d(walk, i, j) - self.delta0, 2) / pow(self.delta0, 2)
                elif self.theta1 < self.data.theta_matrix[i][j]:
                    result += self.beta * (1 - math.exp(-(pow(self.d(walk, i, j) - self.delta(i, j), 2) / self.mu1)))
                else:
                    result += self.tau * (1 - 1 / (1 + math.exp(-(self.d(walk, i, j) - (self.delta1 - self.rho)) / self.phi)))
        
        return result
        
    def save_to_file(self, file_path):
        with open(file_path, 'wb') as file:
            pickle.dump(self, file)

    @staticmethod
    def load_from_file(file_path):
        with open(file_path, 'rb') as file:
            return pickle.load(file)
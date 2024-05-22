import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import pickle
import math
import copy
from scipy.stats import chi2

import model_init

def load(file_path):
    with open(file_path, 'rb') as file:
        return pickle.load(file)

class scData():

    def __init__(self, filepath, sep='\t', mu2=2):
        self.base_df = pd.read_csv(filepath, sep=sep)
        self.resolution = None
        self.chromosome = None
        self.contact_matrix = None
        self.theta_matrix = None
        self.mu2 = mu2
    
    def prep(self, chromosome=None, resolution=1e7):
        self.isolate_chromosome(chromosome)
        if not chromosome:
            self.chromosome = "ALL"
        else:
            self.chromosome = chromosome
        self.resolution = int(resolution)
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

    def scHiC_matrix(self, matrix=None):
        if not matrix:
            matrix = self.contact_matrix
        plt.figure(figsize=(8, 6))
        plt.imshow(matrix, cmap='Reds', interpolation='nearest')
        plt.colorbar()
        plt.xlabel('Genomic Bins')
        plt.ylabel('Genomic Bins')
        plt.show()
    
    def get_chromosome_list(self, df=None):

        if df is None:
            df = self.base_df

        chromosomes1 = df["chrom1"].unique()
        chromosomes2 = df["chrom2"].unique()
        chrom_list = np.unique(np.concatenate((chromosomes1, chromosomes2), axis=0))
        sorted_chromosomes = sorted(chrom_list, key=lambda x: int(x) if x.isdigit() else float('inf'))
        return sorted_chromosomes

    def get_all_chromosomes(self):
        chromosomes = self.get_chromosome_list(self.df)
        records = []
        end_bin = -1
        for chromosome in chromosomes:
            min_coord, max_coord = self.get_coordinate_range(chromosome)
            bin_count = self.get_bin_count(min_coord, max_coord)
            record = {'chromosome': chromosome, 'first_bin': end_bin+1, 'last_bin': end_bin+bin_count, 'min_coord': min_coord, 'max_coord': max_coord}
            records.append(record)
            end_bin += bin_count
        
        return pd.DataFrame(records)

    def isolate_chromosome(self, chromosome):
        if not chromosome:
            self.df = self.base_df
        else:
            self.df = self.base_df[(self.base_df["chrom1"] == chromosome) & (self.base_df["chrom2"] == chromosome)]

    def get_coordinate_range(self, chromosome):
        df_for_chromosome = self.base_df[(self.base_df['chrom1'] == chromosome) & (self.base_df["chrom2"] == chromosome)]
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

        print("preparing theta matrix...")
        contacts = self.get_contacts()

        for i in range(self.bin_count-1):
            print(i, " of ", self.bin_count)
            for j in range(i+1, self.bin_count):
                value = 0
                for bin_index1, bin_index2 in contacts:
                    if abs(bin_index1 - i) > d0 or abs(bin_index2 - j) > d0:
                        continue
                    value += self.pair_score(i, j, bin_index1, bin_index2)
                theta_matrix[i][j] = min(1, value)

        return theta_matrix

    def pair_score(self, i, j, x, y):
        return math.exp(-(pow(x-i, 2) / self.mu2 + pow(y-j, 2) / self.mu2))

    def is_raw(self):
        return self.resolution is None or self.contact_matrix is None or self.theta_matrix is None

    def save(self, file_path=None):
        if not file_path:
            if self.is_raw():
                file_path = f"data/raw_data.pkl"
            else:
                file_path = f"data/chrom{self.chromosome}_res{self.resolution}.pkl"
        with open(file_path, 'wb') as file:
            pickle.dump(self, file)


class Model():

    def __init__(self, data, delta0=8, theta1=0.7, beta=1, tau=1, mu1=20, rho=1, phi=0.1):

        if data.is_raw():
            raise ValueError("Provided data has not been processed")

        self.data = data
        self.n = data.bin_count
        self.walk = model_init.SARW(self.n, 100)

        self.delta0 = delta0
        self.theta1 = theta1
        self.beta = beta
        self.tau = tau
        self.mu1 = mu1
        self.rho = rho
        self.phi = phi

        self.evaluation_score = self.evaluate(self.walk.walk)

    def score_to_prob(self, score):
        percentage = (score - self.evaluation_score) / self.evaluation_score * 100
        x = max(0, 3 + percentage)
        print(f"score: {score}\tcurrent score: {self.evaluation_score}\tpercentage: {percentage}\tx: {x}")
        return chi2.pdf(x, 2)
        

    def evolve(self, iterations=500, step=5):

        for i in range(iterations):
            print(f"after iteration {i}: {round(self.evaluation_score,2)}")
            candidate, changed_index = self.generate_sibling_walks(count=1, step=step)
            candidate_score = self.reevaluate(old_walk=self.walk.walk, new_walk=candidate, changed_index=changed_index)
            candidate_prob = self.score_to_prob(candidate_score)
            print(f"candidate prob: {round(candidate_prob,2)}")
            acceptance_prob = min(1, candidate_prob)
            if random.uniform(0, 1) < acceptance_prob:
                self.walk.walk = candidate
                self.evaluation_score = candidate_score


    def generate_sibling_walks(self, count=10, index_to_modify=None, step=5):
        if not index_to_modify:
            index_to_modify = random.randint(0, self.data.bin_count-1)

        new_walks = []
        for i in range(count):
            new_walk = copy.deepcopy(self.walk.walk)
            new_walk[index_to_modify].x += random.uniform(-step,step)
            new_walk[index_to_modify].y += random.uniform(-step,step)
            new_walk[index_to_modify].z += random.uniform(-step,step)
            new_walks.append(new_walk)

        if len(new_walks) == 1:
            return new_walks[0], index_to_modify

        return new_walks, index_to_modify
    
    def delta(self, i, j):
        return self.delta0 / pow(min(1, self.data.theta_matrix[i][j]), 1/3)

    @property
    def delta1(self):
        return self.delta0 / pow(min(1, self.theta1), 1/3)

    def d(self, walk, i, j):
        return model_init.Field.get_distance(walk[i], walk[j])

    def evaluate(self, walk):
        result = 0
        for i in range(self.n-1):
            for j in range(i+1, self.n):
                result += self.pair_score(walk, i, j)
        return result
        
    def reevaluate(self, old_walk, new_walk, changed_index):
        result = self.evaluation_score
        i = changed_index
        for j in range(self.n):
            if j == i:
                continue
            result -= self.pair_score(old_walk, i, j)
            result += self.pair_score(new_walk, i, j)
        return result

    def pair_score(self, walk, i, j):

        if i > j:
            i, j = j, i

        if self.data.contact_matrix[i][j] >= 1 or self.data.theta_matrix[i][j] == 1:
            return pow(self.d(walk, i, j) - self.delta0, 2) / pow(self.delta0, 2)
        elif self.theta1 < self.data.theta_matrix[i][j]:
            return self.beta * (1 - math.exp(-(pow(self.d(walk, i, j) - self.delta(i, j), 2) / self.mu1)))
        else:
            return self.tau * (1 - 1 / (1 + math.exp(-(self.d(walk, i, j) - (self.delta1 - self.rho)) / self.phi)))

    def plot(self):
        self.walk.plot()

    def save(self, file_path):
        with open(file_path, 'wb') as file:
            pickle.dump(self, file)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

import model_init

class scData():
    
    def __init__(self, filepath, sep='\t', resolution=1e6):
        self.df = pd.read_csv(filepath, sep=sep)
        self.resolution = resolution
        self.chromosomes = self.get_all_chromosomes()
        print(self.chromosomes)
        self.bin_count = self.chromosomes.iloc[-1]["last_bin"] + 1
        self.contact_matrix = self.generate_contact_matrix()
        self.plot_contact_matrix()

    def generate_contact_matrix(self):
        contact_matrix = np.zeros((self.bin_count, self.bin_count), dtype=int)
        
        for _, row in self.df.iterrows():
            bin_index1 = (self.get_bin(row["chrom1"], row["coord1"]))
            bin_index2 = (self.get_bin(row["chrom2"], row["coord2"]))

            contact_matrix[bin_index1, bin_index2] += 1
            contact_matrix[bin_index2, bin_index1] += 1

        return contact_matrix

    def plot_contact_matrix(self):

        plt.figure(figsize=(8, 6))
        plt.imshow(self.contact_matrix, cmap='Reds', interpolation='nearest')
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

class Model():

    def __init__(self, data):
        self.data = data

        # dependant on the data
        # determine bin count -> will be the node count for sarw
        # determine grid_size -> dependant on i don't know what yet
        
        # self.path = model_init.AdvancedSARW()


    def evolve(self, iterations=1e3):
        pass

        # generate siblings
        # evaluate siblings
        # choose the best sibling model


    def next_generation(self):
        pass

    def generate_siblings(self, count=1e2):
        pass

    def evaluate(self):
        pass
    
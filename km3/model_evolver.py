import pandas as pd
import numpy as np

import model_init

class scData():
    
    def __init__(self, filepath, sep='\t'):
        self.df = pd.read_csv(filepath, sep=sep)
        self.chromosomes = self.get_all_chromosomes()

    def get_all_chromosomes(self):
        chromosomes1 = self.df["chrom1"].unique()
        chromosomes2 = self.df["chrom2"].unique()
        chromosomes = np.unique(np.concatenate((chromosomes1, chromosomes2), axis=0))
        print(chromosomes)

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
    
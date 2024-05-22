import model_init
import model_evolver

node_count = int(1e4)
grid_size = int(1e2)
data_file_path = "data/GSM1173493_cell-1.txt"
chrom_num = "5"
resolution = int(1e6)
chrom5_data = f"data/chrom{chrom_num}_{resolution}.pkl"

def prep_data(resolution=1e7):
    data = model_evolver.scData(data_file_path)
    data.save()
    chromosomes = data.get_chromosome_list()
    for chromosome in chromosomes:
        data.prep(resolution=resolution, chromosome=chromosome)
        data.save()
    data.prep(resolution=resolution)
    data.save()

if __name__ == "__main__":
    # prep_data(resolution=1e6)
    # data = model_evolver.scData.load("data/chrom1_res1000000.pkl")
    # data.scHiC_matrix()
    # model = model_evolver.Model(data)
    model = model_evolver.Model.load("models/new.pkl")
    model.plot()
    model.evolve(iterations=1000, step=1)
    model.save("models/new.pkl")
    model.plot()
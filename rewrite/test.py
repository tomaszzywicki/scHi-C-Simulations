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
    # data = model_evolver.scData(data_file_path)
    # data.prep(chromosome="1", resolution=5e5)
    # data.save()
    # data.scHiC_matrix()
    # data = model_evolver.load("data/chrom1_res500000.pkl")
    # model = model_evolver.Model(data)
    model = model_evolver.load("models/high_res3test.pkl")
    model.plot()
    model.evolve(iterations=10000, step=0.01)
    model.save("models/high_res3test.pkl")
    model.plot()
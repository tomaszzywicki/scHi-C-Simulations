import model_init
import model_evolver

node_count = int(1e4)
grid_size = int(1e2)
data_file_path = "data/GSM1173493_cell-1.txt"
chrom_num = "5"
resolution = int(1e6)
chrom5_data = f"data/chrom{chrom_num}_{resolution}.pkl"

if __name__ == "__main__":
    # data = model_evolver.scData(data_file_path)
    # data.save()
    data = model_evolver.scData.load("data/chromALL_res10000000.pkl")
    # chromosomes = data.get_chromosome_list()
    # for chromosome in chromosomes:
    #     data.prep(resolution=1e7, chromosome=chromosome)
        # data.save()
    data.scHiC_matrix()
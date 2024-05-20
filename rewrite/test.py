import model_init
import model_evolver

node_count = int(1e4)
grid_size = int(1e2)
data_file_path = "data/GSM1173493_cell-1.txt"
chrom5_data = "data/chrom5_1e6.pkl"
chrom_num = "5"

if __name__ == "__main__":
    # data = model_evolver.scData(data_file_path)
    # data.prep(resolution=1e6, chromosome=chrom_num)
    # data.save(chrom5_data)
    data = model_evolver.scData.load(chrom5_data)
    data.scHiC_matrix()
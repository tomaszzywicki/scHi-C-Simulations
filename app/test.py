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
    # data = model_evolver.load("data/chrom1_res1000000.pkl")
    # data.scHiC_matrix()
    # model = model_evolver.Model(data)
    # model = model_evolver.load("models/2M_1e6_c1.pkl")
    # model.plot()
    # model.evolve_simulated_annealing(iterations=10000, step=5)
    # model.evolve(iterations=50000, step=5)
    # model.save("models/sim_ann_test.pkl")
    # model.plot()
    # model.show_info()   

    # ----- TESTING ----- 
    # data = model_evolver.load("data/chrom1_res1000000.pkl")
    # model = model_evolver.Model(data)
    model = model_evolver.load("models/oskar_pjastri.pkl")
    # model.evolve_simulated_annealing(iterations=10000, step=5)
    model.save("models/oskar_pjastri.pkl")
    model.plot_score_history(start_iter=3000, end_iter=4000)
    model.plot_walk_history()
    # model.plot()  
    
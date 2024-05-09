import model_evolver

scHIC_file_path = "data/GSM1173493_cell-1.txt"
model_file_path = "test_model.pkl"

def initialize_model(chromosome=None, resolution=1e6):
    data = model_evolver.scData(scHIC_file_path)
    if chromosome:
        data.isolate_chromosome(chromosome)
    data.prep(resolution=resolution)
    return model_evolver.Model(data)

def load_model(file_path="test_model.pkl"):
    return model_evolver.Model.load_from_file(file_path)


if __name__ == "__main__":
    model = load_model(file_path="model_delta_5_1500.pkl")
    model.path.plot2()
    # model.evolve(iterations=500)
    # model.path.plot2()
    # model.save_to_file("model_delta_5_2000.pkl")
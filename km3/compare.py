import model_evolver

scHIC_file_path = "data/GSM1173493_cell-1.txt"

def prepare_data(chromosome=None, resolution=1e6):
    data = model_evolver.scData(scHIC_file_path)
    if chromosome:
        data.isolate_chromosome(chromosome)
    data.prep(resolution)
    return data

def generate_models(data, model_count, evolution_steps=100):

    for i in range(model_count):
        init_path = "comparison/init_" + str(i+1) + ".pkl"
        model_path = "comparison/model_" + str(i+1) + ".pkl"

        model = model_evolver.Model(data)
        model.save_to_file(init_path)

        print("==========================")
        print("==========================")
        print("==========================")
        print(f"Training model nr {i+1}")
        print("==========================")
        print("==========================")
        print("==========================")

        model.evolve(iterations=evolution_steps)
        model.save_to_file(model_path)

def compare_models(model_count):
    for i in range(model_count):
        init_path = "comparison/init_" + str(i+1) + ".pkl"
        model = model_evolver.Model.load_from_file(init_path)
        print("==========================")
        print(f"Showing initial model nr {i+1}")
        print(f"Model evaluation: {round(model.evaluate(model.path.walk), 2)}")
        model.path.plot2()

    for i in range(model_count):
        model_path = "comparison/model_" + str(i+1) + ".pkl"
        model = model_evolver.Model.load_from_file(model_path)
        print("==========================")
        print(f"Showing trained model nr {i+1}")
        print(f"Model evaluation: {round(model.evaluate(model.path.walk), 2)}")
        model.path.plot2()

if __name__ == "__main__":

    model_count = 6
    evolution_steps = 10000

    # data = prepare_data(chromosome="1")
    # generate_models(data, model_count=model_count, evolution_steps=evolution_steps)
    compare_models(model_count=model_count)

import model_evolver

model = model_evolver.Model.load_from_file("test_model_evolved_5.pkl")
model.path.plot2()
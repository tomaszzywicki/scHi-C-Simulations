import model_evolver

model = model_evolver.Model.load_from_file("test_model.pkl")
model.path.plot2()
model.evolve(iterations=5)
model.path.plot2()
model.save_to_file("test_model_evolved_5.pkl")
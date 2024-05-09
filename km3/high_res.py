import model_evolver
import test

# model = test.initialize_model(chromosome="1", resolution=1e5)
# model.path.plot2()
# model.save_to_file("init_high_res.pkl")

model = model_evolver.Model.load_from_file("trained_high_res.pkl")
model.evolve(iterations=10)
model.path.plot2()
model.save_to_file("trained_high_res.pkl")

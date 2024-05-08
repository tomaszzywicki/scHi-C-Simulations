import model_init
import validator
import model_evolver

# validator = validator.WalkerValidator(model_init.AdvancedSARW)
# validator = validator.WalkerValidator(model_init.SARW)
# validator.run()

# walk = model_init.AdvancedSARW(int(1e4), 40)
# print(walk)
# walk.show_info()
# walk.plot2()

# walk = model_init.SARW(10000,40)
# print(walk)
# walk.show_info()
# walk.plot()

# data = model_evolver.scData('data/GSM1173493_cell-1.txt', chromosome="2")
# data.prep()

data = model_evolver.scData('data/GSM1173493_cell-1.txt')
data.isolate_chromosome("1")
data.prep(resolution=5e5)
# data.plot_matrix(data.theta_matrix)

model = model_evolver.Model(data)
model.path.plot2()
model.evolve()
print(model.evaluate(model.path.walk))
model.path.plot2()


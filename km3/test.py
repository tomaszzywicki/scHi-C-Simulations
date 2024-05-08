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
data.isolate_chromosome("5")
data.prep(resolution=1e6)
# data.plot_matrix(data.contact_matrix)
# model = model_evolver.Model(data)
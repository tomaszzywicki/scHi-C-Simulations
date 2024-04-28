import model_init
import validator

node_count = 1000
grid_size = 20

# validator = validator.WalkerValidator(model_init.AdvancedSARW)
# validator.run()

# walk = model_init.SARW(node_count=node_count, grid_size=grid_size)
# walk.show_info()
# walk.plot_path()

walk = model_init.AdvancedSARW(node_count=node_count, grid_size=grid_size)


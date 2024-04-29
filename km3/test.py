import model_init
import validator

# validator = validator.WalkerValidator(model_init.AdvancedSARW)
# validator = validator.WalkerValidator(model_init.SARW)
# validator.run()

walk = model_init.AdvancedSARW(10000, 40)
print(walk)
walk.show_info()
walk.plot2()

walk = model_init.SARW(10000,40)
print(walk)
walk.show_info()
walk.plot()
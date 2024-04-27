import model_init
import validator

validator = validator.WalkerValidator(model_init.SARW)
validator.run()
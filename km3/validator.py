class WalkerValidator:

    def __init__(self, Walker):
        self.Walker = Walker

    def run(self):
        print(self.Walker(1000, 100))
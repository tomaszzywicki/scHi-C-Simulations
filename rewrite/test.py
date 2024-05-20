import model_init

if __name__ == "__main__":
    walk = model_init.AdvancedSARW(10, 10)
    print(walk)
    walk.show_info()
    walk.plot()
from PIL import Image
import customtkinter as ctk
import matplotlib.pyplot as plt
import os
import io

import model_evolver
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class App(ctk.CTk):

    def __init__(self):
        super().__init__()
        self.title('Chromatine visualizer')
        self.geometry('1200x800')
        self.resizable(False, False)

        self.orange = "#FF8000"
        self.blue = "#47C7FC"

        self.file_path = "data/GSM1173493_cell-1.txt"
        self.data = model_evolver.scData(self.file_path)

        self.available_chromosomes = self.data.get_chromosome_list()
        self.available_resolutions = ["1e5", "5e5", "1e6", "5e6", "1e7"]
        
        self.generate_widgets()

    def generate_widgets(self):
        self.label = ctk.CTkLabel(self, text='Chromatine visualizer', font=("Open Sans", 50, "bold"), text_color=self.orange)
        self.label.grid(row=0, column=0, columnspan=2, padx=40, pady=40, sticky='w')

        self.data_frame = ctk.CTkFrame(self)
        self.data_frame.grid(row=1, column=0, padx=20, sticky='wn')

        self.new_model_label = ctk.CTkLabel(self.data_frame, text='New model', font=("Open Sans", 20, "bold"), text_color=self.blue)
        self.new_model_label.grid(row=0, column=0, columnspan=2, padx=20, pady=10, sticky='w')

        self.chromosome_label = ctk.CTkLabel(self.data_frame, text='Select chromosome')
        self.chromosome_label.grid(row=1, column=0, padx=20, pady=10, sticky='w')
        self.chromosome_box = ctk.CTkComboBox(self.data_frame, values=self.available_chromosomes, command=self.on_user_choice_changed)
        self.chromosome_box.grid(row=1, column=1, padx=20, pady=10)

        self.resolution_label = ctk.CTkLabel(self.data_frame, text='Select resolution')
        self.resolution_label.grid(row=2, column=0, padx=20, pady=10, sticky='w')
        self.resolution_box = ctk.CTkComboBox(self.data_frame, values=self.available_resolutions, command=self.on_user_choice_changed)
        self.resolution_box.grid(row=2, column=1, padx=20, pady=10)
        self.resolution_box.set("1e6")

        self.chromosome_box.bind("<Return>", self.on_user_choice_changed)
        self.resolution_box.bind("<Return>", self.on_user_choice_changed)

        self.generate_data_button = ctk.CTkButton(self.data_frame, command=self.generate_data, font=("Open Sans", 14, "bold"))
        self.generate_data_button.grid(row=3, column=0, columnspan=2, padx=20, pady=10, sticky="ew")

        self.data_plotted = ctk.CTkLabel(self.data_frame, text=None, width=300, height=300)
        self.data_plotted.grid(row=4, column=0, columnspan=2, sticky='we')

        self.create_model_button = ctk.CTkButton(self.data_frame, text="Create model", command=self.create_model, font=("Open Sans", 14, "bold"))
        self.create_model_button.grid(row=5, column=0, columnspan=2, padx=20, pady=10, sticky="ew")

        self.on_user_choice_changed(None)

        self.generate_model_frame()

    def generate_model_frame(self):
        self.model_frame = ctk.CTkFrame(self)
        self.model_frame.grid(row=1, column=1, padx=20, sticky='wn')

        self.button_frame = ctk.CTkFrame(self.model_frame)
        self.plot_frame = ctk.CTkFrame(self.model_frame)
        self.button_frame.grid(row=1, column=0, padx=20, pady=10, sticky='wn')
        self.plot_frame.grid(row=1, column=1, padx=20, pady=10, sticky='wn')

        self.train_model_label = ctk.CTkLabel(self.model_frame, text='Train model', font=("Open Sans", 20, "bold"), text_color=self.blue)
        self.train_model_label.grid(row=0, column=0, columnspan=2, padx=20, pady=10, sticky='w')
        self.plot_initial_model_button = ctk.CTkButton(self.button_frame, text="Plot initial model", font=("Open Sans", 14, "bold"), fg_color=self.orange, text_color="black", command=self.plot_initial_model)
        self.plot_initial_model_button.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        self.plot_current_model_button = ctk.CTkButton(self.button_frame, text="Plot current model", font=("Open Sans", 14, "bold"), fg_color=self.orange, text_color="black", command=self.plot_model)
        self.plot_current_model_button.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.plot_score_history_button = ctk.CTkButton(self.button_frame, text="Plot score history", font=("Open Sans", 14, "bold"), fg_color=self.orange, text_color="black", command=self.plot_score_history)
        self.plot_score_history_button.grid(row=2, column=0, padx=10, pady=10, sticky="ew")

        self.train_model_button = ctk.CTkButton(self.button_frame, text="Train model", font=("Open Sans", 14, "bold"), fg_color=self.orange, text_color="black", command=self.train_model)
        self.train_model_button.grid(row=3, column=0, padx=10, pady=10, sticky="ew")

        self.train_iter_entry = ctk.CTkEntry(self.button_frame, placeholder_text="Number of iterations", font=("Open Sans", 14))
        self.train_iter_entry.grid(row=4, column=0, padx=10, pady=10, sticky="ew")
        
        self.canvas = ctk.CTkCanvas(self.plot_frame, width=300, height=225)
        self.canvas.grid(row=1, column=0, columnspan=2, padx=20, pady=10, sticky="ew")
        self.zoom_button = ctk.CTkButton(self.plot_frame, text="Zoom", font=("Open Sans", 14, "bold"), fg_color=self.blue, text_color="black", command=self.zoom_model)
        self.zoom_button.grid(row=2, column=0, columnspan=2, padx=20, pady=10, sticky="ew")

    
    def zoom_model(self):
        self.model.plot()

    def on_user_choice_changed(self, event):
        chromosome = self.chromosome_box.get()
        resolution = self.resolution_box.get()
        file_path = f"data/chrom{chromosome}_res{resolution}.pkl"
        if not os.path.exists(file_path):
            self.data_not_ready()
        else:
            self.data = model_evolver.load(file_path=file_path)
            self.data_ready()

    def generate_data(self):
        chromosome = self.chromosome_box.get()
        resolution = self.resolution_box.get()
        if chromosome not in self.available_chromosomes:
            raise ValueError(f"Chromosome {chromosome} not available")
        if not resolution.replace("e", "", 1).isdigit():
            raise ValueError(f"Resolution {resolution} not available. Must be a number.")
        file_path = f"data/chrom{chromosome}_res{resolution}.pkl"
        resolution = int(float(resolution))
        self.data.prep(resolution=resolution, chromosome=chromosome)
        self.data.save(file_path=file_path)
        self.data_ready()

    def data_ready(self):
        self.data_plotted.configure(text="")
        self.generate_data_button.configure(state="disabled", fg_color="lightgray", text_color="white", text="Data ready")
        self.plot_data()
        self.create_model_button.configure(state="normal", fg_color=self.blue, hover_color="lightblue", text_color="black")

    def data_not_ready(self):
        self.generate_data_button.configure(state="normal", fg_color=self.orange, hover_color="darkorange", text_color="black", text="Generate data")
        empty_photo = ctk.CTkImage(Image.new("RGB", (300, 225), "#444444"), size=(300,225))
        self.data_plotted.configure(image=empty_photo, text="No data file found")
        self.data_plotted.image = empty_photo
        print("not ready")
        self.create_model_button.configure(state="disabled", fg_color="lightgray", text_color="white")

    def plot_data(self):

        buf = io.BytesIO()
        fig = self.data.scHiC_matrix(show=False)
        fig.savefig(buf, format='png')
        plt.close()
        buf.seek(0)

        image = Image.open(buf)
        photo = ctk.CTkImage(image, size=(300,225))
        self.data_plotted.configure(image=photo)
        self.data_plotted.image = photo

    def create_model(self):
        self.model = model_evolver.Model(self.data)

    def run(self):
        self.mainloop()
    
    def plot_model(self, fig=None):
        if fig is None:
            fig = self.model.plot(show=False)
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.canvas_fig = fig
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=1, column=0, columnspan=2, padx=20, pady=10, sticky="ew")
        plt.close()

    def plot_initial_model(self):
        self.plot_model(self.model.init_walk.plot(show=False))

    def plot_score_history(self):
        self.plot_model(self.model.plot_score_history(show=False))

    def train_model(self):
        try:
            iterations = int(self.train_iter_entry.get())
        except ValueError:
            print("Invalid number of iterations")
            iterations = 100
        self.model.evolve_simulated_annealing(iterations=iterations)
        self.plot_model()
    
# class AskForName(ctk.CTk):

#     def __init__(self, master=None):
#         super().__init__()
#         self.master = master
#         self.title("Create new model")
#         screen_width = self.winfo_screenwidth()
#         screen_height = self.winfo_screenheight()
#         width = 300
#         height = 200
#         x = (screen_width - width) // 2
#         y = (screen_height - height) // 2
#         self.geometry(f"{width}x{height}+{x}+{y}")
#         self.resizable(False, False)

#         self.canelled = True

#         self.name_entry = ctk.CTkEntry(self, font=("Open Sans", 14), placeholder_text="Model name")
#         self.name_entry.grid(row=0, column=0, columnspan=2, sticky="ew")

#         self.cancel_button = ctk.CTkButton(self, text="Cancel", command=self.after(100, self.destroy), font=("Open Sans", 14, "bold"))
#         self.cancel_button.grid(row=1, column=0, sticky="ew")
#         self.create_button = ctk.CTkButton(self, text="Create", command=self.create_model, font=("Open Sans", 14, "bold"))
#         self.create_button.grid(row=1, column=1, sticky="ew")
    
#     def create_model(self):

#         name = self.name_entry.get()
#         print(name)

#         self.after(100, self.destroy)

#     def run(self):
#         self.mainloop()

# class ConfirmPopup(ctk.CTk):
    
#     def __init__(self, message="Are you sure?"):
#         super().__init__()
#         self.title("Confirm")
#         self.geometry("300x200")
#         self.resizable(False, False)
#         self.confirmation = False

#         self.label = ctk.CTkLabel(self, text=message, font=("Open Sans", 14))
#         self.label.grid(row=0, column=0, columnspan=2, padx=20, pady=20)

#         self.no_button = ctk.CTkButton(self, text="No", command=self.no)
#         self.no_button.grid(row=1, column=0)
#         self.yes_button = ctk.CTkButton(self, text="Yes", command=self.yes)
#         self.yes_button.grid(row=1, column=1)

#     def no(self):
#         self.destroy()
    
#     def yes(self):
#         self.confirmation = True
#         self.destroy()

if __name__ == '__main__':
    app = App()
    app.run()

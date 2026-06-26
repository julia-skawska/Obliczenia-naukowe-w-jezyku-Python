#Imports---------------------------------------------------------------------

import tkinter as tk
from tkinter import messagebox as ms
from tkinter import ttk
from pathlib import Path
import scanpy as sc
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import pyplot as plt
from tkinter import PhotoImage

#Tests-----------------------------------------------------------------------
def set_item(app, path, text):
    path = Path(path)

    if path.exists():
        try:
            app.iconbitmap(str(path))
        except Exception as e:
            ms.showwarning("Icon Error", str(e))
    else:
        ms.showwarning(text, f"Not found:\n{path}")

def check_paths(path_input, text):
    path = Path(path_input.strip().strip('"').strip("'"))

    if path.exists():
        return path
    else:
        ms.showwarning(f"{text}", f"Not found:\n{path}")
        return None
    
#Plotting-----------------------------------------------------------------------
class PlotFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master, bg=controller.bg)

        self.controller = controller
        self.canvas = None

        self.current_index = 0
        self.group_keys = []

        nav = tk.Frame(self, bg=controller.bg)
        nav.pack(side="top", pady=10)

        tk.Button(nav, text="<- Prev",
                  command=self.prev_plot,
                  **controller.button_style).pack(side="left", padx=10)

        tk.Button(nav, text="Next ->",
                  command=self.next_plot,
                  **controller.button_style).pack(side="left", padx=10)

        tk.Button(self, text="End program",
                  command=self.confirm_exit,
                  **controller.button_style).pack(side="bottom", pady=10)
        
        tk.Button(self, text="Run again",
                  command=self.go_back,
                  **controller.button_style).pack(side="bottom", pady=10)

    def prepare_groups(self):
        adata = self.controller.adata
        sort_col = self.controller.sort_col

        self.group_keys = sorted(adata.obs[sort_col].astype(str).unique())

    def render_plot(self):
        adata = self.controller.adata
        genes = self.controller.valid_genes
        sort_col = self.controller.sort_col

        if not genes:
            return

        if sort_col not in adata.obs:
            ms.showerror("Error", f"'{sort_col}' not found in adata.obs")
            return

        if not self.group_keys:
            self.prepare_groups()

        group = self.group_keys[self.current_index]

        if "x_centroid" in adata.obs and "y_centroid" in adata.obs:
            x = adata.obs["x_centroid"].values
            y = adata.obs["y_centroid"].values
        else:
            ms.showerror("Error", "Missing x_centroid / y_centroid")
            return

        mask_group = (adata.obs[sort_col].astype(str) == group).values
        cmap = plt.get_cmap("tab10")
        gene_colors = {g: cmap(i % 10) for i, g in enumerate(genes)}

        gene_expr = {}
        for g in genes:
            if g in adata.var_names:
                v = adata[:, g].X
                v = v.toarray().ravel() if hasattr(v, "toarray") else np.asarray(v).ravel()
                gene_expr[g] = v
            else:
                gene_expr[g] = np.zeros(adata.n_obs)

        fig = Figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

        ax.scatter(
            x[mask_group],
            y[mask_group],
            s=0.8,
            c="lightgray",
            alpha=0.4)

        for g in genes:
            expr = gene_expr[g]
            mask_gene = mask_group & (expr > 0)

            ax.scatter(
                x[mask_gene],y[mask_gene],
                s=3,color=gene_colors[g],label=g)

        ax.set_title(f"{sort_col}: {group}")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.axis("off")
        fig.legend(loc="center right", frameon=False)

        if self.canvas:
            self.canvas.get_tk_widget().destroy()

        self.canvas = FigureCanvasTkAgg(fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack()

        adata_name = Path(self.controller.input_path).stem
        gene_string = "_".join(genes)
        safe_group = str(group).replace(" ", "_").replace("/", "_")

        filename = f"{adata_name}_{sort_col}_{safe_group}_{gene_string}.png"
        save_path = Path(self.controller.output_path) / filename

        fig.savefig(save_path,dpi=300,bbox_inches="tight")
        plt.close(fig)

    def next_plot(self):
        if not self.group_keys:
            return
        self.current_index = (self.current_index + 1) % len(self.group_keys)
        self.render_plot()

    def prev_plot(self):
        if not self.group_keys:
            return
        self.current_index = (self.current_index - 1) % len(self.group_keys)
        self.render_plot()

    def go_back(self):
        self.controller.valid_genes = []

        self.controller.xenium_frame.genes.clear()
        self.controller.xenium_frame.gene_listbox.delete(0, tk.END)

        self.controller.show_frame(self.controller.xenium_frame)

    def confirm_exit(self):
        if ms.askyesno("Confirm exit", "Are you sure you want to quit?"):
            self.controller.destroy()

#Acessing Xenium-----------------------------------------------------------------------
class XeniumFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master, bg=controller.bg)

        self.controller = controller

        self.is_human = tk.BooleanVar(value=True)
        
        tk.Checkbutton(
            self,
            text="is Homo Sapiens",
            variable=self.is_human,
            fg=controller.accent, bg=controller.bg,
            activeforeground=controller.accent, activebackground=controller.bg,
            selectcolor=controller.bg,).pack(pady=10)
        
        tk.Label(
            self,
            text="Genes of interest:",
            fg=controller.accent, bg=controller.bg).pack(pady=(10, 0))
        
        self.gene_entry = tk.Entry(self, width=25)
        self.gene_entry.pack(pady=5)

        tk.Button(
            self,
            text="Add gene",
            command=self.add_gene,
            **controller.button_style).pack(pady=5)

        self.gene_listbox = tk.Listbox(self, height=5, width=25)
        self.gene_listbox.pack(pady=10)
        self.genes = []

        tk.Button(
            self,
            text="End program",
            command=self.confirm_exit,
            **controller.button_style).pack(side="bottom", pady=10)

        tk.Button(
            self,
            text="Back to Setup",
            command=self.go_to_setup,
            **controller.button_style).pack(side="bottom", pady=10)

        tk.Button(
            self,
            text="Create plot",
            command=self.create_plot,
            **controller.button_style).pack(side="bottom", pady=10)

    def add_gene(self):
        gene = self.gene_entry.get().strip()

        if gene and gene not in self.genes:
            self.genes.append(gene)
            self.gene_listbox.insert(tk.END, gene)

        self.gene_entry.delete(0, tk.END)

    def check_genes(self, adata, genes):
        if self.is_human.get():
            genes = [g.upper() for g in genes]
        else:
            genes = [g.lower() for g in genes]

        valid = [g for g in genes if g in adata.var_names]
        missing = [g for g in genes if g not in adata.var_names]

        if missing:
            ms.showwarning(
                "Missing genes",
                "\n".join(missing))

        return valid
    
    def create_plot(self):
        adata = self.controller.adata
        valid_genes = self.check_genes(adata, self.genes)

        if not valid_genes:
            ms.showerror("Error", "No valid genes selected")
            return

        self.controller.valid_genes = valid_genes
        self.controller.adata = adata

        self.controller.plot_frame.current_index = 0
        self.controller.show_frame(self.controller.plot_frame)
        self.controller.plot_frame.render_plot()

    def go_to_setup(self):
        self.controller.show_frame(self.controller.setup_frame)

    def confirm_exit(self):
        if ms.askyesno(
            "Confirm exit",
            "Are you sure you want to quit?"
        ):
            self.controller.destroy()

#HomeScreen------------------------------------------------------------------
class StartFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master, bg=controller.bg)

        self.controller = controller

        tk.Label(
            self,
            text="Welcome to Xenium Visualisation",
            font=controller.title_font,
            fg=controller.accent,
            bg=controller.bg).pack(pady=5)

        tk.Label(
            self,
            text="Let's visualise your gene expression together",
            font=controller.normal_font,
            fg=controller.accent,
            bg=controller.bg).pack(pady=10)

        tk.Button(
            self,
            text="Begin Setup",
            command=lambda: self.controller.go_to(self.controller.input_frame),
            **controller.button_style).pack()

#Setup-----------------------------------------------------------------------
class SetupFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master, bg=controller.bg)

        self.controller = controller

        tk.Label(
            self,
            text="Sorting / grouping column:",
            font=controller.normal_font,
            fg=controller.accent,
            bg=controller.bg
        ).pack(pady=(20, 5))

        self.sort_var = tk.StringVar()

        self.sort_box = ttk.Combobox( self, textvariable=self.sort_var, state="readonly", width=40 ) 
        self.sort_box.pack()

        tk.Button(
            self,
            text="Begin plotting",
            command=self.go_next,
            **controller.button_style).pack(pady=20)

    def set_default_sort(self, cols):
        preferred = ["cell_type_cluster","celltype","annotation","condition","sample","patient"]

        for p in preferred:
            if p in cols:
                self.sort_var.set(p)
                return

        if cols:
            self.sort_var.set(cols[0])

    def go_next(self):
        self.controller.sort_col = self.sort_var.get()

        self.controller.plot_frame.group_keys = []
        self.controller.plot_frame.current_index = 0

        self.controller.show_frame(self.controller.xenium_frame)

#User input-----------------------------------------------------------------------
class InputFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master, bg=controller.bg)

        self.controller = controller

        tk.Label(
            self,
            text="Dataset Path:",
            font=controller.normal_font,
            fg=controller.accent,
            bg=controller.bg).pack(pady=(20, 5))

        self.dataset_entry = tk.Entry(self,width=40, fg=controller.accent, font=controller.normal_font)
        self.dataset_entry.pack()

        tk.Label(
            self,
            text="Output Directory:",
            font=controller.normal_font,
            fg=controller.accent,
            bg=controller.bg).pack(pady=(20, 5))

        self.output_entry = tk.Entry(self,width=40, fg=controller.accent, font=controller.normal_font)
        self.output_entry.pack()

        tk.Button(
            self,
            text="Run",
            command=self.user_input,
            **controller.button_style).pack(pady=20)

    def user_input(self):
        dataset = self.dataset_entry.get()
        output = self.output_entry.get()

        dataset_path = check_paths(dataset, "Path to dataset: ")
        output_path = check_paths(output, "Path to output: ")

        if dataset_path is None or output_path is None:
            return

        self.controller.input_path = str(dataset_path)
        self.controller.output_path = str(output_path)

        try:
            self.controller.adata = sc.read_h5ad(self.controller.input_path)
        except Exception as e:
            ms.showerror("Load Error", str(e))
            return
        
        obs_cols = list(self.controller.adata.obs.columns)

        self.controller.setup_frame.sort_box["values"] = obs_cols
        self.controller.setup_frame.set_default_sort(obs_cols)

        self.controller.show_frame(self.controller.setup_frame)

#Tkinter-----------------------------------------------------------------------
class Xenium_Visualiser(tk.Tk):
    def __init__(self, assets_dir):
        super().__init__()

        self.title("Xenium Visualisation")
        self.geometry("600x400")
        icon_path = assets_dir / "Xvis.ico"
        png_path = assets_dir / "Xvis.png"
        self.iconbitmap(icon_path)
        icon = PhotoImage(file=png_path)

        self.iconphoto(True, icon)
        self._icon = icon
        
        self.input_path = None
        self.output_path = None
        self.adata = None

        self.setup_style()
        self.container = tk.Frame(self, bg=self.bg)
        self.container.pack(fill="both", expand=True)
        self.start_frame = StartFrame(self.container, self)
        self.input_frame = InputFrame(self.container, self)
        self.setup_frame = SetupFrame(self.container, self)
        self.xenium_frame = XeniumFrame(self.container, self)
        self.plot_frame = PlotFrame(self.container, self)
        self.start_frame.pack(fill="both", expand=True)

    def setup_style(self):
        self.bg = "#ADC7ED"
        self.fg = "#1F3A5F"
        self.accent = "#7C5AED"

        self.title_font = ("Calibri", 20, "bold")
        self.normal_font = ("Calibri", 12)

        self.button_style = {
            "fg": "#7C5AED",
            "bg": "#ADDBED",
            "relief": "ridge",
            "activeforeground": "#7C5AED",
            "activebackground": "#ADDBED",
            "bd": 5,
            "font": self.normal_font}

        self.configure(bg=self.bg)
        style = ttk.Style(self)
        style.theme_use("default")
        style.configure(".", font=self.normal_font)
        style.configure("TLabel", foreground=self.fg)
        style.configure("TButton", font=self.normal_font)

    def show_frame(self, frame_to_show):
        for frame in (self.start_frame, self.input_frame, self.setup_frame, self.xenium_frame, self.plot_frame):
            frame.pack_forget()

        frame_to_show.pack(fill="both", expand=True)
    
    def go_to(self, frame):
        self.show_frame(frame) 

#main-----------------------------------------------------------------------
if __name__ == "__main__":
    base_dir = Path(__file__).resolve().parent
    assets_dir = base_dir / "Images_For_UI"
    app = Xenium_Visualiser(assets_dir)
    app.mainloop()
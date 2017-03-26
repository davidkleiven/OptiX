import numpy as np
import matplotlib as mpl
mpl.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib.backends import backend_tkagg as btk
import tkinter as tk

class PlotFrame( tk.Frame ):
    def __init__( self, parent, controller ):
        tk.Frame.__init__( self, parent )
        self.data = []

        self.fig1 = plt.figure(figsize=(5,5))
        self.ax1 = self.fig1.add_subplot(1,1,1)
        self.canvas = btk.FigureCanvasTkAgg(self.fig1, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=0,column=0)


def main():
    root = tk.Tk()
    root.wm_title("Extract Unit Cell")
    mainframe = tk.Frame()
    pFrame = PlotFrame(mainframe, root)
    pFrame.tkraise()
    root.mainloop()

if __name__ == "__main__":
    main()

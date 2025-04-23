import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import PySimpleGUI as psg
import matplotlib
matplotlib.use('TkAgg')


def draw_figure(canvas, figure):
   tkcanvas = FigureCanvasTkAgg(figure, canvas)
   tkcanvas.draw()
   tkcanvas.get_tk_widget().pack(side='top', fill='both', expand=1)
   return tkcanvas



def get_figure(a,b):
    groups = [np.arange(0,40), np.arange(0,10), np.arange(10,20), np.arange(20,30), np.arange(30,40)]
    
    xmin, xmax = np.min([np.min(a),np.min(b)]), np.max([np.max(a),np.max(b)])
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].hist(a[groups[1]])
    axs[0, 0].hist(b[groups[1]])
    axs[0, 0].set_title('aSYN')
    axs[0, 1].hist(a[groups[2]])
    axs[0, 1].hist(b[groups[2]])
    axs[0, 1].set_title('comb.')
    axs[1, 0].hist(a[groups[3]])
    axs[1, 0].hist(b[groups[3]])
    axs[1, 0].set_title('INFg')
    axs[1, 1].hist(a[groups[4]])
    axs[1, 1].hist(b[groups[4]])
    axs[1, 1].set_title('UT')
    for i in range(2):
        for j in range(2):
            axs[i,j].set_xlim([xmin,xmax])
    fig.tight_layout()
    return fig

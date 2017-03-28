import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg 
mpl.rcParams['pdf.fonttype'] = 42
import os
import numpy as np

def density_plot(x, y, nbins=42, log=False):
    mask = (~np.isnan(x)) & (~np.isnan(y))
    x = x[mask]
    y = y[mask]
    H, xedges, yedges = np.histogram2d(x,y,bins=nbins)
    ix = np.searchsorted(xedges, x)
    ix[ix == nbins] = nbins - 1
    iy = np.searchsorted(yedges, y)
    iy[iy == nbins] = nbins - 1
    v = H[ix, iy]
    i = v.argsort()
    cc = v[i]
    if log:
        cc = np.log(cc + 1)
    plt.scatter(x[i], y[i], c=cc, s=3, edgecolor='')


def plot_better(width=10, height=5, grid='xy', legend=False, visible_axes=True):
    plt.figure(figsize=(width, height))
    ax = plt.subplot(111)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",  
            labelbottom="on", left="off", right="off", labelleft="on")
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False) 
    if visible_axes:
        ax.spines["bottom"].set_visible(True) 
        ax.spines["bottom"].set_color('gray') 
        ax.spines["left"].set_visible(True)   
        ax.spines["left"].set_color('gray')
    else:
        ax.spines["bottom"].set_visible(False)  
        ax.spines["left"].set_visible(False) 
    
    if grid == 'xy':
        ax.xaxis.grid(True) 
        ax.yaxis.grid(True) 
    if grid == 'x':
        ax.xaxis.grid(True) 
    if grid == 'y':
        ax.yaxis.grid(True) 
    if legend:
        plt.legend(loc=2, bbox_to_anchor=(1.05, 1), frameon=False)

    return ax
    


def improve_plot(ax, grid='xy', legend=False, visible_axes=True):
    plt.tick_params(axis="both", which="both", bottom="off", top="off",  
            labelbottom="on", left="off", right="off", labelleft="on")
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False) 
    if visible_axes:
        ax.spines["bottom"].set_visible(True) 
        ax.spines["bottom"].set_color('gray') 
        ax.spines["left"].set_visible(True)   
        ax.spines["left"].set_color('gray')
    else:
        ax.spines["bottom"].set_visible(False)  
        ax.spines["left"].set_visible(False) 
    
    if grid == 'xy':
        ax.xaxis.grid(True) 
        ax.yaxis.grid(True) 
    if grid == 'x':
        ax.xaxis.grid(True) 
    if grid == 'y':
        ax.yaxis.grid(True) 
    if legend:
        plt.legend(loc=2, bbox_to_anchor=(1.05, 1), frameon=False)

    return ax


class Counter:

    def __init__(self, initial_count=0):
        self.count = initial_count


    def get_number(self, as_string=True):
        self.count += 1
        if as_string:
            return '%02d' % self.count
        return self.count

    def reset(self, initial_count=0):
        self.count = initial_count


def save_image(image_counter, title, folder, prefix):
    figure_name = '%s_img%s_%s.png' %(prefix, image_counter.get_number(), '_'.join(title.split()))
    plt.savefig(os.path.join(folder, figure_name), dpi=300)



def get_colors(number_of_colors, colormap):
    step = 255./(number_of_colors-1)
    return [colormap(int(i*step)) for i in range(number_of_colors)]


def hex_to_rgb(value):
    """Return (red, green, blue) for the color given as #rrggbb."""
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)


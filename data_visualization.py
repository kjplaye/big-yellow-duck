import os
import tempfile
import ctypes
import numpy as np

import distinctipy
import matplotlib.patches as Patch
from ggobi import *
from imagesc import imagesc
from matplotlib.colors import ListedColormap
from wview import wview


def plotn(x):
    """
    Hopefully we get a better version soon
    """
    if type(x) == str:
        os.system('plotn '+x)
    else:
        f = tempfile.NamedTemporaryFile()
        for L in x:
            f.write(' '.join(map(str,L)) + '\n')
        f.flush()
        os.system('plotn '+f.name)
        f.close()
    
        
def imshow_legend(X, colors, names, aspect_ratio):
    """Displays an image X with values from 0 to n-1, using a list of n colors and a list of n names as the legend.
    
    Args:
    - X: 2D numpy array with values in the range 0 to n-1.
    - colors: n by 3 array of RGB colors (each value should be in [0, 1]).
    - names: List of n strings representing the names for each color in the legend.
    - aspect_ratio: Aspect ratio for the plot.
    """
    if isinstance(colors, int):
        colors = distinctipy.get_colors(colors)

    # Create a colormap from the colors
    cmap = ListedColormap(colors)

    # Plot the array using imshow
    plt.imshow(X, cmap=cmap, aspect=aspect_ratio * X.shape[1] / X.shape[0], interpolation = 'none')

    # Create a legend with patches for each color
    legend_patches = [Patch(color=colors[i], label=names[i]) for i in range(len(names))]
    plt.legend(handles=legend_patches, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    # Display the plot
    plt.show()
    


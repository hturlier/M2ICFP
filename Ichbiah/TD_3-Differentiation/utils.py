#Quelques librairies utiles pour le TD
from functools import partial 
from collections import defaultdict
import math
import numpy as np # Numerical computing library
import matplotlib
import matplotlib.pyplot as plt # Plotting library
import scipy.integrate #Integration library
from mpl_toolkits.mplot3d import axes3d #Used for the 3d bifurcation plot
import matplotlib.patches as mpatches #used to write custom legends
import sympy as sm
from tqdm import tqdm
import pylab as p
import sys
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection



from matplotlib.colors import LinearSegmentedColormap

colors = [(0, 0, 1), (1, 0, 0), (0, 1, 0)]  # B -> R -> G

n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins

# https://matplotlib.org/3.1.0/gallery/color/custom_cmap.html

cm = LinearSegmentedColormap.from_list(
    'rgb', colors, N=15)

def add_arrow(line,  direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    start_ind = 150
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )
    
    start_ind = 200
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )
    
    start_ind = 250
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )


def plot_voronoi_setup(figsize=(10,10),alpha = 0.3,plot_bassins=False,markersize = 200):
    ds = np.array([[-0.86156, -0.49742],
       [-0.     ,  0.99485],
       [ 0.86156, -0.49742]])
    d1,d2,d3 = ds

    points = np.array([d1,d2,d3])

    u1 = d2 - d1
    u2 = d3 - d2 
    u3 = d3 - d1
    def orth(u):
        return(np.array([u[0],-u[1]]))
    n1 = orth(u1)
    n2 = orth(u2)
    n3 = orth(u3)

    meanp = np.mean(points,axis=0)
    l = (np.sqrt(2)-meanp[1])/n2[1]
    vorr = (n2[0]*l + meanp[0])

    fig, ax = plt.subplots(figsize=figsize)

    patches = []
    polygon_1 = Polygon([meanp,[np.sqrt(2),vorr],[np.sqrt(2),np.sqrt(2)],[-np.sqrt(2),np.sqrt(2)],[-np.sqrt(2),vorr]], True)
    polygon_2 = Polygon([[0,0],meanp,[np.sqrt(2),vorr],[np.sqrt(2),-np.sqrt(2)],[0,-np.sqrt(2)]], True)
    polygon_3 = Polygon([[0,0],meanp,[-np.sqrt(2),vorr],[-np.sqrt(2),-np.sqrt(2)],[0,-np.sqrt(2)]], True)
    #polygon_3 = Polygon([[0.5,0],meanp,[0,vorr],[0,0]], True)

    patches.append(polygon_1)
    patches.append(polygon_2)
    patches.append(polygon_3)

    colors = np.array([0.5,1,0])#np.array([89.20561225, 43.00103067, 23.64260471])
    p = PatchCollection(patches, alpha=alpha,cmap=matplotlib.cm.jet)
    p.set_array(colors)
    ax.add_collection(p)

    if plot_bassins : 
        ax.scatter(meanp[0],meanp[1],color = 'black',s=markersize)
        ax.scatter(d1[0],d1[1],color = 'black',s=markersize)
        ax.scatter(d2[0],d2[1],color = 'black',s=markersize)
        ax.scatter(d3[0],d3[1],color = 'black',s=markersize)

    ax.set_xlim(-np.sqrt(2),np.sqrt(2))
    ax.set_ylim(-np.sqrt(2),np.sqrt(2))
    plt.axis("off")
    return fig,ax
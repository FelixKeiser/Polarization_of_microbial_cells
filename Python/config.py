# This file contains the configuration for the project
import numpy as np
import pandas as pd
import copy
import os
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.markers as mmarkers
import matplotlib.cm as cm


# Using seaborn's style
plt.style.use('seaborn')
# With LaTex fonts
plt.style.use('tex')
tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 12,
    "font.size": 12,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10
}
width = "thesis"
plt.rcParams['axes.grid'] = True
marker_style = mmarkers.MarkerStyle('o', fillstyle='none')
msize = 15
figdir = r"..\Latex\Polarization_of_Bacterial_Cells\Images"
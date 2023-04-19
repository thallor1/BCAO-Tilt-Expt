import numpy as np
import os
import sys
import h5py
import numpy as np
from numpy import array as npa
import pandas as pd
from itertools import product
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
from matplotlib.patches import RegularPolygon
from scipy.interpolate import interp2d, griddata
from matplotlib.tri import Triangulation, CubicTriInterpolator, LinearTriInterpolator, TriInterpolator
import string
sys.path.append('python scripts/')
from ezhang_funcs import *

# First read in the MD calculation, test
# read data from hdf5 file
f = h5py.File("MD Calcultions/1K_data/configuration_0.h5")
#f = h5py.File("MD Calcultions/0.1K_data/configuration_0.h5")
Suv = npa(f["spin_correlations/S_qw"])
freq = npa(f["spin_correlations/freq"])
ks = npa(f["spin_correlations/momentum"])
pc = npa(f["spin_correlations/pc"])
T = f.attrs["T"]
J1x = f.attrs["J1x"]
J1z = f.attrs["J1z"]
J3x = f.attrs["J3x"]
J3z = f.attrs["J3z"]
D = f.attrs["D"]
E = f.attrs["E"]
print(np.shape(freq))
betaomega = 1/T * freq
f.close()
# mask singularities
omit=[12, 19, 20, 22, 44, 51]
shift = 1
for i, p in enumerate(omit):
    cond = p < pc
    pc[cond] -= shift

DSF = structure_factor(Suv, np.c_[ks, np.zeros(ks.shape[0])])
x = ks[:,0]/(2*np.pi)
xs = ks[~np.isin(np.arange(ks.shape[0]),omit),0] / (2*np.pi)
DSF = DSF[:,~np.isin(np.arange(ks.shape[0]),omit)]

# splice data
d = (DSF.T*betaomega).T
d = d[freq>0.5, :]
xs = np.arange(xs.shape[0])

# Make a plot
fig, ax = plt.subplots()
c = "Spectral_r"
cmap = copy.copy(mpl.cm.get_cmap(c))
c1 = ax.pcolormesh(xs, freq[freq>0.5], (d/d.max()), cmap=cmap, shading="auto")
fig.show()
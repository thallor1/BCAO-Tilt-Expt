import glob

import numpy as np
sys.path.append('python scripts/')

import ezhang_funcs as ez
from make_spgfigax import make_bcao_spgfigax

fig, axs = make_bcao_spgfigax()

files = glob.glob('Spag_ascii/*7T.txt')
# order is m1, k1, g2, m2, k1, m3, g3, m1
files = [files[4], files[2], files[0], files[5], files[3], files[6], files[1]]

# iterate through the files and plot
for i, ax in enumerate(axs):
    file = files[i]
    dat = np.genfromtxt(file, skip_header=2)
    f = open(file)
    lines = f.readlines()
    f.close()
    shape = lines[1].split(' ')[-1]
    shape = shape.replace('\n', '')
    shape = shape.split('x')

    # File format is I, Err, HKL, Omega
    hkl = np.unique(dat[:, 2])
    energy = np.unique(dat[:, 3])
    # Change to bin edges
    hkl = np.append(hkl - np.abs(hkl[1] - hkl[0]) / 2, np.nanmax(hkl) + np.abs(hkl[1] - hkl[0]) / 2)
    energy = np.append(energy - np.abs(energy[1] - energy[0]) / 2,
                       np.nanmax(energy) + np.abs(energy[1] - energy[0]) / 2)

    Q, E = np.meshgrid(hkl, energy)
    Intensity = np.reshape(dat[:, 0], (int(shape[0]), int(shape[1])))
    labels = ez.get_formatted(["M1", "K1", "G2", "M2", "K1", "M3", "G3", "M1"])
    ax.set_xticks([np.nanmin(dat[:, 2])])
    ax.set_xticklabels([labels[i]])
    if i == len(axs) - 1:
        ax.set_xticks([np.nanmin(dat[:, 2]), np.nanmax(dat[:, 2])])
        ax.set_xticklabels([labels[i], labels[i + 1]])
    ax.set_xlim(np.nanmin(dat[:, 2]), np.nanmax(dat[:, 2]))
    ax.set_ylim(0, 7)
    ax.pcolormesh(Q, E, Intensity.T, vmin=0, vmax=1, cmap='Spectral_r')
fig.show()

import numpy as np
import matplotlib.pyplot as plt
import ezhang_funcs as ez
import glob
from make_spgfigax import make_bcao_spgfigax, make_compare_bcao_spgfigax
from importMD import importMD
import scipy
import scipy.interpolate as interpolate

# First we will plot the measurement at 7 T on top.

fig, axsA, axsB = make_compare_bcao_spgfigax()

files = glob.glob("Spag_ascii/*7T.txt")
# order is m1, k1, g2, m2, k1, m3, g3, m1
files = [files[4], files[2], files[0], files[5], files[3], files[6], files[1]]

# Iterate through the files and plot
for i, ax in enumerate(axsA):
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
    Intensity[:, np.where(energy < 0.5)[0]] = np.nan
    labels = ez.get_formatted(["M1", "K1", "G2", "M2", "K1", "M3", "G3", "M1"])
    ax.set_xticks([np.nanmin(dat[:, 2])])
    ax.set_xticklabels([labels[i]])
    if i == len(axsA) - 1:
        ax.set_xticks([np.nanmin(dat[:, 2]), np.nanmax(dat[:, 2])])
        ax.set_xticklabels([labels[i], labels[i + 1]])
    ax.set_xlim(np.nanmin(dat[:, 2]), np.nanmax(dat[:, 2]))
    ax.set_ylim(0, 7)
    c1 = ax.pcolormesh(Q, E, Intensity.T, vmin=0, vmax=1, cmap='Spectral_r')
cbar = fig.colorbar(c1, ax=axsA[6], use_gridspec=True, label='I (a.u.)')
axsA[0].text(0.2, 0.95, '(a)', fontsize=10, transform=axsA[0].transAxes, horizontalalignment='left',
             verticalalignment='top')
# Import the MD calculation.

filename = "MD Calculations/1K_data/configuration_0.h5"
xs, freq, sqw, pc = importMD(filename)

# Import the LSWT calculation
spinw_dir = 'SpinW/Tilted_spectra_files/'
hkl_sw = np.genfromtxt(spinw_dir + 'BCAO_j1j3_hkl_tilt.csv', delimiter=',')
omega_sw = np.genfromtxt(spinw_dir + 'BCAO_j1j3_omega_tilt.csv', delimiter=',')
sqw_sw = np.genfromtxt(spinw_dir + 'BCAO_j1j3_swConv_tilt.csv', delimiter=',').T  # For some reason in transposed shape

# Indices of bin edges for spinW, done manually:
spinw_bin_i = [0, 299, 598, 897, 1196, 1495, 1794, 2093]
# Convert the hkl vectors to indices to match MD
hkl_sw_i = np.arange(0, np.shape(hkl_sw)[1], 1)

# Iterate through the slices and use an interpolation scheme to add the calculations.
for i, ax in enumerate(axsB):
    init_pci = int(pc[i])
    final_pci = int(pc[i + 1])
    init_spinwi = spinw_bin_i[i]
    final_spinwi = spinw_bin_i[i + 1]
    spinw_q_i = hkl_sw_i[init_spinwi:final_spinwi]
    spinw_sqw = sqw_sw[init_spinwi:final_spinwi, :]
    MD_Qpts = xs[init_pci:final_pci]
    MD_omegas = freq[freq > 0.5]
    MD_sqw = np.array(sqw)[:, init_pci:final_pci]
    # Make the spinw Qpts match the MD, both of length 1 and start from 0
    MD_Qpts = MD_Qpts - np.min(MD_Qpts)
    MD_Qpts = MD_Qpts / np.max(MD_Qpts)
    spinw_q_i = spinw_q_i - np.min(spinw_q_i)
    spinw_q_i = spinw_q_i / np.max(spinw_q_i)
    # Generate an interpolation scheme to add the lswt to the MD
    interp = interpolate.RegularGridInterpolator((spinw_q_i,
                                                  omega_sw[1:] - np.abs(omega_sw[1] - omega_sw[0]) / 2.0),
                                                 spinw_sqw,
                                                 bounds_error=False)
    hkl_interp, omega_interp = np.meshgrid(MD_Qpts, MD_omegas)
    Z = interp((hkl_interp, omega_interp))
    Z += MD_sqw * 1.6
    c1 = ax.pcolormesh(MD_Qpts, MD_omegas, Z * 1.0, cmap='Spectral_r', shading="auto", vmin=0, vmax=1.0,
                       rasterized=True)
    # c1 =ax.pcolormesh(hkl_interp,omega_interp,Z,vmin=0,vmax=0.5,cmap='Spectral_r',rasterized=True)
    ax.set_ylim(0, 7)
    labels = ez.get_formatted(["M1", "K1", "G2", "M2", "K1", "M3", "G3", "M1"])
    ax.set_xticks([0])
    ax.set_xticklabels([labels[i]])
    if i == len(axsB) - 1:
        ax.set_xticks([0, 1])
        ax.set_xticklabels([labels[i], labels[i + 1]])
    ax.set_xlim(0, 1)

cbar = fig.colorbar(c1, ax=ax, use_gridspec=True, label='I (a.u.)')
axsB[0].set_ylabel(r'$\hbar\omega$ (meV)')
axsB[0].text(0.2, 0.95, '(b)', fontsize=10, transform=axsB[0].transAxes, horizontalalignment='left',
             verticalalignment='top')

fig.savefig('bcao_vs_md_lswt.pdf',bbox_inches='tight',dpi=300)
fig.show()

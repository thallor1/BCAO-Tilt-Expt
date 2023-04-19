import sys

sys.path.append('python scripts/')
import ezhang_funcs as ez
import numpy as np
import matplotlib.pyplot as plt
from make_spgfigax import make_bcao_spgfigax

plt.close('all')
import matplotlib

# Read in the spinW calculation
spinw_dir = 'SpinW/Tilted_spectra_files/'
hkl_sw = np.genfromtxt(spinw_dir + 'BCAO_j1j3_hkl_tilt.csv', delimiter=',')
omega_sw = np.genfromtxt(spinw_dir + 'BCAO_j1j3_omega_tilt.csv', delimiter=',')
sqw_sw = np.genfromtxt(spinw_dir + 'BCAO_j1j3_swConv_tilt.csv', delimiter=',').T  # For some reason in transposed shape

# Indices of bin edges for spinW, done manually:
spinw_bin_i = [0, 299, 598, 897, 1196, 1495, 1794, 2093]
# Convert the hkl vectors to indices to match MD notation
hkl_sw_i = np.arange(0, np.shape(hkl_sw)[1], 1)

fig,axs = make_bcao_spgfigax()

for i, ax in enumerate(axs):
    init_spinwi = spinw_bin_i[i]
    final_spinwi = spinw_bin_i[i + 1]
    spinw_q_i = hkl_sw_i[init_spinwi:final_spinwi + 1]
    spinw_sqw = sqw_sw[init_spinwi:final_spinwi, :]
    HKL, Omega = np.meshgrid(spinw_q_i, omega_sw)
    c1 = ax.pcolormesh(HKL, Omega, spinw_sqw.T, vmin=0, vmax=0.5, cmap='Spectral_r', rasterized=True)
    ax.set_ylim(0, 7)
    labels = ez.get_formatted(["M1", "K1", "G2", "M2", "K1", "M3", "G3", "M1"])
    ax.set_xticks([init_spinwi])
    ax.set_xticklabels([labels[i]])
    if i == len(axs) - 1:
        ax.set_xticks([init_spinwi, final_spinwi])
        ax.set_xticklabels([labels[i], labels[i + 1]])
    ax.set_xlim(init_spinwi, final_spinwi)
cbar = fig.colorbar(c1, ax=axs[-1], use_gridspec=True, label='I (a.u.)')
axs[0].set_ylabel('$\hbar\omega$ (meV)')
fig.show()

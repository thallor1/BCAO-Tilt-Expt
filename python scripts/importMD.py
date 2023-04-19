import h5py
from numpy import array as npa
from ezhang_funcs import *


def importMD(filename):
    """

    :rtype: object
    """
    f = h5py.File(filename, mode='r')
    # f = h5py.File("MD Calcultions/0.1K_data/configuration_0.h5")
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
    betaomega = 1 / T * freq
    f.close()
    # mask singularities
    omit = [12, 19, 20, 22, 44, 51]
    shift = 1
    for i, p in enumerate(omit):
        cond = p < pc
        pc[cond] -= shift

    DSF = structure_factor(Suv, np.c_[ks, np.zeros(ks.shape[0])])
    x = ks[:, 0] / (2 * np.pi)
    xs = ks[~np.isin(np.arange(ks.shape[0]), omit), 0] / (2 * np.pi)
    DSF = DSF[:, ~np.isin(np.arange(ks.shape[0]), omit)]

    # splice data
    d = (DSF.T * betaomega).T
    d = d[freq > 0.5, :]
    xs = np.arange(xs.shape[0])

    # Return the relevant quantities.
    return xs, freq, d/d.max(), pc
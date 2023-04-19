import numpy as np


def structure_factor(spins, qs, basis=np.eye(3)):
    #     diag = spins[...,(0, 4, 8)].sum(axis=-1)
    cond = (qs[:, 0] == 0) & (qs[:, 1] == 0) & (qs[:, 2] == 0)
    qs[cond] = np.ones(3) * 1e-8
    qsquared = qs[:, 0] ** 2 + qs[:, 1] ** 2 + qs[:, 2] ** 2
    SF = np.zeros(spins.shape[:-1])

    for a in range(3):
        for b in range(3):
            projector = np.dot(basis[:, a], basis[:, b]) - (
                    np.einsum("ij,j", qs, basis[:, a]) * np.einsum("ij,j", qs, basis[:, b])) / qsquared
            SF += projector * spins[..., a + b + (2 * a)]

    return SF


def get_formatted(keys):
    for i, k in enumerate(keys):
        key = list(k)
        if key[0] == "G":
            keys[i] = r"$\Gamma_{}$".format(key[1])
        else:
            keys[i] = r"${}_{}$".format(key[0], key[1])
    return keys

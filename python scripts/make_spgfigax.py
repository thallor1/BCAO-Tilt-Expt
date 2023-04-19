import matplotlib
import matplotlib.pyplot as plt


def make_bcao_spgfigax():
    fig = plt.figure(figsize=(3.54 * 2, 9.489 / 4))
    spec = matplotlib.gridspec.GridSpec(1, 7, width_ratios=[0.28867, 0.5777, 0.5, 0.288675, 0.288675, 0.5, 0.5],
                                        hspace=0,
                                        wspace=0)
    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])
    ax2 = fig.add_subplot(spec[2])
    ax3 = fig.add_subplot(spec[3])
    ax4 = fig.add_subplot(spec[4])
    ax5 = fig.add_subplot(spec[5])
    ax6 = fig.add_subplot(spec[6])
    axs = [ax0, ax1, ax2, ax3, ax4, ax5, ax6]
    for i, ax in enumerate(axs):
        if i > 0:
            ax.tick_params(labelleft=False, labelright=False, labeltop=False, labelbottom=True, direction='in', left=True, right=True)
        else:
            ax.tick_params(labelleft=True, labelright=False, labeltop=False, labelbottom=True, direction='in', left=True, right=True)
    return fig, axs


def make_compare_bcao_spgfigax():
    fig = plt.figure(figsize=(3.54 * 2, 9.489 / 2))
    spec = matplotlib.gridspec.GridSpec(2, 7, width_ratios=[0.28867, 0.5777, 0.5, 0.288675, 0.288675, 0.5, 0.5],
                                        hspace=0.4, wspace=0)
    ax0a = fig.add_subplot(spec[0, 0])
    ax1a = fig.add_subplot(spec[0, 1])
    ax2a = fig.add_subplot(spec[0, 2])
    ax3a = fig.add_subplot(spec[0, 3])
    ax4a = fig.add_subplot(spec[0, 4])
    ax5a = fig.add_subplot(spec[0, 5])
    ax6a = fig.add_subplot(spec[0, 6])
    axs1 = [ax0a, ax1a, ax2a, ax3a, ax4a, ax5a, ax6a]
    ax0 = fig.add_subplot(spec[1, 0])
    ax1 = fig.add_subplot(spec[1, 1])
    ax2 = fig.add_subplot(spec[1, 2])
    ax3 = fig.add_subplot(spec[1, 3])
    ax4 = fig.add_subplot(spec[1, 4])
    ax5 = fig.add_subplot(spec[1, 5])
    ax6 = fig.add_subplot(spec[1, 6])
    axs2 = [ax0, ax1, ax2, ax3, ax4, ax5, ax6]
    for i, ax in enumerate(axs1):
        if i > 0:
            ax.tick_params(labelleft=False, labelright=False, labeltop=False, labelbottom=True, direction='in', left=True, right=True)
        else:
            ax.tick_params(labelleft=True, labelright=False, labeltop=False, labelbottom=True, direction='in', left=True, right=True)
    for ax in axs2:
        ax.tick_params(labelleft=False, labelright=False, labeltop=False, labelbottom=True, direction='in', left=True, right=True)
    ax0a.set_ylabel(r'$\hbar\omega$ (meV)')
    return fig, axs1, axs2

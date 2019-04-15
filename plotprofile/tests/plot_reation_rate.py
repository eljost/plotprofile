#!/usr/bin/env python3

# [1] 10.1021/acs.organomet.8b00456
#     https://pubs.acs.org/doi/10.1021/acs.organomet.8b00456

import numpy as np
import matplotlib.pyplot as plt

from plotprofile.main import reaction_rate


def plot_reaction_rate():
    """See Fig. 9 in [1]."""
    temps = np.array((273.15, 298.15, 373.15, 473.15, 573.15))
    acts = np.linspace(10, 45, 100)
    fig, ax = plt.subplots()
    for T in temps:
        acts_J_mol = acts * 1000 * 4.1868
        rrs = reaction_rate(acts_J_mol, T) * 3600# / 6.022e23
        ax.semilogy(acts, rrs, label=f"T={T:.2f} K")
    ax.set_xlim(10, 45)
    ax.set_ylim(1e-20, 1e+15)
    ax.set_xlabel("$E_A \quad / \quad kcal \cdot mol^{-1}$")
    ax.set_ylabel("$k \quad / \quad 1 \cdot s^{-1}$")
    ax.legend()
    plt.show()
    return


plot_reaction_rate()

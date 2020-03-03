import numpy as np

from plotprofile.constants import KB, H_PLANCK, R


class Reaction:

    def __init__(self, educts, ts, products, energies):
        self.educts = educts
        self.ts = ts
        self.products = products
        self.energies = energies

    def reaction_rate(self, temp=298.15, energy_key="G_solv_alt"):
        energies = self.energies[energy_key]
        educts, ts, products = energies
        barrier = ts - educts
        return KB*temp/H_PLANCK * np.exp(-barrier/(R*temp))

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

    def factor(self, mol):
        factor = 0

        # mol is used in reaction
        if mol in self.educts:
            factor = -1
        # mol is formed in reaction
        elif mol in self.products:
            factor = 1
        return factor

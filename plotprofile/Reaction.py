import numpy as np
import sympy as sym

from plotprofile.constants import AU2KJMOL, KB, H_PLANCK, R
from plotprofile.helpers import to_list


class Reaction:

    def __init__(self, name, educts, ts, products, energies=None, k=None,
                 add=None, label=None):
        self.name = name
        self.educts = to_list(educts)
        self.ts = to_list(ts)
        self.products = to_list(products)
        self.energies = energies
        self.k = k

    def get_back_reaction(self):
        educts = self.products
        products = self.educts

        energies = self.energies.copy()
        for k, v in energies.items():
            eds, ts, prod = v
            energies[k] = np.array((prod, ts, eds))
        name = f"{self.name} back"
        back_rx = Reaction(name, educts, self.ts, products, energies)
        return back_rx

    def reaction_rate(self, temp=298.15, barrier=None, energy_key="G_solv_alt"):
        """barrier must be in [J/mol] if given."""
        if barrier is None:
            energies = self.energies[energy_key].copy()
            energies -= energies.min()
            energies *= AU2KJMOL * 1000
            educts, ts, products = energies
            barrier = ts - educts
        return KB*temp/H_PLANCK * np.exp(-barrier/(R*temp))

    def get_expr(self, mol, mol_map, temp=298.15, k=None, energy_key="G_solv_alt",
                 k_thresh=1e-8):
        # Take stochiometry into account
        formed = sum([mol == reactant for reactant in self.products])
        used = sum([mol == reactant for reactant in self.educts])
        stoch = formed - used

        # No need to do further calculations if the reaction does not change
        # the concentration of the molecule.
        if stoch == 0:
            return 0

        if k == None:
            try:
                k = self.reaction_rate(temp, energy_key=energy_key)
            # Raised when self.energies is None
            except TypeError:
                k = self.k

        if k < k_thresh:
            return 0

        symbols = [str(mol_map[educt]) for educt in self.educts]
        j = "*".join(symbols)
        expr = sym.sympify(f"{stoch}*{k}*{j}")
        return expr

    def __str__(self):
        return self.name

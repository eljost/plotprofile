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
        if barrier is None:
            energies = self.energies[energy_key].copy()
            energies -= energies.min()
            energies *= AU2KJMOL * 1000
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

    def get_expr(self, mol, mol_map, temp=298.15, k=None, energy_key="G_solv_alt",
                 k_thresh=1e-8):
        factor = self.factor(mol)

        if factor == 0:
            return 0

        if k == None:
            try:
                k = self.reaction_rate(temp, energy_key=energy_key)
            # Raised when self.energies is None
            except TypeError:
                k = self.k

        # if k < k_thresh:
            # return 0

        # Take stochiometry into account
        look_at = {
            1: self.products,
            -1: self.educts,
        }
        stoch = sum([reactant == mol for reactant in look_at[factor]])
        factor *= stoch

        symbols = [str(mol_map[educt]) for educt in self.educts]
        j = "*".join(symbols)
        expr = sym.sympify(f"{k}*{factor}*{j}")
        return expr


    def __str__(self):
        return self.name

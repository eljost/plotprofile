import itertools as it

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import sympy as sym


def kinetics(reactions, c0s, t_span, temp=298.15):
    # Gather all educts
    mols = set(it.chain(*[rx.educts for rx in reactions]))
    mols = list(mols)

    cs = sym.symbols(f"c:{len(mols)}")
    
    mol_map = {mol: c for mol, c in zip(mols, cs)}

    exprs = list()
    for mol in mols:
        mol_exprs = sum([rx.get_expr(mol, mol_map) for rx in reactions])
        exprs.append(mol_exprs)

    print()
    print("Rate constants")
    max_width = max([len(rx.name) for rx in reactions])
    for i, rx in enumerate(reactions):
        if rx.energies is None:
            continue
        k = rx.reaction_rate(temp=temp)
        print(f"\t{rx.name: >{max_width}}: {k: >10.4e} 1/s")
        if i > 0 and (i-1) % 2 == 0:
            print()

    # print()
    # print("Rate equations")
    # for c_, mol_, exprs_ in zip(cs, mols, exprs):
        # print(mol_, c_)
        # print("\t", exprs_)

    # Create Jacobian
    jac_ = sym.Matrix(exprs).jacobian(cs)
    jac = sym.lambdify(cs, jac_)
    def jac_func(t, cs):
        return jac(*cs)

    # Set starting concentrations
    c0 = np.zeros(len(mols))
    for mol, c0_ in c0s.items():
        try:
            i = mols.index(mol)
            c0[i] = c0_
        except ValueError:
            print(f"Initial concentration for '{mol}' given, but it does not "
                   "participate in any reaction!")

    # Create RHS of ODEs
    rhs = sym.lambdify(cs, exprs)
    def func(t, cs):
        return rhs(*cs)

    ivp_kwargs = {
        # "method": "LSODA",
        # "method": "Radau",
        "method": "BDF",
        "rtol": 1e-8,
        "atol": 1e-11,
    }
    res = solve_ivp(func, t_span, c0, jac=jac_func, **ivp_kwargs)

    return mols, res


def plot_kinetics(mols, res):
    t = res.t
    ys = res.y
    colors = cm.tab20(np.linspace(0, 1, len(mols)))
    fig, ax = plt.subplots()
    print()
    print("Final concentrations [mol/l]")
    for mol, y, color in zip(mols, ys, colors):
        ax.plot(t, y, color=color, label=mol, linewidth=3)
        print(f"\t{mol: >30}: {y[-1]: .12f}")
    # ax.set_ylim(0, 1)
    ax.set_ylim(1e-8, max(ys[:,0]))
    ax.set_yscale("log")
    ax.legend()
    plt.show()

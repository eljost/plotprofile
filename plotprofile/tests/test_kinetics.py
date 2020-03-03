import matplotlib.pyplot as plt
import numpy as np

from plotprofile.constants import KB, H_PLANCK, R
from plotprofile.kinetics import kinetics, plot_kinetics
from plotprofile.Reaction import Reaction


def test_kinetics():
    educts = ["no", "no", "br2"]
    ts = ""
    products = ["nobr", "nobr"]
    kf = 0.42
    kb = 0.17
    energies = None

    frx = Reaction("forward", educts, ts, products, energies, k=kf)
    k = frx.reaction_rate(barrier=25 * 4.184*1000)
    brx = Reaction("forward", products, ts, educts, energies, k=kb)
    c0s = {"no": 1, "br2": 1}
    reactions = (frx, brx)
    t_span = (0, 10)
    mols, res = kinetics(reactions, c0s, t_span)
    # plot_kinetics(mols, res)

    ref = {
        "br2":  0.714712997785,
        "no": 0.429425995571,
        "nobr": 0.570574004429,
    }
    ys = res.y
    ref_ys = [ref[mol] for mol in mols]
    np.testing.assert_allclose(ys[:,-1], ref_ys)


def test_eyring():
    def k(barrier, temp):
        return KB*temp/H_PLANCK * np.exp(-barrier/(R*temp))

    fig, ax = plt.subplots()
    # kcal/mol to J/mol
    kcal_barriers = np.linspace(10, 45, 250)
    joule_barriers = kcal_barriers* 4.184 * 1000
    for temp in (273.15, 298.15, 373.15, 473.15, 573.15):
        ks = k(joule_barriers, temp) * 60 * 60
        ax.plot(kcal_barriers, ks, label=temp)
    ax.legend()
    ax.set_yscale("log")
    ax.set_xlim(10, 45)
    ax.set_ylim(1e-20, 1e15)
    # plt.show()

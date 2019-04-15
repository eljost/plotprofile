#!/usr/bin/env python3

# [1] 10.1021/acs.organomet.8b00456
#     https://pubs.acs.org/doi/10.1021/acs.organomet.8b00456

import argparse
from pprint import pprint
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from plotprofile.parser import parse_orca


AU2KJMOL = 2625.499638 # PT pleaaaase
KB = 1.380649e-23 # J/K, Boltzmann constant
R = 8.3144598 # J/(K*mol), ideal gas constant
H_PLANCK = 6.62607e-34 # Js, Planck constant


def load_molecule_energies(molecules, program):
    parser_funcs = {
        "orca": parse_orca,
    }
    parser = parser_funcs[program]

    energies = dict()
    for mol_name, as_dict in molecules.items():
        print(f"Molecule '{mol_name}'")
        freq_fn = as_dict["freq"]
        scf_fn = as_dict.get("scf", None)
        solv_fn = as_dict.get("solv", None)
        thermo = parser(freq_fn, solv_fn, scf_fn)
        total_energy = thermo.single_point_alt + thermo.g_solv
        energies[mol_name] = total_energy
        print()
    return energies


def to_list(inp):
    if isinstance(inp, str):
        inp = [inp, ]
    elif not isinstance(inp, list):
        inp = list(inp)
    return inp


def make_reactions(reactions, mol_energies):
    def energy_sum(keys):
        return sum([mol_energies[k] for k in keys])

    rx_energies = dict()
    all_rx_energies = list()
    for rx_name, reagents in reactions.items():
        print(f"Calculating energies for reaction '{rx_name}'")
        pprint(reagents)
        educts = to_list(reagents["educts"])
        ts = reagents["ts"]
        products = to_list(reagents["products"])

        educt_en = energy_sum(educts)
        ts_en = mol_energies[ts]
        product_en  = energy_sum(products)
        energies = np.array((educt_en, ts_en, product_en))

        add = reagents.get("add", None)
        if add:
            add_energy = mol_energies[add]
            energies += add_energy
        rx_energies[rx_name] = energies
        all_rx_energies.extend(energies.copy())
        print()
    all_rx_energies = np.array(all_rx_energies)
    return rx_energies, all_rx_energies


def reaction_rate(activation_energy, temperature=298.15):
    """Eyring equation, returns k in units of 1/s."""
    T = temperature
    return KB*T/H_PLANCK * np.exp(-activation_energy/(R*T))


def set_labels(ax, xs, ys, label_strs):
    min_y = min(ys)
    for x, y, lbl in zip(xs, ys, label_strs):
        y_shifted = max(min_y, y-10)
        ax.annotate(lbl, (x, y_shifted), ha="center", fontweight="bold")


def plot_barrier(ax, x, y, barrier):
    props = {
        "arrowstyle": "<->",
    }
    ax.annotate(f"", (x, y+barrier), xytext=(x, y), arrowprops=props)
    ax.text(x+.05, y+(barrier/2), f"{barrier:.1f} kJ/mol")


def plot_rx_energies(rx_energies, all_rx_energies, rx_labels, temperature):
    xs = [0, 1, 2]
    plot_kwargs = {
        "ms": 20,
        "color": "k",
        "marker": "_",
        "linestyle": "--",
    }
    all_label_strs = list()
    for rx_name, energies in rx_energies.items():
        print(rx_name)
        labels = [v for k, v in rx_labels[rx_name].items() if k!= "add"]
        label_strs = [", ".join(to_list(lbl)) for lbl in labels]
        all_label_strs.extend(label_strs)
        ed_lbl, _, prod_lbl = label_strs

        ens = energies - energies.min()
        ens *= AU2KJMOL
        educt, ts, product = ens
        barrier = ts - educt
        print(f"\tbarrier = {barrier:.1f} kJ/mol")
        k = reaction_rate(barrier*1000, temperature=temperature)
        print(f"\tTST rate constant k = {k:.4e} 1/s "
              f"(using T={temperature:.2f} K)")
        fig, ax = plt.subplots()
        ax.plot(xs, ens, **plot_kwargs)
        ax.set_ylabel("$\Delta E / kJ \cdot mol^{-1}$")
        set_labels(ax, xs, ens, label_strs)
        plot_barrier(ax, 0, educt, barrier)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(bottom=False, labelbottom=False)
        fig.suptitle(f"{rx_name}: {ed_lbl} -> {prod_lbl}, k={k:.4e} 1/s")
        pdf_name = f"{rx_name}.pdf"
        fig.savefig(pdf_name)
        print(f"\tsaved PDF to '{pdf_name}'")
        print()
        plt.show()

    if len(rx_energies.keys()) == 1:
        return

    fig, ax = plt.subplots()
    all_rx_energies -= all_rx_energies.min()
    all_rx_energies *= AU2KJMOL
    xs = np.arange(all_rx_energies.size)
    ax.plot(xs, all_rx_energies, **plot_kwargs)
    ax.set_ylabel("$\Delta E / kJ \cdot mol^{-1}$")
    set_labels(ax, xs, all_rx_energies, all_label_strs)
    pdf_name = f"overview.pdf"
    fig.savefig(pdf_name)
    print(f"saved PDF to 'overview.pdf'")
    plt.show()


def dump_energies(rx_energies):
    cols = "educts ts products".split()
    df = pd.DataFrame.from_dict(rx_energies, orient="index", columns=cols)
    csv_fn = "path_energies.csv"
    df.to_csv(csv_fn)
    print(f"Dumped energies to {csv_fn}.")


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("yaml",
        help=".yaml file containing the description of all relevant molecules "
             "and reactions"
    )
    parser.add_argument("-T", default=298.15, type=float,
        help="Temperature to use for TST reaction rate calculation."
    )

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])
    with open(args.yaml) as handle:
        inp_dict = yaml.load(handle)

    mol_energies = load_molecule_energies(inp_dict["molecules"], inp_dict["program"])
    print("Using these G(sol) values:")
    pprint(mol_energies)
    print()

    rx_energies, all_rx_energies = make_reactions(inp_dict["reactions"], mol_energies)
    print("Calculated reaction energies:")
    pprint(rx_energies)
    print()

    dump_energies(rx_energies)
    plot_rx_energies(rx_energies, all_rx_energies, inp_dict["reactions"], args.T)


if __name__ == "__main__":
    run()

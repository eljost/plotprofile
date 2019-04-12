#!/usr/bin/env python3

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import yaml

from plotprofile.parser import parse_orca

AU2KJMOL = 2625.50


def load_molecule_energies(molecules, program):
    parser_funcs = {
        "orca": parse_orca,
    }
    parser = parser_funcs[program]

    energies = dict()
    for mol_name, as_dict in molecules.items():
        freq_fn = as_dict["freq"]
        scf_fn = as_dict.get("scf", None)
        solv_fn = as_dict.get("solv", None)
        thermo = parser(freq_fn, solv_fn, scf_fn)
        total_energy = thermo.single_point_alt + thermo.g_solv
        energies[mol_name] = total_energy
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
        print(rx_name, reagents)
        educts = to_list(reagents["educts"])
        ts = reagents["ts"]
        products = to_list(reagents["products"])

        # import pdb; pdb.set_trace()
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
    all_rx_energies = np.array(all_rx_energies)
    return rx_energies, all_rx_energies


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


def plot_rx_energies(rx_energies, all_rx_energies, rx_labels):
    xs = [0, 1, 2]
    plot_kwargs = {
        "ms": 20,
        "color": "k",
        "marker": "_",
        "linestyle": "--",
    }
    all_label_strs = list()
    for rx_name, energies in rx_energies.items():
        labels = [v for k, v in rx_labels[rx_name].items() if k!= "add"]
        label_strs = [", ".join(to_list(lbl)) for lbl in labels]
        all_label_strs.extend(label_strs)

        ens = energies - energies.min()
        ens *= AU2KJMOL
        educt, ts, product = ens
        barrier = ts - educt
        print(f"{rx_name}: barrier = {barrier:.1f} kJ/mol")
        fig, ax = plt.subplots()
        ax.plot(xs, ens, **plot_kwargs)
        ax.set_ylabel("$\Delta E / kJ \cdot mol^{-1}$")
        set_labels(ax, xs, ens, label_strs)
        plot_barrier(ax, 0, educt, barrier)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(bottom=False, labelbottom=False)
        fig.suptitle(rx_name)
        plt.show()

    fig, ax = plt.subplots()
    all_rx_energies -= all_rx_energies.min()
    all_rx_energies *= AU2KJMOL
    xs = np.arange(all_rx_energies.size)
    ax.plot(xs, all_rx_energies, **plot_kwargs)
    set_labels(ax, xs, all_rx_energies, all_label_strs)
    plt.show()


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("yaml")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])
    with open(args.yaml) as handle:
        inp_dict = yaml.load(handle)

    mol_energies = load_molecule_energies(inp_dict["molecules"], inp_dict["program"])
    from pprint import pprint
    pprint(mol_energies)

    rx_energies, all_rx_energies = make_reactions(inp_dict["reactions"], mol_energies)
    pprint(rx_energies)

    plot_rx_energies(rx_energies, all_rx_energies, inp_dict["reactions"])


if __name__ == "__main__":
    run()

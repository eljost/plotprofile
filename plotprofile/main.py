#!/usr/bin/env python3

# [1] 10.1021/acs.organomet.8b00456
#     https://pubs.acs.org/doi/10.1021/acs.organomet.8b00456

import argparse
from collections import namedtuple
import itertools as it
from pathlib import Path
from pprint import pprint
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from plotprofile.parser import parse_orca
from thermoanalysis.QCData import QCData
from thermoanalysis.thermo import thermochemistry


AU2KJMOL = 2625.499638 # PT pleaaaase
KB = 1.380649e-23 # J/K, Boltzmann constant
R = 8.3144598 # J/(K*mol), ideal gas constant
H_PLANCK = 6.62607e-34 # Js, Planck constant
CAL2J = 4.1868 # Calories to J

GibbsEnergies = namedtuple("GibbsEnergies",
                           "G_gas G_solv G_gas_alt G_solv_alt"
)

GFIELDS = GibbsEnergies._fields


def load_thermos(molecules, program):
    parser_funcs = {
        "orca": parse_orca,
    }
    parser = parser_funcs[program]

    energies = dict()
    thermos = dict()
    for mol_name, as_dict in molecules.items():
        base_path = Path(as_dict.get("base", "./"))
        # Frequency calculation is mandatory, alternative single point (scf)
        # and solvated calculation are optional.
        freq_fn = base_path / as_dict["freq"]
        scf_fn = as_dict.get("scf", None)
        scf_fn = base_path / scf_fn if scf_fn else None
        solv_fn = as_dict.get("solv", None)
        solv_fn = base_path / solv_fn if solv_fn else None
        thermo = parser(freq_fn, solv_fn, scf_fn)
        thermos[mol_name] = thermo
    return thermos


def load_molecule_energies(thermos, no_alt):
    energies = dict()
    for mol_name, thermo in thermos.items():
        G_gas = thermo.single_point + thermo.gibbs_corr
        G_solv = G_gas + thermo.dG_solv
        G_gas_alt = thermo.single_point_alt + thermo.gibbs_corr
        G_solv_alt = G_gas_alt + thermo.dG_solv

        gibbs_energies = GibbsEnergies(
                            G_gas=G_gas,
                            G_solv=G_solv,
                            G_gas_alt=G_gas_alt,
                            G_solv_alt=G_solv_alt,
        )
        energies[mol_name] = gibbs_energies
        print(f"Loaded molecule '{mol_name}'")
    return energies


def to_list(inp):
    if isinstance(inp, str):
        inp = [inp, ]
    elif not isinstance(inp, list):
        inp = list(inp)
    return inp


def make_reactions(reactions, mol_energies):

    def energy_sum(keys):
        return {
            field: sum([getattr(mol_energies[k], field) for k in keys])
            for field in GFIELDS
        }

    rx_energies = dict()
    all_rx_energies = list()
    for rx_name, reagents in reactions.items():
        print(f"Calculating energies for reaction '{rx_name}'")
        pprint(reagents)
        educts = to_list(reagents["educts"])
        ts = [reagents["ts"], ]
        assert len(ts) == 1
        products = to_list(reagents["products"])
        add = reagents.get("add", [])
        add = [add, ] if isinstance(add, str) else list(add)

        educt_ens = energy_sum(educts + add)
        ts_ens = energy_sum(ts + add)
        product_ens  = energy_sum(products + add)

        energies = {}
        for field in GFIELDS:
            ens = np.array(
                [dct[field] for dct in (educt_ens, ts_ens, product_ens)]
            )
            energies[field] = ens

        rx_energies[rx_name] = energies
    return rx_energies


ForDash = namedtuple("ForDash", "single_point single_point_alt dG_solv dG")


def make_dash_data(reactions, thermos, thermochem_df):

    unique_Ts = thermochem_df["T"].unique()

    def sum_(keys, field):
        return sum([getattr(thermos[k], field) for k in keys])

    df_cols = "name reactant single_point single_point_alt dG_solv T dG".split()
    rows = list()
    for rx_name, reagents in reactions.items():
        add = reagents.get("add", [])
        add = [add, ] if isinstance(add, str) else list(add)
        for reactant in ("educts", "ts", "products"):
            reactants = to_list(reagents[reactant])
            if isinstance(reactants, str):
                reactants = [reactants, ]
            reactants = reactants + add
            sp = sum_(reactants, "single_point")

            reactant_thermo = thermochem_df[thermochem_df["name"].isin(reactants)]
            dGs = reactant_thermo.groupby("T").dG.sum().to_list()

            sp_alt = sum_(reactants, "single_point_alt")
            dG_solv = sum_(reactants, "dG_solv")
            for T, dG in zip(unique_Ts, dGs):
                row = [rx_name, reactant, sp, sp_alt, dG_solv, T, dG]
                rows.append(row)
    df = pd.DataFrame(rows, columns=df_cols)
    return df


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


def plot_rx_energies(rx_energies, rx_labels, temperature, plot_kwargs):
    print("Reactions")
    for rx_name, energies in rx_energies.items():
        labels = rx_labels[rx_name]
        plot_rx(rx_name, energies, labels, temperature, plot_kwargs)

def plot_rx(rx_name, energies, labels, temperature, plot_kwargs):
    xs = [0, 1, 2]
    ed_lbl, _, prod_lbl = labels

    ens = energies - energies.min()
    ens *= AU2KJMOL
    educt, ts, product = ens
    barrier = ts - educt
    k = reaction_rate(barrier*1000, temperature=temperature)
    fig, ax = plt.subplots()
    ax.plot(xs, ens, **plot_kwargs)
    ax.set_ylabel("$\Delta E / kJ \cdot mol^{-1}$")
    set_labels(ax, xs, ens, labels)
    plot_barrier(ax, 0, educt, barrier)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)
    fig.suptitle(f"{rx_name}: {ed_lbl} -> {prod_lbl}, k={k:.4e} 1/s")
    pdf_name = f"{rx_name}.pdf"
    fig.savefig(pdf_name)
    # print(f"\t{rx_name}: Saved PDF to '{pdf_name}'")
    print(f"\tSaved PDF to '{pdf_name}'")
    plt.show()


def plot_paths(rx_energies, paths, rx_labels, plot_kwargs):
    print("Paths")
    for path_name, rx_names in paths.items():
        path_energies = list()
        path_labels = list()
        for rx_name in rx_names:
            path_energies.extend(rx_energies[rx_name])
            path_labels.extend(rx_labels[rx_name])
        path_energies = np.array(path_energies)
        plot_path(path_energies, path_name, rx_names, path_labels, plot_kwargs)


def plot_path(path_energies, path_name, rx_names, path_labels, plot_kwargs):
    start_energy = path_energies[0]

    fig, ax = plt.subplots()
    path_energies -= path_energies[0]
    path_energies *= AU2KJMOL
    xs = np.arange(path_energies.size)
    ax.plot(xs, path_energies, **plot_kwargs)
    ax.set_title(path_name)
    ax.set_ylabel("$\Delta E / kJ \cdot mol^{-1}$")
    set_labels(ax, xs, path_energies, path_labels)
    path_fn = f"path_{path_name}.pdf"
    fig.savefig(path_fn)
    print(f"\tSaved PDF to '{path_fn}'")
    plt.show()


def dump_energies(rx_energies):
    cols = "educts ts products".split()
    col_prods = ["_".join(_) for _ in it.product(GFIELDS, cols)]
    # Flatten dict
    rx_energies_flat = {
        rx_name:  list(it.chain(*[ens[f] for f in GFIELDS]))
	for rx_name, ens in rx_energies.items()
    }
    df = pd.DataFrame.from_dict(rx_energies_flat, orient="index", columns=col_prods)
    csv_fn = "path_energies.csv"
    df.to_csv(csv_fn)
    print(f"Dumped energies to {csv_fn}.")


def print_rx_energies(rx_energies, temperature):
    print("Calculated reaction energies:")
    for rx_name, rx_ens in rx_energies.items():
        print(f"{rx_name}:")
        for field in GFIELDS:
            ens = rx_ens[field].copy()
            # ens -= ens.min()
            ens -= ens[0]
            ens *= AU2KJMOL
            ens_str = np.array2string(ens, precision=1)
            barrier = ens[1] - ens[0]
            print(f"\t{field: >10s}: {ens_str: >24} kJ/mol, barrier: {barrier: >6.1f} kJ/mol")

        print()
        kcal_barrier = barrier / CAL2J
        print(f"\t(best) barrier = {barrier:.1f} kJ/mol ({kcal_barrier:.1f} kcal/mol)")
        k = reaction_rate(barrier*1000, temperature=temperature)
        k_day = k * 3600 * 24
        print(f"\tTST rate constant k = {k:.4e} 1/s ({k_day:.4e} 1/d) "
              f"@ T={temperature:.2f} K")
        print()


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("yaml",
        help=".yaml file containing the description of all relevant molecules "
             "and reactions"
    )
    parser.add_argument("-T", default=298.15, type=float,
        help="Temperature to use for TST reaction rate calculation."
    )
    parser.add_argument("--parsedash", action="store_true")
    parser.add_argument("--no_alt", action="store_true")
    parser.add_argument("--norxs", action="store_true")
    parser.add_argument("--nopaths", action="store_true")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    no_alt = args.no_alt
    plot_rxs = not args.norxs
    plot_paths = not args.nopaths
    temperature = args.T

    with open(args.yaml) as handle:
        inp_dict = yaml.load(handle)

    thermos = load_thermos(inp_dict["molecules"], inp_dict["program"])
    mol_energies = load_molecule_energies(thermos, no_alt=no_alt)

    if args.parsedash:
        freq_logs = {k: v["freq"] for k, v in inp_dict["molecules"].items()}
        qc_datas = {k: QCData(v) for k, v in freq_logs.items()}
        temps = np.arange(298.00, 373.00, 5)
        rows = list()
        for mol_name, temp in it.product(qc_datas, temps):
            qc = qc_datas[mol_name]
            tc = thermochemistry(qc, temp)
            row = [mol_name, ] + list(tc)
            rows.append(row)
        cols = ["name", ] + list(tc._fields)
        df = pd.DataFrame(rows, columns=cols)
        df.to_pickle("thermochem_data")
        thermochem_df = pd.read_pickle("thermochem_data")
        df = make_dash_data(inp_dict["reactions"], thermos, thermochem_df)
        df.to_pickle("dash_data")

    print("Using these G-values:")
    for key, val in mol_energies.items():
        print(key)
        pprint(val._asdict())
    print()

    rx_energies = make_reactions(inp_dict["reactions"], mol_energies)
    print_rx_energies(rx_energies, temperature)
    print()

    # dump_energies(rx_energies)

    rx_labels = {}
    rxs = inp_dict["reactions"]
    # Create label strings for plotting
    for rx_name in rxs:
        labels = [v for k, v in rxs[rx_name].items() if k!= "add"]
        label_strs = [", ".join(to_list(lbl)) for lbl in labels]
        rx_labels[rx_name] = label_strs

    # Try to use the 'best' energies for plotting. That is with alternative
    # single point and solvation.
    best_rx_energies = {
        rx_name: rx_energies[rx_name]["G_solv_alt"] for rx_name in rx_labels
    }

    plot_kwargs = {
        "ms": 20,
        "color": "k",
        "marker": "_",
        "linestyle": "--",
    }
    if plot_rxs:
        plot_rx_energies(best_rx_energies, rx_labels, temperature, plot_kwargs)
    else:
        print("Skipped plotting of reactions!")

    if plot_paths:
        paths = inp_dict.get("paths", None)
        if paths is None:
            print("Found no defined paths in .yaml file. Using all defined "
                  "reactions instead.")
            paths = {"reaction path": list(rxs.keys())}

        plot_paths(best_rx_energies, paths, rx_labels, plot_kwargs)
    else:
        print("Skipped plotting of reaction paths!")


if __name__ == "__main__":
    run()

#!/usr/bin/env python3

# [1] 10.1021/acs.organomet.8b00456
#     https://pubs.acs.org/doi/10.1021/acs.organomet.8b00456

import argparse
from collections import namedtuple
import itertools as it
import os
from pathlib import Path
from pprint import pprint
import sys
import textwrap

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

PLOT_KWARGS = {
    "ms": 30,
    "mew": 2,
    "color": "k",
    "marker": "_",
    "ls": "--",
}


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


def set_labels(ax, xs, ys, label_strs, y_shift=10, ts_above=False,
               fused=False):
    min_y = min(ys)
    max_y = max(ys)
    ts_step = 2 if fused else 3
    ts_inds = range(1, len(label_strs), ts_step)
    for i, (x, y, lbl) in enumerate(zip(xs, ys, label_strs)):
        # Print labels below the markes so we decrease the y position a bit
        y_shifted = max(min_y, y-y_shift)
        # Put TS label above marker, instead of below
        if ts_above and (i in ts_inds):
            y_shifted = min(max_y, y+.75*y_shift)
        ax.annotate(lbl, (x, y_shifted), ha="center", fontweight="bold")


def savefig(fig, base_name, out_dir=None, save_as=(".pdf", ".svg")):
    if out_dir is None:
        out_dir = "."
    out_path = Path(out_dir)
    try:
        os.mkdir(out_path)
    except FileExistsError:
        pass

    names = [base_name + ext for ext in save_as]
    full_paths = [out_path / name for name in names]
    [fig.savefig(p) for p in full_paths]
    print(f"\tSaved {', '.join(names)}")
    return full_paths


def subplots():
    fig, ax = plt.subplots()
    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    return fig, ax


def expand_points(ax, xs, ys, ms):
    yy = np.repeat(ys, 2)
    xx = np.arange(len(ys))[:,None]

    trans = ax.transData.transform
    trans_inv = ax.transData.inverted().transform
    ms_half = ms / 2
    # Transform to display coordinates
    xt, yt = trans((0, 0))
    xt_left = xt - ms_half
    xt_right = xt + ms_half
    # Transform back to data coordinates
    x_left, _ = trans_inv((xt_left, yt))
    x_right, _ = trans_inv((xt_right, yt))
    dx = np.array((x_left, x_right))[None,:] * 1.325
    xx = xx + dx
    xx = xx.flatten()
    return xx, yy


def plot(ax, xs, ys):
    marker_kwargs = PLOT_KWARGS.copy()
    marker_update = {
        "ls": "",
        "zorder": 10,
    }
    marker_kwargs.update(marker_update)

    line_kwargs = PLOT_KWARGS.copy()
    line_update = {
        "ls": "--",
        "marker": "",
        "zorder": 1,
    }
    line_kwargs.update(line_update)

    # Plot only marker
    ax.plot(xs, ys, **marker_kwargs)

    # Plot lines. Expand points for this
    xx, yy = expand_points(ax, xs, ys, marker_kwargs["ms"])
    ax.plot(xx, yy, **line_kwargs)


def plot_barrier(ax, x, y, barrier):
    props = {
        "arrowstyle": "<->",
    }
    ax.annotate(f"", (x, y+barrier), xytext=(x, y), arrowprops=props)
    ax.text(x+.05, y+(barrier/2), f"{barrier:.1f} kJ/mol")


def plot_reactions(rx_energies, rx_labels, rx_titles, temperature):
    print("Reactions")
    for rx_name, energies in rx_energies.items():
        labels = rx_labels[rx_name]
        title = rx_titles[rx_name]
        plot_reaction(rx_name, energies, labels, title, temperature)


def plot_reaction(rx_name, energies, labels, title, temperature):
    xs = [0, 1, 2]
    ed_lbl, _, prod_lbl = labels

    ens = energies - energies.min()
    ens *= AU2KJMOL
    educt, ts, product = ens
    barrier = ts - educt
    k = reaction_rate(barrier*1000, temperature=temperature)
    fig, ax = subplots()
    plot(ax, xs, ens)
    ax.set_ylabel("$\Delta E / kJ \cdot mol^{-1}$")
    set_labels(ax, xs, ens, labels, ts_above=True)
    plot_barrier(ax, 0, educt, barrier)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)
    ax.set_title(f"{title}, k={k:.4e} 1/s")
    savefig(fig, rx_name, out_dir="reactions")

    # plt.show()
    plt.close()


def plot_paths(rx_energies, paths, rx_labels, rx_titles):
    print("Paths")
    for path_name, rx_names in paths.items():
        path_energies = list()
        path_labels = list()
        barriers = list()
        titles = list()
        for rx_name in rx_names:
            educts, ts, _ = rx_energies[rx_name]
            barrier = ts - educts
            barriers.append(barrier)
            path_energies.extend(rx_energies[rx_name])
            path_labels.extend(rx_labels[rx_name])
            titles.append(rx_titles[rx_name])
        path_energies = np.array(path_energies)
        barriers = np.array(barriers)
        plot_path(path_energies, path_name, barriers, rx_names, titles, path_labels)
        plot_path(path_energies, path_name, barriers, rx_names, titles, path_labels,
                  fuse=True)


def plot_path(path_energies, path_name, barriers, rx_names, rx_titles,
              path_labels, fuse=False):
    path_energies = path_energies.copy()
    path_labels = path_labels.copy()

    if fuse:
        assert len(path_energies) % 3 == 0
        # Keep the last entry, don't drop it
        drop_inds = range(2, len(path_energies), 3)[:-1]
        keep_inds = [i for i, _ in enumerate(path_energies)
                     if i not in drop_inds]
        path_energies = path_energies[keep_inds]
        path_labels = [path_labels[i] for i in keep_inds]

    start_energy = path_energies[0]

    fig, ax = subplots()
    path_energies -= path_energies[0]
    path_energies *= AU2KJMOL
    xs = np.arange(path_energies.size)
    plot(ax, xs, path_energies)
    ax.set_title(path_name)
    ax.set_ylabel("$\Delta E / kJ \cdot mol^{-1}$")
    set_labels(ax, xs, path_energies, path_labels, y_shift=35, ts_above=True,
               fused=fuse)

    # Here also markup characters like $$ etc. are counted, so this may not
    # be very accurate...
    max_width = max([len(title) for title in rx_titles])
    pad_titles = [" " * (max_width - len(title)) + title for title in rx_titles]
    barriers = barriers.copy()
    barriers *= AU2KJMOL
    barrier_str = "Barriers in kJ mol⁻¹:\n\n" \
        + "\n".join([
            f"{rx_title}: {barrier: >6.1f}" for rx_title, barrier in zip(pad_titles,
                                                                         barriers)
        ])
    ax.annotate(barrier_str, (0, -175), family="monospace")
    fused_str = "_fused" if fuse else ""
    base_name = path_name + fused_str
    out_dir = "paths" + fused_str
    savefig(fig, base_name, out_dir=out_dir)

    plt.close()
    # plt.show()


def get_reactants(rx):
    reactants = list()
    for key in ("educts", "ts", "products"):
        reactants.append(to_list(rx[key]))
    # 'add' may not be present so we have to handle it separately
    reactants.append(to_list(rx.get("add", [])))
    return reactants


def rx_to_string(rx):
    print(rx)
    educts, ts, products, add = get_reactants(rx)

    eds = " + ".join(educts)
    prods = " + ".join(products)
    ts = ts[0]

    # Don't print the 'add' species as they don't participate in the reaction.
    # add = " "*10 + f"+({', '.join(add)})" if add else ""
    # rx_str = f"{eds} ---> [{ts}]‡ ---> {prods}" + add
    # rx_str = f"{eds} ---> [{ts}]‡ ---> {prods}"

    add = f"\nConsidering additionally:  {', '.join(add)}" if add else ""
    # sep = "\n" + " "*(len(eds) // 2) + "↓" + "\n"
    sep = " "*(len(eds) // 2) + "↓"
    rx_str = f"{eds}\n{sep}\n[{ts}]‡\n{sep}\n{prods}{add}"
    return rx_str


def print_path_rx_energies(paths, rx_energies, rx_strs, temperature):
    for path, rxs in paths.items():
        print(f"# Reaction energies for {path}")
        print()
        path_energies = {rx_name: rx_energies[rx_name] for rx_name in rxs}
        print_rx_energies(path_energies, rx_strs, temperature)
        print()


def print_rx_energies(rx_energies, rx_strs, temperature):
    for rx_name, rx_ens in rx_energies.items():
        print(f"\t## {rx_name.capitalize()}:")
        print(textwrap.indent(rx_strs[rx_name], "\t"))
        print()
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
        print(f"\tG_solv_alt barrier = {barrier:.1f} kJ/mol ({kcal_barrier:.1f} kcal/mol)")
        k = reaction_rate(barrier*1000, temperature=temperature)
        k_day = k * 3600 * 24
        print(f"\tTST rate constant k = {k:.4e} 1/s ({k_day:.4e} 1/d) "
              f"@ T={temperature:.2f} K")
        print()


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


def plot_compare(name, molecules, energies, ylabel):
    fig, ax = subplots()
    xs = [i for i, _ in enumerate(molecules)]
    plot(ax, xs, energies)
    set_labels(ax, xs, energies, molecules, y_shift=1)
    ax.set_ylabel(f"$\Delta${ylabel} / kJ mol⁻¹")
    ax.set_title(f"Comparison: '{name}'")
    return fig, ax


def compare_molecules(to_compare, mol_energies, attr="G_solv_alt"):
    for name, molecules in to_compare.items():
        print(f"Comparison of '{name}' using '{attr}'")
        energies = np.array([getattr(mol_energies[mol], attr) for mol in molecules])
        energies -= energies.min()
        energies *= AU2KJMOL
        inds = np.argsort(energies)

        most_stable_ind = inds[0]
        print(f"\t'{molecules[most_stable_ind]}' has the lowest energy.")
        # Skip first entry as this was already handled above (most_stable)
        for ind in inds[1:]:
            print(f"\t'{molecules[ind]}' is {energies[ind]:.2f} higher "
                   "kJ mol⁻¹ in energy.")

        fig, ax = plot_compare(name, molecules, energies, ylabel=attr)
        # Only use attr in name if it is not the default
        attr_str = "" if attr == "G_solv_alt" else f"_{attr}"
        base_name = f"{name}{attr_str}_comparison"
        savefig(fig, base_name, out_dir="comparisons")

        plt.close()
        # plt.show()
        print()


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("yaml", nargs="?",
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
    parser.add_argument("-q", "--quick",
        help="Quick plotting of one profile. Educts, TS and products are separated "
            "by ';'. Multiple educts and/or products are separated by ','. Example: "
            "plotprofile --quick 'ed.out;ts.out;prod.out'.")
    parser.add_argument("-i", "--interactive", action="store_true")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    no_alt = args.no_alt
    show_rxs = not args.norxs
    show_paths = not args.nopaths
    temperature = args.T

    if args.quick:
        show_paths = False
        educts, ts, products = [reactants.split(",")
                                for reactants in args.quick.split(";")]
        assert len(ts) == 1
        ed_keys = [f"ed{i}" for i in range(len(educts))]
        prod_keys = [f"prod{i}" for i in range(len(products))]

        inp_dict = {
            "molecules": dict(),
            "program": "orca",
        }
        def add_mols(keys, fns):
            for key, fn in zip(keys, fns):
                inp_dict["molecules"][key] = {
                    "freq": fn,
                }

        add_mols(ed_keys, educts)
        add_mols(["ts", ], ts)
        add_mols(prod_keys, products)

        # Dummy reaction
        inp_dict["reactions"] = {
            "quick": {
                "educts": ed_keys,
                "ts": "ts",
                "products": prod_keys,
            }
        }
    elif args.yaml:
        with open(args.yaml) as handle:
            inp_dict = yaml.load(handle)
    else:
        print("Specify either a YAML file or use --quick! See also 'plotprofile --help'")
        sys.exit()

    reactions = inp_dict["reactions"]

    if args.interactive:
        rx_keys = list(reactions.keys())
        for i, rx in enumerate(rx_keys):
            print(f"{i:02d}: {rx}")
        ind = int(input("Show reaction: ", ))
        rx_key = rx_keys[ind]
        reaction = reactions[rx_key]
        reactions = {rx_key: reaction, }
        # Modify molecules so only the ones actually needed are loaded
        rx_molecules = list(it.chain(*[to_list(mols) for mols in reaction.values()]))
        inp_dict["molecules"] = {key: val for key, val in inp_dict["molecules"].items()
                                 if key in rx_molecules}

    molecules = inp_dict["molecules"]
    thermos = load_thermos(molecules, inp_dict["program"])
    mol_energies = load_molecule_energies(thermos, no_alt=no_alt)

    rx_strs = {rx_name: rx_to_string(rx)
               for rx_name, rx in reactions.items()
    }

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

    print()
    print("Using these G-values:")
    for key, val in mol_energies.items():
        print(key)
        pprint(val._asdict())
    print()

    rx_energies = make_reactions(reactions, mol_energies)

    paths = inp_dict.get("paths", None)
    if paths is None:
        print("Found no defined paths in .yaml file. Using all defined "
              "reactions instead.")
        paths = {"autogenerated": list(reactions.keys())}
    if args.interactive:
        paths = dict()

    print_path_rx_energies(paths, rx_energies, rx_strs, temperature)
    dump_energies(rx_energies)

    mol_labels = {
        # Use molecule key as fallback if no label is defined
        mol: values.get("label", mol) for mol, values in molecules.items()
    }
    rx_labels = {}
    rx_titles = {}
    # Create label strings for plotting
    for rx_name in reactions:
        labels = [v for k, v in reactions[rx_name].items()
                  if k in ("educts", "ts", "products")]

        pretty_labels = list()
        for lbl in labels:
            if isinstance(lbl, str):
                lbl = [lbl, ]
            # Replace with pretty labels
            lbl = [mol_labels[mol] for mol in lbl]
            lbl = ",\n".join(lbl)
            pretty_labels.append(lbl)
        rx_labels[rx_name] = pretty_labels

        rx_title = reactions[rx_name].get("label", rx_name)
        rx_titles[rx_name] = rx_title

    # Try to use the 'best' energies for plotting. That is with alternative
    # single point and solvation.
    best_rx_energies = {
        rx_name: rx_energies[rx_name]["G_solv_alt"] for rx_name in rx_labels
    }

    if show_rxs:
        plot_reactions(best_rx_energies, rx_labels, rx_titles, temperature)
    else:
        print("Skipped plotting of reactions!")

    if show_paths and not args.interactive and (len(best_rx_energies) > 1):
        plot_paths(best_rx_energies, paths, rx_labels, rx_titles)
    else:
        print("Skipped plotting of reaction paths!")

    path_reactions = set(list(it.chain(*paths.values())))
    all_reactions = set(list(reactions.keys()))
    remainder = all_reactions - path_reactions
    if len(remainder) > 0:
        print()
        print("Warning!".upper())
        print("Not all defined reactions appeared in (defined) paths!")
        for i, rx in enumerate(remainder):
            print(f"\t{i:02d}: {rx}")
        print("Warning!".upper())
        print()

        remainder_path = {"remainder": remainder, }
        print_path_rx_energies(remainder_path, rx_energies, rx_strs, temperature)

    to_compare = inp_dict.get("compare", dict())
    compare_molecules(to_compare, mol_energies)
    # compare_molecules(to_compare, mol_energies, attr="G_gas_alt")


if __name__ == "__main__":
    run()

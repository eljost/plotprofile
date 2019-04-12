#!/usr/bin/env python3

from collections import namedtuple
import re

FLOAT_RE = "([-\d\.]+)"
ThermoTuple = namedtuple("ThermoTuple",
                         "single_point gibbs"
)
ThermoSolvTuple = namedtuple("ThermoSolvTuple",
                             "single_point gibbs g_solv single_point_alt"
)


def parse_g_solv(text):
    """Returns Corrected G(solv) in Hartree."""
    g_solv_re = re.compile("Corrected G\(solv\)\s+:\s+" + FLOAT_RE + " Eh")
    # mobj = g_solv_re.search(text)
    # g_solv_tot = float(mobj[1])
    all_g_solv_tots = g_solv_re.findall(text)
    print(f"Found {len(all_g_solv_tots)} solvated calculation(s). Using the last one.")
    g_solv_tot = float(all_g_solv_tots[-1])
    return g_solv_tot


def parse_thermo(text):
    single_point_re = re.compile("FINAL SINGLE POINT ENERGY\s+" + FLOAT_RE)
    all_ens = single_point_re.findall(text)
    print(f"Found {len(all_ens)} single point energies. Using the last one.")
    single_point = float(all_ens[-1])

    gibbs_re = re.compile("Final Gibbs free enthalpy\s+\.\.\.\s+" + FLOAT_RE + " Eh")
    all_gibbs = gibbs_re.findall(text)
    print(f"Found {len(all_gibbs)} thermochemistry calculation(s). Using the last one.")
    gibbs = float(all_gibbs[-1])

    thermo = ThermoTuple(
                single_point=single_point,
                gibbs=gibbs,
    )
    return thermo


def parse_orca(freq_fn, solv_fn=None, scf_fn=None):
    print(f"Reading thermochemistry data from '{freq_fn}'.")
    with open(freq_fn) as handle:
        freq_text = handle.read()
    thermo = parse_thermo(freq_text)

    single_point_alt = thermo.single_point
    g_solv = 0.0
    if solv_fn:
        print(f"Reading solvation data from '{solv_fn}'.")
        with open(solv_fn) as handle:
            solv_text = handle.read()
        g_solv_tot = parse_g_solv(solv_text)
        g_solv = g_solv_tot - thermo.single_point
    if scf_fn:
        raise Exception("Ignoring this right now!")

    thermo_solv = ThermoSolvTuple(
                    *thermo,
                    g_solv=g_solv,
                    single_point_alt=single_point_alt,
    )
    return thermo_solv


def run():
    freq = "01_alkinaddition_alkin_edukt_opt_freq.out"
    solv = "02_alkinaddition_alkin_edukt_smd_dioxan.out"
    thermo = parse_orca(freq, solv)
    print(thermo)

if __name__ == "__main__":
    run()

#!/usr/bin/env python3

from collections import namedtuple
import re

import numpy as np

FLOAT_RE = "([-\d\.]+)"
ThermoTuple = namedtuple("ThermoTuple",
                         "single_point gibbs gibbs_corr"
)
ThermoSolvTuple = namedtuple("ThermoSolvTuple",
                             "single_point gibbs gibbs_corr dG_solv single_point_alt"
)


# def parse_g_solv(text):
    # """Returns Corrected G(solv) in Hartree."""
    # g_solv_re = re.compile("Corrected G\(solv\)\s+:\s+" + FLOAT_RE + " Eh")
    # all_g_solv_tots = g_solv_re.findall(text)
    # print(f"\t\tFound {len(all_g_solv_tots)} solvated calculation(s). Using the last one.")
    # g_solv_tot = float(all_g_solv_tots[-1])
    # return g_solv_tot


def parse_final_single_point(text):
    single_point_re = re.compile("FINAL SINGLE POINT ENERGY\s+" + FLOAT_RE)
    all_ens = single_point_re.findall(text)
    print(f"\t\tFound {len(all_ens)} final single point energies. Using the last one.")
    single_point = float(all_ens[-1])
    return single_point


def parse_thermo(text):
    single_point = parse_final_single_point(text)
    gibbs_re = re.compile("Final Gibbs free (?:enthalpy|energy)\s+\.\.\.\s+" + FLOAT_RE + " Eh")
    all_gibbs = gibbs_re.findall(text)
    print(f"\t\tFound {len(all_gibbs)} thermochemistry calculation(s). Using the last one.")
    gibbs = float(all_gibbs[-1])
    gibbs_corr_re = re.compile("G\-E\(el\)\s+\.\.\.\s+" + FLOAT_RE + " Eh")
    gibbs_corr_refs = gibbs_corr_re.findall(text)
    gibbs_corr_ref = float(gibbs_corr_refs[-1])
    gibbs_corr = gibbs - single_point
    np.testing.assert_allclose(gibbs_corr, gibbs_corr_ref, atol=1e-8, rtol=1e-4)

    thermo = ThermoTuple(
                single_point=single_point,
                gibbs=gibbs,
                gibbs_corr=gibbs_corr,
    )
    return thermo


def parse_orca(freq_fn, solv_fn=None, scf_fn=None):
    print(f"\tReading thermochemistry data from '{freq_fn}'.")
    with open(freq_fn) as handle:
        freq_text = handle.read()
    thermo = parse_thermo(freq_text)

    single_point_alt = thermo.single_point
    dG_solv = 0.0
    if solv_fn:
        print(f"\tReading solvation data from '{solv_fn}'.")
        with open(solv_fn) as handle:
            solv_text = handle.read()
        solvated_single_point = parse_final_single_point(solv_text)
        dG_solv = solvated_single_point - thermo.single_point
    if scf_fn:
        print(f"\tReading alternative single point energy from '{scf_fn}'")
        with open(scf_fn) as handle:
            single_point_alt_text = handle.read()
        single_point_alt = parse_final_single_point(single_point_alt_text)
        print(f"\t\tAlternative single point energy: {single_point_alt:.8f} au")

    thermo_solv = ThermoSolvTuple(
                    *thermo,
                    dG_solv=dG_solv,
                    single_point_alt=single_point_alt,
    )
    return thermo_solv


def run():
    thermo = parse_orca(freq, solv)
    print(thermo)

if __name__ == "__main__":
    run()

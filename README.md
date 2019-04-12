# plotprofile
Quickly generate plots of solution-phase free energies (G(sol)) for reaction paths.

## Installation
```
git clone https://github.com/eljost/plotprofile.git
cd plotprofile
python setup.py develop
```

## Usage
```
plotprofile [yaml]
```

Calculated data is dumped to a CSV file.


## Overview idea

All relevant energy values are directly parsed from logfiles. In this regard no manual
work is required. The necessary program input is given in a YAML file. An example can
be found in the .plotprofile/tests subfolder.

Three sections are required in the YAML file: `program`, `molecules`, and `reactions`.

### program
Right now only `orca` is supported.

### molecules
Three different logfiles can be provided: freq (thermochemistry data), solv (solvated calculation) and scf (alternative single point energy, e.g. at a higher level/bigger basis set etc.)

```
molecules:
 hcn:
  freq: 05_hcn_iso_educt_freq.out
  solv: 06_hcn_iso_educt_solv.out
  scf: 10_hcn_iso_educt_dlpnoccsdt.out
```

### reactions
Every reaction is comprised of three parts (required): educts, ts, and products. The key `add` can
be used to add a constant energy by a molecule to all reagents of a reaction. This can be useful
if e.g. a water molecule is not present in the first reaction step but introduced later on.

```
reactions:
 alkin_addition:
  educts: [alkin, pbu3]
  ts: ts1
  products: addukt
  add: water
 water_addition:
  educts: [addukt, water]
  ts: ts2
  products: addukt_h_oh
 ...
```


G(sol) is calculated as described in [10.1021/acs.organomet.8b00456](https://pubs.acs.org/doi/10.1021/acs.organomet.8b00456) (Eq. 1 - 4) from the gas-phase Gibbs free energy G(gas) and the solvation free energy G(solv).

The plots are surely not of publication quality but allow a quick check of reaction energetics.


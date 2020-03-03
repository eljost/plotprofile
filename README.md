# plotprofile
Quickly generate plots of reaction solution-phase free energies (G(sol)) and
integrate the corresponding rate-equations.

## Installation
```
git clone https://github.com/eljost/plotprofile.git
cd plotprofile
python setup.py install
```

## Usage
```
plotprofile [yaml]
```

Calculated data is dumped to a CSV file.


## Overview

All relevant energy values are directly parsed from logfiles. In this regard no manual
work is required. The necessary program input is read from a YAML file. An example can
be found in the .plotprofile/tests subfolder.

There are three mandatory sections in the YAML file: `program`, `molecules`, and `reactions`.
Optional sections include `paths` and `kinetics`.

### program
Right now only `orca` is supported.

### molecules
Three different logfiles can be provided: freq (thermochemistry data), solv (solvated calculation) and scf (alternative single point energy, e.g. at a higher level/bigger basis set etc.)
For improved printing a label can be given, using matplotib syntax ($cis$-alkene etc.).

```
molecules:
 hcn:
  freq: 05_hcn_iso_educt_freq.out
  solv: 06_hcn_iso_educt_solv.out
  scf: 10_hcn_iso_educt_dlpnoccsdt.out
  label: HCN
```

### reactions
Every reaction is comprised of three parts (required): educts, ts, and products. The key `add` can
be used to add a constant energy by a molecule to all reagents of a reaction. This can be useful
if e.g. a water molecule is not present in the first reaction step but introduced later on in
the remaining reactions. Similar to `molecules` a label can be given.

```
reactions:
 alkin_addition:
  educts: [alkin, pbu3]
  ts: ts1
  products: addukt
  add: water
  label: Alkye addition
 water_addition:
  educts: [addukt, water]
  ts: ts2
  products: addukt_h_oh
  label: Water addition
 ...
```

### paths
Paths are made up from subsequently occuring reactions. 


```
paths:
 path_name:
  - reaction0
  - reaction1
  - ...
```

### kinetics
If a `kinetics` section is given the corresponding rate equations will be set up and integrated.
Two subsections are required: `c0s` specifying the inital concentrations and `t_span` specifying
the integration interval in seconds. 

```
kinetics:
 # in mol/l
 c0s:
  alkin: 1
  pbu3: 1
  water: 10
 # in seconds
 t_span: [0, 86400]
```

## Remarks
G(sol) is calculated as described in [10.1021/acs.organomet.8b00456](https://pubs.acs.org/doi/10.1021/acs.organomet.8b00456) (Eq. 1 - 4) from the gas-phase Gibbs free energy G(gas) and the solvation free energy G(solv).

The plots are surely not of publication quality but allow a quick check of reaction energetics.


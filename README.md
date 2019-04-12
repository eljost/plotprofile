# plotprofile
Quickly generate plots of solution-phase free energies (G(sol)) for reaction paths.

The plots are surely not of publication quality but allow a quick check of reaction energetics.

G(sol) is calculated as described in [10.1021/acs.organomet.8b00456](https://pubs.acs.org/doi/10.1021/acs.organomet.8b00456) (Eq. 1 - 4) from the gas-phase Gibbs free energy G(gas) and the solvation free energy G(solv).

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

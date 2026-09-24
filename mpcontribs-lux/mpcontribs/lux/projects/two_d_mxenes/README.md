# two_d_mxenes

Relaxed structures and computed properties of 2D MXenes, M<sub>n+1</sub>X<sub>n</sub>T<sub>x</sub>, systematically varying

- the transition metal **M**: Ti, Mo, Hf, Re
- the non-metal **X**: C, N
- the thickness **n**: 1, 2, 3 (M<sub>2</sub>X, M<sub>3</sub>X<sub>2</sub>, M<sub>4</sub>X<sub>3</sub>)
- the surface termination **T**: none (pristine), F, O
- the stacking of layers: octahedral (O) or trigonal-prismatic (P) coordination of each layer

**Reference:** N. Oyeniran, O. Chowdhury, C. Hu, T. Dumitrica, P. Ganesh, J. Jakowski, Z. Chen, R. R. Unocic, M. Naguib, V. Meunier, Y. Gogotsi, P. R. C. Kent, B. G. Sumpter, J. Huang, "A Panoramic View of MXenes via a New Design Strategy", *Adv. Funct. Mater.* (2025), [doi:10.1002/adfm.202508047](https://doi.org/10.1002/adfm.202508047); preprint [arXiv:2501.15390](https://doi.org/10.48550/arXiv.2501.15390).

**Raw data:** [Materials Data Facility](https://www.materialsdatafacility.org/detail/a65168f7-8f13-4552-b660-c1565f6d093e-1.0)

**Contributors:** Aidan Gesch (University of Alabama, Hu group)

## Layout

```
two_d_mxenes/
├── schemas/
│   ├── labels.py        M, X, T, n, stacking label, termination site (+ consistency rules)
│   ├── structure.py     relaxed structure + descriptors derived from it
│   ├── properties.py    energetics and elastic properties from the authors' spreadsheet
│   ├── calculation.py   project-level DFT settings (identical for all entries)
│   └── mxene.py         MXeneEntry: one record per structure, cross-validated
└── pipelines/
    └── build_contributions.py   CONTCAR tree + spreadsheet -> entries -> contributions
```

## Stacking labels

Each atomic layer that lies between two others is **O** (octahedral) if the layers below and above it are staggered, or **P** (trigonal prismatic) if they are eclipsed. `StructureDescriptors.coordinationSequence` measures this directly from each structure. The authors' labels map onto the sequence of layers between the two outer metal layers as follows:

| label | n | core sequence |
|---|---|---|
| `t` | 1 | O |
| `h` | 1 | P |
| `t` | 2, 3 | O-O-O… |
| `h2` | 2, 3 | P-P-P… |
| `h1a` | 2, 3 | O-P-O… |
| `h1b` | 2, 3 | P-O-P… |

For terminated MXenes the suffix `-1`/`-2` selects the termination site, which sets the coordination of the outer metal layers (first and last entries of the sequence). `MXeneEntry` rejects any record whose structure does not match its label.

## Record contents

| group | fields | source |
|---|---|---|
| `labels` | metal, nonmetal, termination, n, stacking, terminationSite | CONTCAR composition + folder name |
| `structure` | lattice, species, fracCoords | CONTCAR |
| `descriptors` | a, b, gamma, cellArea, areaPerFormulaUnit, arealMassDensity, thickness, vacuum, layerSequence, coordinationSequence, M–X and M–T bond lengths, space group | computed from structure |
| `properties.energetics` | totalEnergyPerAtom, formationEnergyPerAtom, relativeStackingEnergy | spreadsheet (relative energy computed) |
| `properties.elastic` | C11, C12, C66, Young's modulus, Poisson's ratio, shear modulus, mechanical stability | spreadsheet (moduli computed) |

## Spreadsheet format (proposed)

One row per structure, one sheet, a single header row with **exactly** these column names:

| column | unit | required | notes |
|---|---|---|---|
| `mxeneId` | – | yes | formula + folder label, e.g. `Ti3C2O2-h1a-2`, `Mo2N-t` |
| `totalEnergyPerAtom` | eV/atom | recommended | used to compute `relativeStackingEnergy` |
| `formationEnergyPerAtom` | eV/atom | optional | state the reference energies in the project description |
| `c11` | N/m | optional | 2D elastic constant |
| `c12` | N/m | optional | 2D elastic constant |
| `c66` | N/m | optional | leave blank to use (C11 − C12)/2 |

- Leave cells **blank** for missing values (no `N/A`, `-`, or `0` placeholders).
- Numbers only in data cells; units live in this table, not the cells.
- Lattice parameters are **not** needed: they are computed from each CONTCAR.
- Young's modulus, Poisson's ratio, shear modulus and stability are **derived** from C11/C12/C66, so they need not be supplied.
- The pipeline refuses unknown columns, duplicate IDs, and rows without a matching structure, so problems surface before upload.

Additional columns (e.g. band gap, work function) can be added by extending `properties.py` and `SPREADSHEET_COLUMNS` together.

## Simulation settings

DFT settings are the same for every entry and are recorded once, at the project level, using `CalculationSettings`.

"""Build validated MXene entries and MPContribs contributions from raw data.

Usage::

    from mpcontribs.lux.projects.two_d_mxenes.pipelines.build_contributions import (
        build_entries, load_properties, to_contribution,
    )

    props = load_properties("mxene_properties.xlsx")
    entries = build_entries("MXENE_DATA", properties=props)
    contributions = [to_contribution(e) for e in entries]

Labels are assigned as follows, which keeps the pipeline independent of how
the dataset's folders are nested (e.g. the extra `ReC/`, `ReN/` levels):

- M, X, T and n are inferred from the composition of each CONTCAR;
- the stacking label and termination site come from the name of the folder
  that directly contains the CONTCAR (e.g. `h1a-2`, `t`).
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterator, Mapping
from pathlib import Path

import pandas as pd
from mpcontribs.lux.projects.two_d_mxenes.schemas import (
    ElasticProperties,
    Energetics,
    MXeneEntry,
    MXeneLabel,
    MXeneProperties,
)
from pymatgen.core import Structure

STRUCTURE_FILENAME = "CONTCAR"

SPREADSHEET_COLUMNS: dict[str, str] = {
    "mxeneId": "",
    "totalEnergyPerAtom": "eV/atom",
    "formationEnergyPerAtom": "eV/atom",
    "c11": "N/m",
    "c12": "N/m",
    "c66": "N/m",
}
"""Expected spreadsheet columns and their units (see the project README)."""


def iter_structure_files(root: str | Path) -> Iterator[Path]:
    """Yield every CONTCAR below `root`, in sorted order."""
    yield from sorted(Path(root).rglob(STRUCTURE_FILENAME))


def entry_from_file(
    path: str | Path, properties: MXeneProperties | None = None
) -> MXeneEntry:
    """Build one validated entry from a CONTCAR and its folder label."""
    path = Path(path)
    stacking, site = MXeneLabel.parse_folder_label(path.parent.name)
    structure = Structure.from_file(path)
    return MXeneEntry.from_structure(
        structure, stacking=stacking, terminationSite=site, properties=properties
    )


def load_properties(path: str | Path) -> dict[str, MXeneProperties]:
    """Read the properties spreadsheet (xlsx or csv) keyed by `mxeneId`."""
    path = Path(path)
    df = pd.read_csv(path) if path.suffix == ".csv" else pd.read_excel(path)
    missing = {"mxeneId"} - set(df.columns)
    if missing:
        raise ValueError(f"Spreadsheet is missing required columns {missing}")
    unknown = set(df.columns) - set(SPREADSHEET_COLUMNS)
    if unknown:
        raise ValueError(f"Spreadsheet has unrecognized columns {sorted(unknown)}")
    if df["mxeneId"].duplicated().any():
        dupes = df.loc[df["mxeneId"].duplicated(), "mxeneId"].tolist()
        raise ValueError(f"Duplicate mxeneId values: {dupes}")

    df = df.astype(object).where(pd.notna(df), None)
    out: dict[str, MXeneProperties] = {}
    for row in df.to_dict(orient="records"):
        c11, c12, c66 = row.get("c11"), row.get("c12"), row.get("c66")
        elastic = (
            ElasticProperties.from_elastic_constants(c11, c12, c66)
            if c11 is not None and c12 is not None
            else None
        )
        energetics = Energetics(
            totalEnergyPerAtom=row.get("totalEnergyPerAtom"),
            formationEnergyPerAtom=row.get("formationEnergyPerAtom"),
        )
        out[row["mxeneId"]] = MXeneProperties(energetics=energetics, elastic=elastic)
    return out


def build_entries(
    root: str | Path, properties: Mapping[str, MXeneProperties] | None = None
) -> list[MXeneEntry]:
    """Build and cross-validate entries for every CONTCAR under `root`.

    Raises if two structures map to the same `mxeneId`, or if the spreadsheet
    lists an `mxeneId` with no matching structure.
    """
    properties = dict(properties or {})
    entries: dict[str, MXeneEntry] = {}
    for path in iter_structure_files(root):
        entry = entry_from_file(path)
        if entry.mxeneId in entries:
            raise ValueError(f"Duplicate mxeneId {entry.mxeneId!r} at {path}")
        entry.properties = properties.pop(entry.mxeneId, None)
        entries[entry.mxeneId] = entry
    if properties:
        raise ValueError(
            f"No structure found for spreadsheet rows {sorted(properties)}"
        )

    result = list(entries.values())
    add_relative_stacking_energies(result)
    return [MXeneEntry.model_validate(e.model_dump()) for e in result]


def add_relative_stacking_energies(entries: list[MXeneEntry]) -> None:
    """Fill `relativeStackingEnergy` (meV/atom) within each composition."""
    groups: dict[str, list[MXeneEntry]] = defaultdict(list)
    for entry in entries:
        energetics = entry.properties and entry.properties.energetics
        if energetics and energetics.totalEnergyPerAtom is not None:
            groups[entry.labels.formula].append(entry)
    for group in groups.values():
        e_min = min(e.properties.energetics.totalEnergyPerAtom for e in group)
        for entry in group:
            energetics = entry.properties.energetics
            energetics.relativeStackingEnergy = 1000.0 * (
                energetics.totalEnergyPerAtom - e_min
            )


def to_contribution(entry: MXeneEntry, project: str = "two_d_mxenes") -> dict:
    """Convert an entry into an MPContribs contribution dictionary.

    Only searchable quantities go into `data` (MPContribs allows at most 50
    flattened keys); values carry units as strings, as MPContribs expects.
    The full record, including arrays, lives in the project's Parquet file.
    """
    lab, desc = entry.labels, entry.descriptors
    data: dict = {
        "mxeneId": entry.mxeneId,
        "M": lab.metal,
        "X": lab.nonmetal,
        "T": lab.termination or "none",
        "n": lab.n,
        "stacking": lab.stacking,
        "terminationSite": lab.terminationSite or "none",
        "coordination": "-".join(desc.coordinationSequence),
        "structure": {
            "a": f"{desc.a:.4f} Å",
            "thickness": f"{desc.thickness:.4f} Å",
            "MX": f"{desc.metalNonmetalBondLength:.4f} Å",
            "spaceGroup": desc.spaceGroupSymbol,
        },
    }
    if desc.metalTerminationBondLength is not None:
        data["structure"]["MT"] = f"{desc.metalTerminationBondLength:.4f} Å"

    props = entry.properties
    if props and props.energetics:
        en = props.energetics
        energy = {
            "formation": _with_unit(en.formationEnergyPerAtom, "eV/atom"),
            "relativeStacking": _with_unit(en.relativeStackingEnergy, "meV/atom"),
        }
        if energy := {k: v for k, v in energy.items() if v is not None}:
            data["energy"] = energy
    if props and props.elastic:
        el = props.elastic
        elastic = {
            "C11": _with_unit(el.c11, "N/m"),
            "C12": _with_unit(el.c12, "N/m"),
            "C66": _with_unit(el.c66, "N/m"),
            "Y": _with_unit(el.youngsModulus, "N/m"),
            "nu": _with_unit(el.poissonRatio, ""),
            "stable": (
                None if el.mechanicallyStable is None else str(el.mechanicallyStable)
            ),
        }
        if elastic := {k: v for k, v in elastic.items() if v is not None}:
            data["elastic"] = elastic

    return {
        "project": project,
        "identifier": entry.descriptors.reducedFormula,
        "data": data,
        "structures": [entry.structure.pymatgen_structure],
    }


def _with_unit(value: float | None, unit: str) -> str | None:
    if value is None:
        return None
    return f"{value:.6g} {unit}".strip()

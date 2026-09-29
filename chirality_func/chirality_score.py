#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Final robust version with correct stereochemistry logic
and corrected q1 computation:
- q1 now sums ONLY positions of stereogenic atoms (R, S, U).
-------------------------------------------------------
- Keeps all original Excel columns
- Preserves row order
- Handles "m.d." cases
- Handles invalid SMILES (no crash)
- Correct CIP logic:
      R  → R
      S  → S
      ?  → U  (true unknown stereocenter)
   no CIP → N (not stereogenic)
- Adds columns:
      inchi
      absolute_configuration
      stereo_positions
      q0
      q1
      q2
      f_chiral
      observations
"""

import argparse
import json
import math
from typing import List, Dict, Tuple, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import inchi as rd_inchi
except ImportError:
    raise SystemExit("RDKit is required in your chemo_env environment.")

try:
    import pandas as pd
except Exception:
    pd = None


SCORE_MAP = {"R": 1.0, "S": 0.5, "U": 0.75, "N": 0.0}


# ---------------------------------------------------------
# Safe SMILES → Mol
# ---------------------------------------------------------

def smiles_to_mol(smiles: str) -> Optional[Chem.Mol]:
    """Convert SMILES to RDKit Mol. Return None if invalid."""
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            return None
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True,
                                   flagPossibleStereoCenters=True)
        return mol
    except Exception:
        return None


def mol_to_inchi(mol: Chem.Mol) -> str:
    try:
        return rd_inchi.MolToInchi(mol)
    except Exception:
        return ""


# ---------------------------------------------------------
# Atom labeling R/S/U/N (correct logic)
# ---------------------------------------------------------

def label_atoms(mol: Chem.Mol) -> Tuple[List[str], List[int], List[int]]:
    """
    Returns:
        labels: per-atom labels R/S/U/N
        center_indices: indices of stereogenic centers (R, S, U)
        stereo_positions_1based: positions in IUPAC 1-based numbering
    """
    n = mol.GetNumAtoms()
    labels = ["N"] * n   # default: NOT stereogenic

    chiral_list = Chem.FindMolChiralCenters(
        mol, includeUnassigned=True, useLegacyImplementation=False
    )
    chiral_set = set(idx for idx, _ in chiral_list)

    for idx, _ in chiral_list:
        atom = mol.GetAtomWithIdx(idx)
        cip = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None

        if cip == "R":
            labels[idx] = "R"
        elif cip == "S":
            labels[idx] = "S"
        elif cip == "?":
            labels[idx] = "U"  # true unknown stereocenter
        else:
            labels[idx] = "N"  # not stereogenic

    # Centers are only R, S, or U
    center_indices = sorted(
        idx for idx in range(n) if labels[idx] in ("R", "S", "U")
    )
    positions_1based = [i + 1 for i in center_indices]

    return labels, center_indices, positions_1based


# ---------------------------------------------------------
# q0, q1, q2 and f_chiral
# ---------------------------------------------------------

def compute_q_values(labels: List[str]) -> Tuple[float, float, float]:
    q0 = q1 = q2 = 0.0

    for i, lab in enumerate(labels):

        # Position is always 1-based:
        x = float(i + 1)

        # Score:
        s = SCORE_MAP.get(lab, 0.0)

        # q0 sums ALL s-values (R/S/U/N), but s(N)=0:
        q0 += s

        # q1 must sum ONLY positions of stereogenic atoms:
        if lab in ("R", "S", "U"):
            q1 += x

        # q2 uses x*s; for N, s=0 so it contributes nothing:
        q2 += x * s

    return q0, q1, q2


def f_chiral_from_q(q0: float, q1: float, q2: float) -> float:
    if q0 == 0.0 or q1 == 0.0:
        return 0.0
    return 1.0 / (1.0 + math.exp(-(q2 / q1)))


# ---------------------------------------------------------
# Process one record safely
# ---------------------------------------------------------

def process_record(chembl_id: str, smiles: str) -> Dict[str, object]:

    # m.d. handling
    if (chembl_id or "").strip().lower() == "m.d." or (smiles or "").strip().lower() == "m.d.":
        return {
            "inchi": "",
            "absolute_configuration": "",
            "stereo_positions": "",
            "q0": 0.0,
            "q1": 0.0,
            "q2": 0.0,
            "f_chiral": 0.25,
            "observations": "m.d."
        }

    mol = smiles_to_mol(smiles)
    if mol is None:
        return {
            "inchi": "",
            "absolute_configuration": "",
            "stereo_positions": "",
            "q0": 0.0,
            "q1": 0.0,
            "q2": 0.0,
            "f_chiral": 0.25,
            "observations": "Invalid SMILES"
        }

    inchi = mol_to_inchi(mol)
    labels, center_indices, pos_1based = label_atoms(mol)
    config_centers = [labels[i] for i in center_indices]

    q0, q1, q2 = compute_q_values(labels)
    f_chiral = 0.0 if q0 == 0.0 else f_chiral_from_q(q0, q1, q2)

    return {
        "inchi": inchi,
        "absolute_configuration": ",".join(config_centers) if config_centers else "",
        "stereo_positions": ",".join(str(x) for x in pos_1based) if pos_1based else "",
        "q0": q0,
        "q1": q1,
        "q2": q2,
        "f_chiral": f_chiral,
        "observations": "OK"
    }


# ---------------------------------------------------------
# Process Excel
# ---------------------------------------------------------

def process_excel(in_path: str, sheet: str, id_col: str, smiles_col: str, out_path: str):
    if pd is None:
        raise SystemExit("pandas is required. Install with: pip install pandas openpyxl")

    df = pd.read_excel(in_path, sheet_name=sheet)

    if id_col not in df.columns or smiles_col not in df.columns:
        raise SystemExit(f"Columns '{id_col}' and '{smiles_col}' must be present in the input file.")

    results = []
    for _, row in df.iterrows():
        chembl_id = str(row[id_col]) if pd.notna(row[id_col]) else ""
        smiles = str(row[smiles_col]) if pd.notna(row[smiles_col]) else ""
        res = process_record(chembl_id, smiles)
        results.append(res)

    res_df = pd.DataFrame(results)
    combined = pd.concat([df.reset_index(drop=True), res_df.reset_index(drop=True)], axis=1)
    combined.to_excel(out_path, index=False)
    print(f"Done. Results written to {out_path} ({len(combined)} rows).")


# ---------------------------------------------------------
# CLI
# ---------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compute chirality function from SMILES.")
    parser.add_argument("--excel", required=True)
    parser.add_argument("--sheet", required=True)
    parser.add_argument("--id-col", required=True)
    parser.add_argument("--smiles-col", required=True)
    parser.add_argument("--out", default="salida.xlsx")
    args = parser.parse_args()

    process_excel(args.excel, args.sheet, args.id_col, args.smiles_col, args.out)


if __name__ == "__main__":
    main()

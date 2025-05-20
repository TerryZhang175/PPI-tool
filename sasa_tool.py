#!/usr/bin/env python3
"""sasa_tool.py

Standalone module providing functions to calculate Solvent Accessible Surface
Area (SASA) for residues in PDB files using the FreeSASA library.  The key
functions – ``process_single_pdb_file`` and ``process_pdb_directory`` – are
imported by ``ppi_suite.py`` but the script can also be executed directly from
the command line.
"""

from __future__ import annotations

import argparse
import os
import re
from typing import Dict

import freesasa

__all__ = [
    "process_single_pdb_file",
    "process_pdb_directory",
    "calculate_cys_sasa",
    "calculate_residue_sasa",
]


class CustomClassifier(freesasa.Classifier):
    """Custom atomic‐radius classifier.

    The default classifier in *FreeSASA* relies on element and residue‐specific
    radii tables.  For our purpose a simplified model is sufficient and keeps
    the dependency footprint small.  The radii values come from standard VDW
    radii frequently used in structural biology.
    """

    purePython = True  # signal that the implementation is in pure Python

    def radius(self, residueName: str, atomName: str) -> float:  # noqa: N802 – names imposed by freesasa
        if re.match(r"^\s*N", atomName, re.IGNORECASE):
            return 1.6
        if re.match(r"^\s*C", atomName, re.IGNORECASE):
            return 1.7
        if re.match(r"^\s*O", atomName, re.IGNORECASE):
            return 1.4
        if re.match(r"^\s*S", atomName, re.IGNORECASE):
            return 1.8
        # Return 0 for unknown atoms so they are ignored by FreeSASA.
        return 0.0


# ---------------------------------------------------------------------------
# Core calculation helpers
# ---------------------------------------------------------------------------

def _sasa_per_atom(pdb_path: str, classifier: freesasa.Classifier | None = None) -> tuple[freesasa.Structure, freesasa.Result]:
    """Return *FreeSASA* Structure and Result for *pdb_path*."""
    classifier = classifier or CustomClassifier()
    structure = freesasa.Structure(pdb_path, classifier)
    result = freesasa.calc(structure)
    return structure, result


def calculate_residue_sasa(pdb_path: str, residue_type: str | None = None) -> Dict[str, float]:
    """Return total SASA per residue.

    Parameters
    ----------
    pdb_path
        Path to a PDB file.
    residue_type
        Three‐letter residue code (e.g. "CYS").  If *None* SASA for **all**
        residues is returned.
    """
    structure, result = _sasa_per_atom(pdb_path)

    residue_sasa: Dict[str, float] = {}
    for i in range(structure.nAtoms()):
        res_name = structure.residueName(i)
        if residue_type is None or res_name == residue_type:
            res_num_str = structure.residueNumber(i)
            # Ensure res_name and res_num_str are stripped and res_num is formatted correctly
            cleaned_res_name = str(res_name).strip()
            cleaned_res_num = str(res_num_str).strip()
            key = f"{cleaned_res_name}_{cleaned_res_num}"
            residue_sasa[key] = residue_sasa.get(key, 0.0) + result.atomArea(i)
    return residue_sasa


def calculate_cys_sasa(pdb_path: str) -> Dict[int, float]:
    """Convenience wrapper returning SASA per *CYS* residue."""
    structure, result = _sasa_per_atom(pdb_path)

    cys_sasa: Dict[int, float] = {}
    for i in range(structure.nAtoms()):
        if structure.residueName(i) == "CYS":
            res_num_str = structure.residueNumber(i) # Get as string or whatever it returns
            try:
                # Convert to int after stripping whitespace
                cleaned_res_num = int(str(res_num_str).strip()) 
                cys_sasa[cleaned_res_num] = cys_sasa.get(cleaned_res_num, 0.0) + result.atomArea(i)
            except ValueError:
                print(f"[WARNING sasa_tool] Could not convert residue number '{res_num_str}' to int for CYS in {pdb_path}. Skipping this atom.")
    return cys_sasa


# ---------------------------------------------------------------------------
# File and directory helpers
# ---------------------------------------------------------------------------

def _write_sasa_to_file(pdb_name: str, residue_sasa: Dict[str, float] | Dict[int, float], output_file: str) -> None:
    """Write *residue_sasa* dictionary to *output_file*."""
    with open(output_file, "w", encoding="utf-8") as fh:
        fh.write(f"{pdb_name}\n")
        for res_id, area in sorted(residue_sasa.items()):
            if isinstance(res_id, int):
                # CYS only – integer residue number
                fh.write(f"CYS {res_id}  SASA {area:.2f}\n")
            else:
                res_type, res_num = res_id.split("_", 1)
                fh.write(f"{res_type} {res_num}  SASA {area:.2f}\n")
        fh.write("\n")


def process_single_pdb_file(pdb_path: str, output_dir: str | None = None, residue_type: str = "CYS", verbose: bool = True) -> bool:
    """Calculate SASA for *pdb_path* and write results to *output_dir*."""
    if verbose:
        print(f"[SASA] Processing single PDB file: {pdb_path}")

    output_dir = output_dir or os.path.join("./sasa_results", "single_pdb_sasa")
    os.makedirs(output_dir, exist_ok=True)

    pdb_name = os.path.basename(pdb_path)
    output_file = os.path.join(output_dir, pdb_name.replace(".pdb", "_sasa.txt"))

    try:
        if residue_type.upper() == "CYS":
            residue_sasa = calculate_cys_sasa(pdb_path)
        else:
            residue_sasa = calculate_residue_sasa(pdb_path, residue_type.upper())
        _write_sasa_to_file(pdb_name, residue_sasa, output_file)
        if verbose:
            print(f"[SASA] Results written to {output_file}")
        return True
    except Exception as exc:  # pylint: disable=broad-except – FreeSASA may raise many exceptions
        print(f"[SASA] Error processing {pdb_name}: {exc}")
        return False


def process_pdb_directory(folder_path: str, output_dir: str | None = None, residue_type: str = "CYS") -> None:
    """Traverse *folder_path* and calculate SASA for each PDB file."""
    print(f"[SASA] Starting batch calculation for directory: {folder_path}")

    output_dir = output_dir or os.path.join("./sasa_results", f"{os.path.basename(os.path.normpath(folder_path))}_sasa")
    os.makedirs(output_dir, exist_ok=True)

    pdb_files = sorted(f for f in os.listdir(folder_path) if f.lower().endswith(".pdb"))

    processed, errors = 0, 0
    for pdb_name in pdb_files:
        pdb_path = os.path.join(folder_path, pdb_name)
        if process_single_pdb_file(pdb_path, output_dir, residue_type, verbose=False):
            processed += 1
        else:
            errors += 1

    print(f"[SASA] Completed: {processed} files processed, {errors} errors.  Results in {output_dir}")


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def _build_cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Calculate Solvent Accessible Surface Area (SASA) for protein structures")
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("-f", "--file", help="Path to a single PDB file")
    grp.add_argument("-d", "--directory", help="Path to a directory containing PDB files")

    parser.add_argument("-o", "--output", help="Output directory path (default: ./sasa_results)")
    parser.add_argument("-r", "--residue", default="CYS", help="Three-letter residue code to analyse (default: CYS)")
    return parser


def main() -> None:  # pragma: no cover – CLI entry-point
    args = _build_cli().parse_args()
    if args.file:
        process_single_pdb_file(args.file, args.output, args.residue.upper())
    else:
        process_pdb_directory(args.directory, args.output, args.residue.upper())


if __name__ == "__main__":
    main() 
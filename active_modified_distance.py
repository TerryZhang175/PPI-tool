import os
import pandas as pd
from Bio import PDB
import math
import sys
import argparse # Added import for argparse for standalone execution

# Global variables to be set by CLI arguments via ppi_suite or __main__
INPUT_CSV_PATH = ""
PDB_FOLDER_PATH = "" # Changed from PDB_FILE_PATH
OUTPUT_FILE = "active_modified_distances.txt"  # Default output filename changed back to .txt

def calc_distance(atom1, atom2):
    """
    Calculate the Euclidean distance between two atoms.
    """
    # diff_vector = atom1.coord - atom2.coord # This causes TypeError with Bio.PDB.Atom.Atom
    # return math.sqrt(sum(x*x for x in diff_vector))
    # Corrected distance calculation for Bio.PDB.Atom objects
    x = atom1.coord[0] - atom2.coord[0]
    y = atom1.coord[1] - atom2.coord[1]
    z = atom1.coord[2] - atom2.coord[2]
    return math.sqrt(x*x + y*y + z*z)

def get_atom_sg_or_ca(structure, residue_number, pdb_file_identifier_for_error):
    """
    Finds the specified atom (SG for CYS, CA otherwise) in a given residue number
    within the first chain of the provided PDB structure (assumed monomer).

    Args:
        structure (Bio.PDB.Structure.Structure): Parsed PDB structure.
        residue_number (int): The residue number to look for.
        pdb_file_identifier_for_error (str): Identifier of the PDB file for error messages.

    Returns:
        tuple: (Bio.PDB.Atom.Atom, str) - The found atom and its type ('SG' or 'CA').

    Raises:
        ValueError: If the model has no chains or the residue/atom is not found.
    """
    model = structure[0]  # Assuming the first model
    if not model.child_list:
        raise ValueError(f"PDB model in {pdb_file_identifier_for_error} contains no chains.")

    chain = model.child_list[0]  # Assuming the first chain for a monomer

    try:
        # Correctly form the residue identifier tuple for BioPython
        # It should be (' ', residue_number, ' ') for hetfield, resseq, icode
        # Assuming no HETATM and no insertion codes for simplicity here.
        residue_id_tuple = (' ', residue_number, ' ') 
        residue = chain[residue_id_tuple]
    except KeyError:
        raise ValueError(f"Residue {residue_number} (ID: {residue_id_tuple}) not found in chain {chain.id} of {pdb_file_identifier_for_error}.")

    if residue.resname == 'CYS' and 'SG' in residue:
        return residue['SG'], 'SG'
    elif 'CA' in residue:
        return residue['CA'], 'CA'
    else:
        raise ValueError(f"Neither SG (for CYS) nor CA atom found in residue {residue_number} (Resname: {residue.resname}) of chain {chain.id} in {pdb_file_identifier_for_error}.")

def process_residue_pairs():
    """
    Main processing logic: reads CSV, iterates through PDBs based on CSV,
    parses PDBs, calculates distances, writes output.
    Uses global variables INPUT_CSV_PATH, PDB_FOLDER_PATH, OUTPUT_FILE.
    """
    if not INPUT_CSV_PATH or not os.path.exists(INPUT_CSV_PATH):
        print(f"Error: Input CSV file not set or found: {INPUT_CSV_PATH}", file=sys.stderr)
        return False
    if not PDB_FOLDER_PATH or not os.path.isdir(PDB_FOLDER_PATH): # Check if it's a directory
        print(f"Error: PDB folder not set or not a valid directory: {PDB_FOLDER_PATH}", file=sys.stderr)
        return False

    pdb_parser = PDB.PDBParser(QUIET=True)
    
    try:
        df_pairs = pd.read_csv(INPUT_CSV_PATH)
        if len(df_pairs.columns) < 3:
            print(f"Error: CSV file '{INPUT_CSV_PATH}' must have at least three columns (PDB_Name, active site, modified site).", file=sys.stderr)
            return False
        pdb_name_col, active_site_col, modified_site_col = df_pairs.columns[0], df_pairs.columns[1], df_pairs.columns[2]
        print(f"Reading PDB names from CSV column: '{pdb_name_col}'.")
        print(f"Reading residue pairs from CSV columns: '{active_site_col}' (Active Site) and '{modified_site_col}' (Modified Site).")
    except Exception as e:
        print(f"Error reading or parsing CSV file '{INPUT_CSV_PATH}': {e}", file=sys.stderr)
        return False

    output_dir = os.path.dirname(OUTPUT_FILE)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            print(f"Error creating output directory '{output_dir}': {e}", file=sys.stderr)
            return False

    processed_count = 0
    error_count = 0
    
    # Keep track of loaded structures to avoid reparsing the same PDB file multiple times
    loaded_structures = {}

    with open(OUTPUT_FILE, "w", encoding="utf-8") as out_f:
        header = ["PDB_Identifier", "ActiveSite_ResNum", "ModifiedSite_ResNum", "Distance(Å)", "ActiveSite_AtomType", "ModifiedSite_AtomType"]
        out_f.write("\t".join(header) + "\n")

        for index, row in df_pairs.iterrows():
            try:
                pdb_name_from_csv = str(row[pdb_name_col]).strip()
                active_res_num_str = str(row[active_site_col])
                modified_res_num_str = str(row[modified_site_col])
                
                active_res_num = int(float(active_res_num_str))
                modified_res_num = int(float(modified_res_num_str))

                # Construct potential PDB file paths
                pdb_file_path_pdb = os.path.join(PDB_FOLDER_PATH, f"{pdb_name_from_csv}.pdb")
                pdb_file_path_cif = os.path.join(PDB_FOLDER_PATH, f"{pdb_name_from_csv}.cif")
                
                actual_pdb_path = None
                pdb_identifier = pdb_name_from_csv

                if os.path.exists(pdb_file_path_pdb):
                    actual_pdb_path = pdb_file_path_pdb
                elif os.path.exists(pdb_file_path_cif):
                    actual_pdb_path = pdb_file_path_cif
                else:
                    raise FileNotFoundError(f"PDB/CIF file for '{pdb_name_from_csv}' not found in {PDB_FOLDER_PATH} (checked for .pdb and .cif extensions).")

                # Load structure if not already loaded
                if actual_pdb_path not in loaded_structures:
                    try:
                        structure = pdb_parser.get_structure(pdb_identifier, actual_pdb_path)
                        loaded_structures[actual_pdb_path] = structure
                    except Exception as e:
                        raise RuntimeError(f"Error parsing PDB/CIF file '{actual_pdb_path}': {e}")
                else:
                    structure = loaded_structures[actual_pdb_path]
                
                atom1, atom1_type = get_atom_sg_or_ca(structure, active_res_num, actual_pdb_path)
                atom2, atom2_type = get_atom_sg_or_ca(structure, modified_res_num, actual_pdb_path)

                distance = calc_distance(atom1, atom2)
                
                # Write tab-separated line
                line_data = [str(pdb_identifier), str(active_res_num), str(modified_res_num), f"{distance:.3f}", str(atom1_type), str(atom2_type)]
                out_f.write("\t".join(line_data) + "\n")
                processed_count += 1

            except FileNotFoundError as fnf_ex:
                print(f"Skipping CSV row {index + 1} (PDB: {pdb_name_from_csv}, Active: {row.get(active_site_col, 'N/A')}, Mod: {row.get(modified_site_col, 'N/A')}): {fnf_ex}", file=sys.stderr)
                error_count += 1
            except RuntimeError as r_ex: # For PDB parsing errors within the loop
                print(f"Skipping CSV row {index + 1} (PDB: {pdb_name_from_csv}, Active: {row.get(active_site_col, 'N/A')}, Mod: {row.get(modified_site_col, 'N/A')}): {r_ex}", file=sys.stderr)
                error_count +=1
            except ValueError as ve:
                print(f"Skipping CSV row {index + 1} (PDB: {pdb_name_from_csv}, Active: {row.get(active_site_col, 'N/A')}, Mod: {row.get(modified_site_col, 'N/A')}): {ve}", file=sys.stderr)
                error_count += 1
            except Exception as ex:
                print(f"Unexpected error processing CSV row {index + 1} (PDB: {pdb_name_from_csv}, Active: {row.get(active_site_col, 'N/A')}, Mod: {row.get(modified_site_col, 'N/A')}): {ex}", file=sys.stderr)
                error_count += 1

    print(f"\nProcessing complete. Results written to '{OUTPUT_FILE}'.")
    print(f"Successfully processed {processed_count} residue pairs from the CSV.")
    if error_count > 0:
        print(f"Encountered {error_count} errors (see messages above).")
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate intramolecular distances between active sites and modified sites for multiple PDB/CIF files in a folder, based on a CSV input.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--input-csv', required=True,
                        help='Path to the input CSV file.\nCSV format should have three columns (headers are optional but recommended):\n  Column 1: PDB file name (without .pdb or .cif extension)\n  Column 2: Active site residue number\n  Column 3: Modified site residue number')
    parser.add_argument('--pdb-folder', required=True, # Changed from --pdb-file
                        help='Path to the folder containing PDB or mmCIF files.')
    parser.add_argument('--output-file', default=OUTPUT_FILE,
                        help=f'Output file path for the results (default: {OUTPUT_FILE})')

    args = parser.parse_args()

    # Set global variables from CLI args
    INPUT_CSV_PATH = args.input_csv
    PDB_FOLDER_PATH = args.pdb_folder # Changed from PDB_FILE_PATH
    OUTPUT_FILE = args.output_file

    if not process_residue_pairs():
        sys.exit(1) 
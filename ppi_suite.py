#!/usr/bin/env python3
import argparse
import importlib

# Import SASA functions from sasa_tool.py
from sasa_tool import process_single_pdb_file, process_pdb_directory

def run_sasa(args):
    """Run SASA calculation using the sasa_tool module."""
    if args.file:
        process_single_pdb_file(args.file, args.output, args.residue)
    elif args.directory:
        process_pdb_directory(args.directory, args.output, args.residue)

def run_ppi_glutathionylation_check(args):
    """Run PPI interface glutathionylation check."""
    # This imports the ppi_interface_glutathionylation.py script
    import ppi_interface_glutathionylation 
    ppi_interface_glutathionylation.process_structures(args.input_csv, args.pdb_folder, args.output_csv, threshold=args.threshold)

def run_active_modified_distance(args):
    """Run active site to modified site distance calculation for monomers using a CSV input."""
    # This imports the new active_modified_distance.py script
    import active_modified_distance 
    
    # Set global vars in active_modified_distance based on args
    active_modified_distance.INPUT_CSV_PATH = args.input_csv
    active_modified_distance.PDB_FOLDER_PATH = args.pdb_folder
    active_modified_distance.OUTPUT_FILE = args.output_file
    
    active_modified_distance.process_residue_pairs() # Ensure this calls the correct main processing function

def main():
    parser = argparse.ArgumentParser(description='PPI Suite - combined tools for protein-protein interaction analysis')
    subparsers = parser.add_subparsers(dest='command', required=True)

    # SASA
    sasa_parser = subparsers.add_parser('sasa', help='Calculate solvent accessible surface area')
    sasa_parser.add_argument('-f', '--file', help='Single PDB file')
    sasa_parser.add_argument('-d', '--directory', help='Directory with PDB files')
    sasa_parser.add_argument('-o', '--output', default='.', help='Output directory')
    sasa_parser.add_argument('-r', '--residue', default='CYS', help='Residue type to analyze')
    sasa_parser.set_defaults(func=run_sasa)

    # PPI Interface Glutathionylation Analysis (uses ppi_interface_glutathionylation.py)
    glut_parser = subparsers.add_parser('ppi-glutathionylation', help='Analyze cysteine glutathionylation proximity to PPI interfaces. Input: CSV (uniprot,chain,residue) & PDB folder.')
    glut_parser.add_argument('input_csv', help='CSV file with uniprot, chain, residue columns for target cysteines')
    glut_parser.add_argument('pdb_folder', help='Folder containing PDB files of protein complexes')
    glut_parser.add_argument('output_csv', help='Output CSV file for results')
    glut_parser.add_argument('--threshold', type=float, default=None, help='Distance threshold (Å) for interface proximity check (e.g., 5.0)')
    glut_parser.set_defaults(func=run_ppi_glutathionylation_check)

    # Active Site to Modified Site Distance (Monomer, CSV input, uses active_modified_distance.py)
    am_parser = subparsers.add_parser('active-modified-distance', 
                                        help='Calculate intramolecular distances (e.g., active to modified sites) for monomers. \nInput: CSV with PDB name & residue pairs, and a folder of PDB/CIF files.',
                                        formatter_class=argparse.RawTextHelpFormatter)
    am_parser.add_argument('--input-csv', required=True, 
                               help='Path to CSV file defining PDBs and residue pairs.\nFormat: Column 1: PDB Name (no extension), Column 2: Active Site ResID, Column 3: Modified Site ResID.')
    am_parser.add_argument('--pdb-folder', required=True, help='Path to the folder containing PDB or mmCIF files.')
    am_parser.add_argument('--output-file', default="active_modified_distances.txt", 
                               help='Output file for results (default: active_modified_distances.txt)')
    am_parser.set_defaults(func=run_active_modified_distance)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()

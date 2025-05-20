import os
import glob
import pandas as pd
import sys
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB import Selection

def calculate_min_distance_sg(residue_a, chain_b_atoms, threshold=None):
    """
    Calculate the minimum Euclidean distance between the SG atom in residue A and all atoms in chain B.
    If residue A does not contain an SG atom, return None.
    If a threshold is provided, a "proximity" screening is first performed using NeighborSearch.
    """
    # Find SG atom in residue A
    sg_atom = None
    for atom in residue_a.get_atoms():
        if atom.get_name() == 'SG':
            sg_atom = atom
            break
    
    # If SG atom is not found, return None
    if sg_atom is None:
        return None
    
    if threshold is not None:
        # Use NeighborSearch to search only within the given threshold
        ns = NeighborSearch(list(chain_b_atoms))
        # Find all neighboring atoms within the threshold range
        neighbors = ns.search(sg_atom.coord, threshold)
        
        if len(neighbors) == 0:
            # If no atoms are found within the specified threshold, return None
            return None
        else:
            # Calculate the distance between the SG atom and all neighboring atoms
            distances = [sg_atom - atom_b for atom_b in neighbors]
            return min(distances)
    else:
        # No threshold limit, directly calculate the minimum distance from SG atom to all atoms in chain B
        distances = [sg_atom - atom_b for atom_b in chain_b_atoms]
        if distances:
            return min(distances)
        return None

def process_structures(input_csv, pdb_folder, output_csv, threshold=None):
    """
    Process each row in input_csv according to the following logic:
    1) Match all files with names matching fold_{uniprot}_*_model_0.pdb in the same folder;
    2) Parse each file, find the specified residue in chain A (chain in CSV);
    3) Calculate the minimum distance from the SG atom in that residue to chain B (optional threshold);
    4) Write the results to output_csv.
    """
    # Validate input parameters
    if not os.path.exists(input_csv):
        print(f"Error: Input CSV file not found: {input_csv}", file=sys.stderr)
        return False
        
    if not os.path.exists(pdb_folder):
        print(f"Error: PDB folder not found: {pdb_folder}", file=sys.stderr)
        return False
    
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        print(f"Error reading CSV file {input_csv}: {e}", file=sys.stderr)
        return False
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_csv)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
            print(f"Created output directory: {output_dir}")
        except OSError as e:
            print(f"Error creating directory {output_dir}: {e}", file=sys.stderr)
            return False
    
    results = []
    parser = PDBParser(QUIET=True)
    
    try:
        total_rows = len(df)
        processed_rows = 0
        
        for idx, row in df.iterrows():
            try:
                uniprot = str(row['uniprot']).strip().lower()
                chain_id = str(row['chain']).strip()  # According to your description, usually 'A'
                residue_id = int(row['residue'])      # Assuming residue is a number
            except (KeyError, ValueError) as e:
                print(f"Error processing row {idx}: {e}", file=sys.stderr)
                print(f"Row data: {row}")
                continue
            
            # Find all pdb files matching the pattern, e.g., {uniprot}.pdb
            pdb_pattern = os.path.join(pdb_folder, f"{uniprot}.pdb")
            pdb_files = glob.glob(pdb_pattern)
            
            if not pdb_files:
                print(f"[WARNING] No PDB files found matching: {pdb_pattern}")
                results.append({
                    'uniprot': uniprot,
                    'chain': chain_id,
                    'residue': residue_id,
                    'pdb_file': None,
                    'min_distance': None,
                    'has_sg': False
                })
                continue
            
            # Process all matching PDB files
            for pdb_file in pdb_files:
                structure_id = os.path.basename(pdb_file)
                
                try:
                    structure = parser.get_structure(structure_id, pdb_file)
                except Exception as e:
                    print(f"[ERROR] Failed to parse PDB file: {pdb_file}, Error: {e}")
                    continue
                
                # PDB structures may have multiple models, simply take the first model here
                model = structure[0]
                
                # Assuming there are only chains A and B
                # First get chain A, then locate the target residue
                if chain_id not in model:
                    print(f"[WARNING] Chain {chain_id} not found in {pdb_file}")
                    results.append({
                        'uniprot': uniprot,
                        'chain': chain_id,
                        'residue': residue_id,
                        'pdb_file': pdb_file,
                        'min_distance': None,
                        'has_sg': False
                    })
                    continue
                
                chain_a = model[chain_id]
                
                # Find the residue with the specified residue_id
                target_res = None
                for res in chain_a.get_residues():
                    if res.id[1] == residue_id:
                        target_res = res
                        break
                
                if target_res is None:
                    print(f"[WARNING] Residue {residue_id} not found in chain {chain_id} of {pdb_file}")
                    results.append({
                        'uniprot': uniprot,
                        'chain': chain_id,
                        'residue': residue_id,
                        'pdb_file': pdb_file,
                        'min_distance': None,
                        'has_sg': False
                    })
                    continue
                
                # Check if the residue contains an SG atom
                has_sg = False
                for atom in target_res.get_atoms():
                    if atom.get_name() == 'SG':
                        has_sg = True
                        break
                        
                if not has_sg:
                    print(f"[WARNING] SG atom not found in residue {residue_id} of chain {chain_id} in {pdb_file}")
                    results.append({
                        'uniprot': uniprot,
                        'chain': chain_id,
                        'residue': residue_id,
                        'pdb_file': pdb_file,
                        'min_distance': None,
                        'has_sg': False
                    })
                    continue
                
                # Get all atoms in chain B
                chain_b_id = 'B' if chain_id == 'A' else 'A'  # The problem states all files have only two chains, A and B
                if chain_b_id not in model:
                    print(f"[WARNING] Chain {chain_b_id} not found in {pdb_file}, cannot calculate distance")
                    results.append({
                        'uniprot': uniprot,
                        'chain': chain_id,
                        'residue': residue_id,
                        'pdb_file': pdb_file,
                        'min_distance': None,
                        'has_sg': has_sg
                    })
                    continue
                
                chain_b = model[chain_b_id]
                chain_b_atoms = list(chain_b.get_atoms())
                
                # Calculate the minimum distance from the SG atom to chain B
                min_dist = calculate_min_distance_sg(target_res, chain_b_atoms, threshold=threshold)
                
                results.append({
                    'uniprot': uniprot,
                    'chain': chain_id,
                    'residue': residue_id,
                    'pdb_file': pdb_file,
                    'min_distance': min_dist,
                    'has_sg': has_sg
                })
                
            processed_rows += 1
            # Print progress every 10 rows
            if processed_rows % 10 == 0 or processed_rows == total_rows:
                print(f"Progress: {processed_rows}/{total_rows} rows processed")
        
        # Output results to CSV
        out_df = pd.DataFrame(results)
        out_df.to_csv(output_csv, index=False)
        print(f"Calculation complete, results saved to: {output_csv}")
        return True
        
    except Exception as e:
        print(f"Error during processing: {e}", file=sys.stderr)
        return False


if __name__ == "__main__":
    # Command-line interface
    import argparse
    
    parser = argparse.ArgumentParser(description='Calculate minimum distances from SG atoms to chain atoms')
    parser.add_argument('input_csv', help='CSV file with uniprot, chain, residue columns')
    parser.add_argument('pdb_folder', help='Folder containing PDB files')
    parser.add_argument('output_csv', help='Output CSV file')
    parser.add_argument('--threshold', type=float, default=None, 
                        help='Distance threshold for neighbor search (optional)')
    
    args = parser.parse_args()
    
    success = process_structures(args.input_csv, args.pdb_folder, args.output_csv, threshold=args.threshold)
    if not success:
        sys.exit(1) 
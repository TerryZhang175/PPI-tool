# PPI-tool

Tools for Protein-Protein Interaction analysis.

## How to Install

1.  **Clone the repository:**
    ```bash
    git clone [<repository_url>](https://github.com/TerryZhang175/PPI-tool.git)
    cd PPI-tool
    ```
2.  **Create a virtual environment (recommended):**
    ```bash
    python -m venv .venv
    ```
    Activate the virtual environment:
    *   Windows:
        ```bash
        .\.venv\Scripts\activate
        ```
    *   macOS/Linux:
        ```bash
        source .venv/bin/activate
        ```
3.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```


## PPI Suite

`ppi_suite.py` provides a single command line entry point that exposes all functionality of this repository. The suite offers three subcommands:

* `sasa` – calculate the solvent accessible surface area of residues using `freesasa` (via the `sasa_tool.py` module).
* `ppi-glutathionylation` – analyze cysteine glutathionylation proximity to protein-protein interfaces (uses `ppi_interface_glutathionylation.py`).
* `active-modified-distance` – measure intramolecular distances in a monomer between specified active site and modified site residue pairs from a CSV file (uses `active_modified_distance.py`).

Run `python ppi_suite.py --help` for an overview of arguments.

## Contents

This repository contains tools for analyzing protein-protein interactions:

### SASA Tool (`sasa_tool.py`)

Calculates the Solvent Accessible Surface Area (SASA) for specified residues (or all residues) in protein structures using the FreeSASA library. SASA is a measure of how much a molecule's surface is accessible to a solvent and is often used to analyze protein folding, stability, and interactions.

This module can be run as a standalone script or used via the `ppi_suite.py sasa` command.

**Input:**

*   **Protein Structure File(s):**
    *   A single PDB file (e.g., `protein.pdb`) using the `-f` or `--file` argument.
    *   A directory containing multiple PDB files (e.g., `my_pdb_collection/`) using the `-d` or `--directory` argument. All files ending with `.pdb` (case-insensitive) in this directory will be processed.
*   **Residue Type (Optional):**
    *   Specify a three-letter amino acid code (e.g., `ALA`, `HIS`) using the `-r` or `--residue` argument.
    *   If not specified, defaults to `CYS` (Cysteine).
    *   To calculate SASA for *all* residue types, this functionality would need to be explicitly added to the script (e.g., by passing a special value like "ALL" or None and modifying `calculate_residue_sasa` to not filter by `residue_type`). *Currently, the script processes either a specific residue or defaults to CYS.*
*   **Output Directory (Optional):**
    *   Specify a directory path for the output files using the `-o` or `--output` argument.
    *   If not specified, defaults to `./sasa_results/single_pdb_sasa` for single files or `./sasa_results/<input_directory_name>_sasa` for directories. The script will create these directories if they don't exist.

**Output:**

*   For each processed PDB file, a text file (`<pdb_filename_prefix>_sasa.txt`) is generated in the specified output directory.
*   Each output file contains:
    *   The name of the PDB file.
    *   Lines detailing the SASA for each instance of the specified residue (or each CYS residue if defaulted). The format is: `RES_TYPE RES_NUM  SASA AREA_VALUE` (e.g., `CYS 42  SASA 15.32`).
    *   If a general residue type (not CYS) is specified, the format is `RES_TYPE RES_NUM SASA AREA_VALUE`.
    *   If CYS is specified (or defaulted), the format is `CYS RES_NUM SASA AREA_VALUE`.

**Atom Radii:**
The script uses a custom classifier with the following VDW radii for common atoms:
*   Nitrogen (N): 1.6 Å
*   Carbon (C): 1.7 Å
*   Oxygen (O): 1.4 Å
*   Sulfur (S): 1.8 Å
Unknown atoms are assigned a radius of 0.0 Å and ignored by FreeSASA.

**CLI Usage (via `ppi_suite.py`):**

*   Process a single PDB file for a specific residue (e.g., Alanine):
    ```bash
    python ppi_suite.py sasa -f path/to/protein.pdb -o path/to/output_dir -r ALA
    ```
*   Process all PDB files in a directory for Cysteine residues (default):
    ```bash
    python ppi_suite.py sasa -d path/to/pdb_directory -o path/to/output_dir
    ```
*   Process all PDB files in a directory for a specific residue (e.g., Leucine), output to default location:
    ```bash
    python ppi_suite.py sasa -d path/to/pdb_directory -r LEU
    ```

**Standalone Usage (`sasa_tool.py`):**

*   Process a single PDB file for Cysteine (default):
    ```bash
    python sasa_tool.py -f path/to/protein.pdb -o path/to/output_dir
    ```
*   Process all PDB files in a directory for a specific residue (e.g., Glycine):
    ```bash
    python sasa_tool.py -d path/to/pdb_directory -o path/to/output_dir -r GLY
    ```

- [SASA Tool Code](./sasa_tool.py)

### PPI Interface Glutathionylation Analysis (`ppi_interface_glutathionylation.py`)

Assesses if specific cysteine residues (potential glutathionylation sites) within a protein are located near a protein-protein interface in a complex. For each target cysteine, it identifies relevant protein complex structures and calculates the minimum distance from the cysteine's SG (sulfur) atom to any atom in other interacting protein chain(s). This helps determine if a glutathionylation site might be involved in or affected by protein interactions.

The script processes a list of target cysteines from an input CSV file and a folder of PDB files representing protein complexes.

**Input:**

1.  **Target Cysteines CSV File:**
    *   A CSV file specifying the UniProt ID, chain ID, and residue number of each cysteine to analyze.
    *   Required columns: `uniprot`, `chain`, `residue`.
    *   Example (`target_cysteines.csv`):
        ```csv
        uniprot,chain,residue
        P12345_P54321,A,42
        P12345_P13452,A,101
        Q67890_Q56567,B,78
        ```
2.  **PDB Folder:**
    *   A folder containing PDB files of protein complexes.
    *   The script searches for PDB files matching the pattern `fold_{uniprot}_*_model_0.pdb` within this folder, where `{uniprot}` is derived from the input CSV (converted to lowercase).
3.  **Output CSV File Path:**
    *   The path where the results CSV file will be saved.
4.  **Distance Threshold (Optional):**
    *   A floating-point value (in Angstroms) for proximity screening using the `--threshold` argument.
    *   If provided, the script first quickly checks for any atoms in the interacting chain(s) within this distance of the target cysteine's SG atom. If no atoms are found within this threshold, detailed distance calculation might be skipped or return `None` more quickly for distant interactions. If not provided, exact minimum distances are calculated to all atoms in the other chain.

**Output:**

*   A CSV file (specified by the `output_csv` argument) containing the analysis results.
*   Columns in the output CSV:
    *   `uniprot`: UniProt ID of the protein.
    *   `chain`: Chain ID of the target cysteine.
    *   `residue`: Residue number of the target cysteine.
    *   `pdb_file`: The PDB file from which the data was derived. If multiple PDB files match the pattern for a UniProt ID, there will be a row for each.
    *   `min_distance`: The minimum calculated distance (in Angstroms) from the SG atom of the target cysteine to any atom in the interacting chain (typically chain B if target is in A, or vice-versa). Will be `None` if the SG atom is not found, the other chain is not found, or no atoms are within the threshold (if specified).
    *   `has_sg`: Boolean (`True`/`False`) indicating if the target residue was found and contained an SG atom.

**Logic:**

*   For each entry in the input CSV, the script finds corresponding PDB files (e.g., for UniProt `P12345`, it looks for `fold_p12345*_model_0.pdb`).
*   It assumes a dimeric interaction (chains A and B). If the target cysteine is in chain A, it calculates distances to chain B, and vice-versa.
*   Only the first model in the PDB file is considered.

**CLI Usage (via `ppi_suite.py`):**

```bash
python ppi_suite.py ppi-glutathionylation --input_csv path/to/your_cysteine_list.csv --pdb_folder path/to/complex_pdb_files/ --output_csv path/to/glutathionylation_interface_results.csv --threshold 5.0
```

**Standalone Usage (`ppi_interface_glutathionylation.py`):**

```bash
python ppi_interface_glutathionylation.py path/to/your_cysteine_list.csv path/to/complex_pdb_files/ path/to/glutathionylation_interface_results.csv --threshold 5.0
```
If the threshold is not needed:
```bash
python ppi_interface_glutathionylation.py path/to/your_cysteine_list.csv path/to/complex_pdb_files/ path/to/glutathionylation_interface_results.csv
```

- [Code for `ppi_interface_glutathionylation.py`](./ppi_interface_glutathionylation.py)

### Active Site to Modified Site Distance Tool (Monomer) (`active_modified_distance.py`)

Calculates intramolecular distances between specified "active site" and "modified site" residues for one or more protein structures (monomers). This is useful for analyzing distances between, for example, catalytic residues and sites of post-translational modifications or mutations within the same protein chain.

It now processes a folder of PDB files based on an input CSV file.

**Input:**

1.  **PDB/Folder:** A folder containing your protein structure files (e.g., `protein1.pdb`).
2.  **Residue Pairs CSV File:** A CSV file that specifies the PDB file and the pairs of residues for distance calculation. The CSV file **must** have three columns. Headers are optional but recommended for clarity. The order of columns is important:
    *   **Column 1:** PDB File Name (the base name of the PDB file, without the `.pdb` extension. E.g., `protein1` if your file is `protein1.pdb`).
    *   **Column 2:** Active Site Residue Number (integer).
    *   **Column 3:** Modified Site Residue Number (integer).

    Example CSV (`residue_pairs.csv`):

    ```csv
    PDB_Name,ActiveSite_ResNum,ModifiedSite_ResNum
    my_protein_model_A,120,250
    my_protein_model_A,120,310
    another_protein,55,180
    ```

**Output:**

*   A tab-separated text file (default: `active_modified_distances.txt`) containing the following columns:
    *   `PDB_Identifier`: The base name of the PDB file from the CSV.
    *   `ActiveSite_ResNum`: Residue number of the active site.
    *   `ModifiedSite_ResNum`: Residue number of the modified site.
    *   `Distance(Å)`: The calculated distance in Angstroms between the specified atoms of the two residues.
    *   `ActiveSite_AtomType`: The type of atom used for the active site residue (CA, or SG for CYS).
    *   `ModifiedSite_AtomType`: The type of atom used for the modified site residue (CA, or SG for CYS).

**Atom Selection Logic:**

*   For Cysteine (CYS) residues, the script attempts to use the SG (sulfur) atom.
*   For all other residue types, it uses the CA (alpha carbon) atom.
*   The script currently assumes analysis of a monomer and considers only the first chain in the PDB file.

**CLI Usage (via `ppi_suite.py`):**

```bash
python ppi_suite.py active-modified-distance --input-csv path/to/your/residue_pairs.csv --pdb-folder path/to/your/pdb_files/ --output-file path/to/active_modified_distances.txt
```

**Standalone Usage (`active_modified_distance.py`):**

```bash
python active_modified_distance.py --input-csv path/to/your/residue_pairs.csv --pdb-folder path/to/your/pdb_files/ --output-file path/to/active_modified_distances.txt
```

This script is optimized for monomers and will use the first model and first chain found in each PDB/CIF file.

## Using the Command Line Interface

### SASA Calculation

Using the PPI Suite:
```bash
python ppi_suite.py sasa -f protein.pdb -o output_dir -r CYS
```
or
```bash
python ppi_suite.py sasa -d pdb_directory -o output_dir -r CYS
```

Standalone SASA tool usage:
```bash
python sasa_tool.py -f protein.pdb -o output_dir -r CYS
```

### PPI Interface Glutathionylation Analysis

Using the PPI Suite:
```bash
python ppi_suite.py ppi-glutathionylation --input_csv path/to/your_cysteine_list.csv --pdb_folder path/to/complex_pdb_files/ --output_csv path/to/glutathionylation_interface_results.csv --threshold 5.0
```
Or directly:
```bash
python ppi_interface_glutathionylation.py path/to/your_cysteine_list.csv path/to/complex_pdb_files/ path/to/glutathionylation_interface_results.csv --threshold 5.0
```

### Active Site to Modified Site Distance (Monomer)

Using the PPI Suite:
```bash
python ppi_suite.py active-modified-distance --input-csv path/to/your/residue_pairs.csv --pdb-folder path/to/your/pdb_files/ --output-file path/to/active_modified_distances.txt
```
Or directly:
```bash
python active_modified_distance.py --input-csv path/to/your/residue_pairs.csv --pdb-folder path/to/your/pdb_files/ --output-file path/to/active_modified_distances.txt
```

## Dependencies

- Python 3.6+
- Biopython
- pandas
- freesasa
- numpy

Install all requirements with:

```bash
pip install -r requirements.txt
```

## Citation

If you use these tools, please cite the original papers:

For SASA calculations this project relies on FreeSASA:
S. Mitternacht (2016)
FreeSASA: An open source C library for solvent accessible surface area.
F1000Research 5:189

## License

These tools are available under GPL V3.

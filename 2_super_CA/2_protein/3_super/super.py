import os
import mdtraj as md

def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def convert_residue_ids(residue_ids):
    # Convert residue IDs for chain 0 and chain 1
    converted_ids = []
    for rid in residue_ids:
        converted_ids.append(rid - 1)  # For chain 0
        converted_ids.append(rid + 98) # For chain 1
    return converted_ids

def calculate_rmsd_and_superimpose(dcd_path, pdb_path, reference_pdb_path, output_directory, residue_ids):
    create_directory_if_not_exists(output_directory)
    output_rmsd_path = os.path.join(output_directory, 'rmsd.txt')
    output_high_rmsd_path = os.path.join(output_directory, 'high_rmsd.txt')
    output_superimposed_pdb_path = os.path.join(output_directory, 'superimposed.pdb')

    try:
        traj = md.load(dcd_path, top=pdb_path)
        reference = md.load(reference_pdb_path)

        # Convert residue IDs to those used in the program
        converted_residue_ids = convert_residue_ids(residue_ids)

        # Select CA atoms of the specified residues for superimposition
        ca_atoms = traj.topology.select(f'name CA and (resid {" or resid ".join(map(str, converted_residue_ids))})')
        ca_atoms_ref = reference.topology.select(f'name CA and (resid {" or resid ".join(map(str, converted_residue_ids))})')

        rmsd_values = md.rmsd(traj, reference, atom_indices=ca_atoms, ref_atom_indices=ca_atoms_ref)
        rmsd_values_angstrom = rmsd_values * 10
        rmsd_output = [f"{i:6d} {rmsd:8.3f}" for i, rmsd in enumerate(rmsd_values_angstrom)]

        with open(output_rmsd_path, 'w') as f:
            f.write('\n'.join(rmsd_output))

        # Write high RMSD values to a separate file
        high_rmsd_output = [f"{i:6d} {rmsd:8.3f}" for i, rmsd in enumerate(rmsd_values_angstrom) if rmsd > 5.0]
        with open(output_high_rmsd_path, 'w') as f:
            f.write('\n'.join(high_rmsd_output))

        # Superimpose the trajectory onto the reference
        superimposed_traj = traj.superpose(reference, atom_indices=ca_atoms, ref_atom_indices=ca_atoms_ref)

        # Select all atoms of specified residues for output
        selected_atoms = traj.topology.select(f'(resid {" or resid ".join(map(str, converted_residue_ids))})')
        superimposed_traj = superimposed_traj.atom_slice(selected_atoms)

        # Save the superimposed coordinates
        superimposed_traj.save_pdb(output_superimposed_pdb_path)

    except Exception as e:
        print(f"Error processing RMSD calculation or superimposition: {e}")

def main():
    dcd_directory = '../2_pack/dcd'
    pdb_path = '../../../1_MD/2_protein_1/3ixo_group1/step4_equilibration.pdb'
    reference_pdb_path = '../../../1_MD/2_protein_1/3ixo_group1/step4_equilibration.pdb'
    output_directory = 'rmsd'
    residue_ids = [8, 23, 25, 27, 28, 29, 30, 32, 47, 48, 49, 50, 81, 82, 84]

    dcd_files = sorted([os.path.join(dcd_directory, f) for f in os.listdir(dcd_directory) if f.endswith('.dcd')])

    for dcd_file in dcd_files:
        print(f"Processing {dcd_file}")
        file_output_directory = os.path.join(output_directory, os.path.basename(dcd_file).replace('.dcd', ''))
        calculate_rmsd_and_superimpose(dcd_file, pdb_path, reference_pdb_path, file_output_directory, residue_ids)

if __name__ == "__main__":
    main()


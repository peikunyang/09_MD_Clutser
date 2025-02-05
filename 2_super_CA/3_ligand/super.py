import os
import mdtraj as md
import numpy as np
import gc
import tqdm

def save_as_pdb(selected_traj, output_file):
    selected_traj.save_pdb(output_file)

def process_trajectory(traj, reference, traj_atoms, reference_atoms, output_file):
    try:
        traj.superpose(reference, atom_indices=traj_atoms, ref_atom_indices=reference_atoms)
        chain_0_atoms = traj.topology.select("chainid 0")
        selected_traj = traj.atom_slice(chain_0_atoms)
        save_as_pdb(selected_traj, output_file)
    except Exception as e:
        print(f"處理軌跡時發生錯誤：{e}")

def Finish_MD():
    pdb_md = [f for f in sorted(os.listdir('../../1_MD/3_ligand/pdb')) 
              if os.path.exists(f'../../1_MD/3_ligand/pdb/{f}/step5.pdb')]
    return pdb_md

def Lig_Name(pdb_files):
    ligs = sorted(os.listdir('../1_pro-lig/3_seperate/2_lig/pdb'))
    lig_name = [lig for pdb_file in pdb_files for lig in ligs if pdb_file == lig[:4]]
    return lig_name

def Main():
    pdb_files = Finish_MD()
    output_dir = 'pdb'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    lig_name = Lig_Name(pdb_files)

    for i in tqdm.tqdm(range(len(pdb_files)), desc="Processing trajectories"):
        traj_file = f'../../1_MD/3_ligand/pdb/{pdb_files[i]}/step5.dcd'
        top_file = f'../../1_MD/3_ligand/pdb/{pdb_files[i]}/step4_equilibration.pdb'
        reference_file = f'../1_pro-lig/3_seperate/2_lig/pdb/{pdb_files[i][:4]}.pdb'
        
        try:
            reference = md.load(reference_file)
            traj = md.load(traj_file, top=top_file)
        except Exception as e:
            print(f"讀取文件時發生錯誤：{e}")
            continue

        traj_atoms = md.load_frame(top_file, index=0).topology.select("chainid 0 and not element H")
        traj_residues = [res for res in md.load_frame(top_file, index=0).topology.chain(0).residues]
        
        if len(traj_residues) != 1:
            print(f"軌跡 {pdb_files[i]} 中的 A 鏈不止一個非氫殘基。跳過。")
            continue

        traj_residue_name = traj_residues[0].name
        ref_chain_id = None
        for chain in reference.topology.chains:
            chain_residues = [res for res in chain.residues if res.name == traj_residue_name]
            if len(chain_residues) == 1:
                ref_chain_id = chain.index
                break

        if ref_chain_id is None:
            print(f"參考結構 {pdb_files[i]} 中找不到相同的單個殘基鏈。跳過。")
            continue

        reference_atoms = reference.topology.select(f"chainid {ref_chain_id} and not element H")
        
        if len(traj_atoms) > 0 and len(reference_atoms) > 0:
            print(f"選擇的原子數量: {len(traj_atoms)}")
            output_pdb_file = f'{output_dir}/{pdb_files[i]}.pdb'
            process_trajectory(traj, reference, traj_atoms, reference_atoms, output_pdb_file)
            print(f"成功處理並儲存 {pdb_files[i]} 的所有原子座標。")
        else:
            print(f"選擇的原子數量不匹配或為0。跳過 {pdb_files[i]}。")

        # 釋放記憶體
        del traj, reference, traj_atoms, reference_atoms
        traj = reference = traj_atoms = reference_atoms = None
        gc.collect()  # 強制進行垃圾回收

if __name__ == "__main__":
    Main()


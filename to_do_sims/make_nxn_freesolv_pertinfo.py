import glob
import csv
import re
import pandas as pd
import itertools
import numpy as np
from tqdm import tqdm
import csv
import matplotlib.pyplot as plt 
import seaborn as sns 

import sklearn
from sklearn.manifold import TSNE
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, rdFMCS, AllChem, rdmolfiles, Descriptors, rdchem

def constructSmarts(lig_mol, mcs_object):
    """
    Given a ligand and MCS (generated with a second ligand), construct an alternative SMARTS that contains
    information on the anchor atom (i.e. the atom in the MCS the perturbed R-group is attached to.)
    
    Get all neighbour indices of the fragment atoms in the original molecule.  
    The (single) index that is in the neighbour indices but not in the fragment 
    indices (set difference) is the atom we want. Anchor atoms and fragments are 
    in the same order because of consistent RDKit indexing.
    """
    # get the fragments by subtracting MCS from ligand.
    lig_fragments = Chem.ReplaceCore(lig_mol, Chem.MolFromSmarts(mcs_object.smartsString))
    
    
    # get the atom indices for the MCS object.
    mcs_indices = lig_mol.GetSubstructMatch(Chem.MolFromSmarts(mcs_object.smartsString))

    # get all the indices for the ligand.
    ligand_indices = set([x for x in range(0, lig_mol.GetNumAtoms())])

    # get all the fragment indices.
    non_mcs_indices = set(ligand_indices) - set(mcs_indices)


    new_smarts = None
    anchor_atoms = []

    for frag_idx in non_mcs_indices:
        # get the neighbours for this fragment atom.
        nghbrs = lig_mol.GetAtomWithIdx(frag_idx).GetNeighbors()

        for nghbr in nghbrs:
            # find the set difference.
            if not nghbr.GetIdx() in non_mcs_indices:
                anchor_atoms.append(lig_mol.GetAtomWithIdx(nghbr.GetIdx()).GetSmarts())

    if not lig_fragments:
    	return ""
    for anchor, molfrag in zip(anchor_atoms, Chem.GetMolFrags(lig_fragments, asMols=True, sanitizeFrags=False)):
        # clean up anchor. We really only care about aromatic vs non-aromatic etc.
        anchor = anchor.replace("@","").replace("[","").replace("]","")
    
        # for each fragment, we construct the smarts as [anchor*]atoms which fits with SMARTS logic.
        # Check https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
        # frag_smiles sometimes contains two dummy attachment points (in cases of ring fragments),
        # but for our purposes it should be OK to only denote the first dummy with the anchor.
        frag_smarts = Chem.MolToSmiles(molfrag)
        
        # we can simply slice the fragment smarts (splitting on "[" breaks in case of multiple dummies).
        # this slice would break if we have >9 fragments in the ligand!
        frag_smarts_anchored = "["+anchor+"*]"+frag_smarts[4:]

        # build the new SMARTS string. Insert a "." when there is already a fragment in the new string.
        if not new_smarts:
            new_smarts = frag_smarts_anchored
        else:
            new_smarts += "."+frag_smarts_anchored
    
    # sometimes the picked ligand is the actual MCS so there are no pert SMARTS.
    if not new_smarts:
        new_smarts = ""
        
    return new_smarts

if __name__ == "__main__":
	with open("nxn_freesolv_pertinfo.csv", "w") as writefile:
		writer = csv.writer(writefile)

		# construct nested list of NxN FreeSolv. 
		ligand_paths = glob.glob("freesolv/*.pdb")
		fully_connected = list(itertools.combinations(ligand_paths, 2))

		# iterate over possible perts.
		for pert in tqdm(fully_connected):
			ligA_name = pert[0].split("/")[-1].replace(".pdb","")
			ligB_name = pert[1].split("/")[-1].replace(".pdb","")
			pertname = ligA_name+"~"+ligB_name

			
			# Read in the molecules.
			ligA = Chem.rdmolfiles.MolFromPDBFile(pert[0])
			ligB = Chem.rdmolfiles.MolFromPDBFile(pert[1])
			mcs = rdFMCS.FindMCS([ligA, ligB], ringMatchesRingOnly=True, completeRingsOnly=True)

			# get the alternate perturbation SMARTS
			pert_smartsA = constructSmarts(ligA, mcs)
			pert_smartsB = constructSmarts(ligB, mcs)
			pert_smarts = pert_smartsA+"~"+pert_smartsB

			# write to file.
			writer.writerow([pertname, pert_smarts])








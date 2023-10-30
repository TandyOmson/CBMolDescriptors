""" Geometric descriptors functions for molecules 
Descriptors:
    - Eccentricity
    - Asphericity
    - Molecular weight
    - Molecular volume (vdW volume)
    - Surface area (TPSA topological polar surface area)
    - WHIM (weighted hollistic invariant molecular descriptor)
    - Atomnum
    - Heavy atomnum
    - Single bondnum
    - Double bondnum
    - Triple bondnum
    - Rotatable bondnum
    - Flexible torsion
    - Ringnum
    - Aromatic ringnum
    - LogP (octanol-water partition coefficient) (better to find hexadecane logL descriptors?)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def eccentricity(mol):
    return rdMolDescriptors.CalcEccentricity(mol)

def asphericity(mol):
    return rdMolDescriptors.CalcAsphericity(mol)

def mol_weight(mol):
    return rdMolDescriptors.CalcExactMolWt(mol)

def mol_volume(mol):
    return AllChem.ComputeMolVolume(mol)

def surface_area(mol):
    return rdMolDescriptors.CalcTPSA(mol)

def whim(mol):
    return rdMolDescriptors.CalcWHIM(mol)

def atomnum(mol):
    return mol.GetNumAtoms()

def heavy_atomnum(mol):
    return rdMolDescriptors.CalcNumHeavyAtoms(mol)

def single_bondnum(mol):
    return len([bond for bond in mol.GetBonds() if bond.GetBondTypeAsDouble()==1.0])

def double_bondnum(mol):
    return len([bond for bond in mol.GetBonds() if bond.GetBondTypeAsDouble()==2.0])

def triple_bondnum(mol):
    return len([bond for bond in mol.GetBonds() if bond.GetBondTypeAsDouble()==3.0])

def rotatable_bondnum(mol):
    return rdMolDescriptors.CalcNumRotatableBonds(mol)

def flexible_torsion(mol):
    """ This will need to have a criterion for what is considered a flexible torision, a list of torsions can be created with rdkit.Chem.TorsionFingerprints.CalculateTorsionLists
    """
    pass

def ringnum(mol):
    return rdMolDescriptors.CalcNumRings(mol)

def aromatic_ringnum(mol):
    return rdMolDescriptors.CalcNumAromaticRings(mol)

def logp(mol):
    return rdMolDescriptors.CalcCrippenDescriptors(mol)[0]

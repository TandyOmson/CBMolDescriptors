""" INPUT METHODS
 1. CSV file with SMILES strings and identifiers
 2. directory with .sdf files (identifiers are filenames, and will be globbed)
    in config file, these could both be in the inputpath field, the program can check if its a directory or file

    OUTPUT METHODS
 1. pkl with dataframe of descriptors and molecules
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from glob import glob

def standardize_mol(mol,smiles):
        """ Standardize molecule with rdkit
        """
        template = Chem.MolFromSmiles(smiles)
        newMol = AllChem.AssignBondOrdersFromTemplate(template,mol)
        newMol_H = Chem.AddHs(newMol,addCoords=True)
        Chem.SanitizeMol(newMol_H)
        Chem.AssignStereochemistry(newMol_H)
        Chem.SetAromaticity(newMol_H)
        Chem.SetHybridization(newMol_H)
        AllChem.ComputeGasteigerCharges(newMol_H)

        return newMol_H

def get_csv_smiles(inputpath):
    """ Get SMILES strings from CSV file
        Format: identifier<tab>SMILES 
     """
    df = pd.read_csv(inputpath, sep='\t', header=None)
    df.columns = ['identifier', 'SMILES']
    df.set_index('identifier', inplace=True)

    # Build SMILES strings into rdkit molecules
    mols = {}
    for smiles, identifier in zip(df['SMILES'],df.index):
        mol = standardize_mol(Chem.MolFromSmiles(smiles),smiles)

        # Embed 3D coordinates
        param = Chem.rdDistGeom.ETKDGv2()
        param.pruneRmsThresh = 0.2
        Chem.rdDistGeom.EmbedMolecule(mol,params=param)

        try:
            AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s")
        except:
            continue
    
        mols[identifier] = mol

    df["guestmol"] = pd.Series(mols)

    return df

def get_sdfs(inputpath):
    """ Get sdf stuctures from a directory as a list of rdkit molecules 
    """
    sdf_files = glob(inputpath + '/*.sdf')
    mols = {}
    smiles = {}
    for i in sdf_files:
        try:
            mol = Chem.SDMolSupplier(i)[0]
        except:
            mol = "NoRecord"

        mols[i[len(inputpath)+1:-4]] = mol
        smiles[i[len(inputpath)+1:-4]] = Chem.MolToSmiles(mol)

    df = pd.DataFrame.from_dict(mols,columns=["guestmol"], orient='index')
    df["SMILES"] = pd.Series(smiles)

    return df

# Function for writing a directory of molecules of a certain filetype, if the a directory with sdfs has not already been provided
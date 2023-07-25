""" Helper functions """

import tblite
from tblite.interface import Structure, Calculator
import numpy as np
from rdkit import Chem
from rdkit.Chem import PropertyMol
import subprocess as sp
import os

def get_flags_as_dict(config, section):
    """ Gets all flags from a section as a dictioanry
    """
    flags = {}
    for flag, value in config.items(section):
        flags[flag] = config.getboolean(section, flag)

    return flags

def read_opt_out(outfile):
    """ Reads properties from an xTB optimisation output into a dictionary
        Propterties must be read backwards from the end of the file
    """
    gen = (i.split() for i in reversed(open(outfile,"r").readlines()))

    propdict = {}

    for i in gen:
        if i:
            # HOMO-LUMO Gap
            if " ".join(i[:2]) == "| HOMO-LUMO":
                propdict["HL_GAP"] = float(i[3])

            if " ".join(i[:3]) == "| Geometry Summary":
                next(gen)
                next(gen)
                propdict["quadrupole"] = next(gen)[1:]

            if i[0] == "full:":
                propdict["dipole"] = i[1:]

            if " ".join(i[:3]) == ":: repulsion energy":
                propdict["repulsion"] = float(i[3]) * 627.509

            if " ".join(i[:3]) == ":: -> Gsolv":
                propdict["gsolv"] = float(i[3]) * 627.509

            if " ".join(i[:3]) == ":: -> dispersion":
                propdict["dispersion"] = float(i[3]) * 627.509

            if " ".join(i[:3]) == ":: SCC energy":
                propdict["SCC"] = float(i[3]) * 627.509

            if i[-1] == "(LUMO)":
                propdict["LUMO"] = float(i[2])

            if i[-1] == "(HOMO)":
                propdict["HOMO"] = float(i[3])

    return propdict

def xtb_opt(df):
    """ Runs xTB optimisations, reads results as properties of the molecule object
    """
    sp.run(["mkdir","xtb_temp"])
    os.chdir("./xtb_temp")
    for mol,mol_name in zip(df["guestmol"],df.index):
        Chem.MolToMolFile(mol,f"{mol_name}_rdkit.sdf")
        sp.run(["obabel","-isdf",f"{mol_name}_rdkit.sdf","-osdf","-O",f"{mol_name}.sdf"])

        # Don't forget to add write SASA to input file
        with open("xtb.inp","w") as fw:
            fw.write("$write\n  gbsa=true\n$end\n")
        sp.run(["xtb","--input","xtb.inp",f"{mol_name}.sdf","--opt","vtight","--alpb","water"],stdout=open(f"{mol_name}.out","w"))
        propdict = read_opt_out(f"{mol_name}.out")

        # Add calculations from xTB to mol properties
        for prop,value in propdict.items():
            if isinstance(value,float):
                mol.SetDoubleProp(prop,value)
            elif isinstance(value,list):
                for count,i in enumerate(value,1):
                    mol.SetDoubleProp(f"{prop}_{count}",float(i))
    
    os.chdir("./../")

    return df


def xtb_calc(df):
    """ Runs xTB calculations, adds all results as properties of the molecule object
        CURRENTLY NOT USED, SEE "xtb_opt" FUNCTION
    """
    for mol,mol_name in zip(df["guestmol"],df.index):
        positions=np.array(mol.GetConformer().GetPositions())
        numbers=np.array([atom.GetAtomicNum() for atom in mol.GetAtoms()])
        charge = Chem.GetFormalCharge(mol)
        calc = Calculator("GFN2-xTB", numbers, positions, charge=charge)
        res = calc.singlepoint()

        # Partial charges are added to each atom
        for i,j in zip(mol.GetAtoms(),res.get("charges")):
            i.SetDoubleProp("partial",j)
        
        # Energy is stored in properties of the first atom
        mol.GetAtoms()[0].SetDoubleProp("energy",float(res.get("energy")))

        # Dipole is stored in the first three atoms
        for i,j in zip(mol.GetAtoms(),res.get("dipole")):
            i.SetDoubleProp("dipole",j)

        # Quadrupole is stored in the first six atoms
        for i,j in zip(mol.GetAtoms(),res.get("quadrupole")):
            i.SetDoubleProp("quadrupole",j)

    return df

class df_item:
    """ Class for storing lists and arrays as single objects in a dataframe
        for ease of insertion and reading
    """
    def __init__(self, listorarray):
        self.contents = listorarray
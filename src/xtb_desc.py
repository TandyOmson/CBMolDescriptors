""" Requires xTB run
Descriptors:
    - Dipole moment
    - Quadrupole moment
    - xTB energy
    - LUMO
    - HOMO
    - HOMO-LUMO Gap
    - repulsion energy
    - solvent energy (using ALPB solvent model)
"""

def dipole_moment(mol):
    # Dipole moment has been stored in the properties of the first three atoms
    dipole_moment = []
    k=1
    while k<4:
        dipole_moment.append(mol.GetDoubleProp(f"dipole_{k}"))
        k+=1
    
    return dipole_moment

def quadrupole_moment(mol):
    # Quadrupole moment has been stored in the properties of the first six atoms
    quad_moment = []
    k=1
    while k<7:
        quad_moment.append(mol.GetDoubleProp(f"quadrupole_{k}"))
        k+=1
    
    return quad_moment

def energy(mol):
    return mol.GetDoubleProp("SCC")

def lumo(mol):
    return mol.GetDoubleProp("LUMO")

def homo(mol):
    return mol.GetDoubleProp("HOMO")

def hl_gap(mol):
    return mol.GetDoubleProp("HL_GAP")

def repulsion(mol):
    return mol.GetDoubleProp("repulsion")

def solv_en(mol):
    return mol.GetDoubleProp("gsolv")
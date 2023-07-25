""" Quantum Mechanical Descriptors
Dscriptors:
    - Polarizability (represented by molar refractivity Crippen calcualted - MR)
    - HOMO
    - LUMO
    - SASA (solvent accessible surface area)
    - Hydrophobicity
    - Radius of gyration
    - Chirality
"""

from mordred import Polarizability

def polarizability(mol):
    desc_instance =  Polarizability.APol()
    return desc_instance(mol)

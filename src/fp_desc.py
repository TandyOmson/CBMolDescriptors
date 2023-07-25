""" Fingerprint descriptors
Note that most of these are multi-component, and must have averages/selected values taken. Either that or special treatments.
Descriptors:
    - Chi0
    - Chi1 (etc?)
    - BCUT (eigenvalues of adjacency matrix, weighting diagonal elements with atoms weights)
    - BalabanJ (number of edges on molecular graph)
    - Barysz (weighted distance matrix accounting for heteroatoms and mulitple bonds)
    - Estrada Index (degree of folding of a protein - could be relevant for docked molecules?)
    - BertzT (Topological index that quantifies complexity)
    - Hall-Kier alpha (electrotopological state index)
    - Weinger index (counts the number of bonds between pairs of atoms and sums the distance between all pairs)
    - Zagreb index (first is sum of squares of the degrees of vertices, second is sum of the products of the degrees of pairs of adjacent vertices, etc.)
    - WalkCount
"""

from mordred import WalkCount,BCUT,Chi,BalabanJ,BaryszMatrix,BertzCT,EState

# note I might have to do more with these after calculating
def walkcount(mol):
    desc_instance =  WalkCount.WalkCount()
    return desc_instance(mol)

def bcut(mol):
    desc_instance = BCUT.BCUT()
    return desc_instance(mol)

def chi(mol):
    desc_instance = Chi.Chi()
    return desc_instance(mol)

def estate(mol):
    # Returns a list
    desc_instance = EState.EState.EState.EStateIndices(mol)
    return desc_instance

def balabanj(mol):
    desc_instance = BalabanJ.BalabanJ()
    return desc_instance(mol)

def baryszmatrix(mol):
    desc_instance = BaryszMatrix.BaryszMatrix()
    return desc_instance(mol)

def bertzct(mol):
    desc_instance = BertzCT.BertzCT()
    return desc_instance(mol)
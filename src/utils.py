#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import numpy as np
import MDAnalysis


def leafletdiv(atomgroup, divatom):
    print("Dividing selection using atomname \"%s\""%divatom)

    # get the average over all atoms
    u = atomgroup.universe
    atoms = atomgroup.residues.atoms.select_atoms("name %s"%divatom)
    zavg = np.mean(atoms.positions[:, 2])
    upper = MDAnalysis.AtomGroup([], u)
    lower = MDAnalysis.AtomGroup([], u)

    for a in atoms:
        if(a.position[2]>zavg):
            upper += a
        else:
            lower += a


    upper = atomgroup.residues.atoms.select_atoms("same residue as group upper", upper=upper) & atomgroup
    lower = atomgroup.residues.atoms.select_atoms("same residue as group lower", lower=lower) & atomgroup

    return upper, lower



def costheta2(vec):
    """
    Calculate the square of the cosine of the angle between the vector and z-axis
     Fortunately the dot product of (x1,x2,x3) and (0,0,1) is simply x3, so the
    calculation is rather trivial
     If vec is a numpy ndarray, then the result will be calculated over the first
    dimension, which should be of length 3 (any additional columns will be ignored
    and less columns will give an error)
    """
    return (vec[2]**2 / (vec[0]**2+vec[1]**2+vec[2]**2))


def get_CtoH_selections(carbons):
    print("Starting to find the hydrogen bonded to %d atoms"%len(carbons))
    u = carbons.universe
    sel1 = MDAnalysis.AtomGroup([], u)
    sel2 = MDAnalysis.AtomGroup([], u)
    for C in carbons:
        selH = C.residue.atoms.select_atoms("name H* and bonded group C", C=MDAnalysis.AtomGroup([C]))
        for a in selH:
            sel1 += C
            sel2 += a
    return sel1, sel2

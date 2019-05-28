#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import numpy as np
import MDAnalysis, sys

outputted1=False


def leafletdiv(atomgroup, divatom):
    global outputted1
    if (not outputted1):
        print("Dividing selection using atomname \"%s\""%divatom)
        outputted1 = True

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
    if(np.all([carbons[0].name==c.name for c in carbons])):
        print("Starting to find the hydrogen bonded to %d atoms (%s)"%(len(carbons), carbons[0].name))
    else:
        print("Starting to find the hydrogen bonded to %d atoms"%(len(carbons)))
    u = carbons.universe
    sel1 = MDAnalysis.AtomGroup([], u)
    sel2 = MDAnalysis.AtomGroup([], u)
    for C in carbons:
        selH = C.residue.atoms.select_atoms("name H* and bonded group C", C=MDAnalysis.AtomGroup([C]))
        for a in selH:
            sel1 += C
            sel2 += a


    print("Found %d atoms"%len(sel2))
    return sel1, sel2



def setup_selections(u, sel1, sel2, leaflets, leafletatom):
        # Get the selections
        print("\nSetting up selections")
        if(isinstance(sel1[0], str)):
            sel1 = [u.select_atoms(s) for s in sel1]
        else:
            sel1 = [MDAnalysis.AtomGroup(s, u) for s in sel1]

        if(len(sel1)>1):
            carbnames = [s.atoms[0].name for s in sel1]
        else:
            carbnames = [""]
        if(leaflets):
            sellower1 = [None for s in sel1]
            for i in range(len(sel1)):
                sel1[i], sellower1[i] = leafletdiv(sel1[i], leafletatom)

        if(sel2==[""]):
            sel2 = [None for s in sel1]
            for i in range(len(sel1)):
                sel1[i], sel2[i] = get_CtoH_selections(sel1[i])
            if(leaflets):
                sellower2 = [None for s in sellower1]
                for i in range(len(sel1)):
                    sellower1[i], sellower2[i] = get_CtoH_selections(sellower1[i])
        else:
            if(isinstance(sel2[0], str)):
                sel2 = u.select_atoms(sel2[0])
            else:
                sel2 = MDAnalysis.AtomGroup(sel2[0], u)
            if(leaflets):
                sellower2 = [None]
                sel2[0], sellower2[0] = leafletdiv(sel2[0], leafletatom)

        # Check the selection sizes
        for i in range(len(sel1)):
            if(len(sel1[i].atoms) != len(sel2[i].atoms)):
                raise ValueError("sel1[%d] and sel2[%d] are different sizes: %d and %d" % (i, i, len(sel1[i].atoms), len(sel2[i].atoms)))
            if(leaflets and len(sellower1[i].atoms) != len(sellower2[i].atoms)):
                raise ValueError("sellower1[%d] and sellower2[%d] are different sizes: %d and %d" % (i, i, len(sel1[i].atoms), len(sel2[i].atoms)))

        if(leaflets):
            print("Calculating order parameter for selections of %d atoms for upper leaflet"%(np.sum([len(s.atoms) for s in sel1])))
            print("Calculating order parameter for selections of %d atoms for lower leaflet"%(np.sum([len(s.atoms) for s in sel1])))
        else:
            print("Calculating order parameter for selections of %d atoms"%(np.sum([len(s.atoms) for s in sel1])))


        if leaflets:
            return sel1, sel2, sellower1, sellower2, carbnames

        return sel1, sel2, carbnames

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



def costheta2(vec, u):
    """
    Calculate the square of the cosine of the angle between the vector and u
     Vec should be a numpy ndarray, then the result will be calculated
    over the last dimension, which should be of length 3
    """
    return (np.dot(vec, u) / np.linalg.norm(vec, axis=-1))**2



def orderNoH(r, u):
    S  = np.full(r.shape[:-1], np.nan)
    # x1 is vec from Cn-1 to Cn
    x1 = r[1:-1]-r[:-2]
    # x2 is vec from Cn to Cn+1
    x2 = r[2:]  -r[1:-1]
    # Sz is vec from Cn-1 to Cn+1
    Sz = x1+x2
    # Sx = x1 x x2
    Sx = np.cross(x1, x2)
    # Sy = Sx x Sz
    Sy = np.cross(Sx, Sz)
    # normalize Sx and Sy
    Sx /= np.linalg.norm(Sx, axis=-1, keepdims=True)
    Sy /= np.linalg.norm(Sy, axis=-1, keepdims=True)
    # get the component of Sx and Sy orhogonal with the membrane normal
    sx = 1.5*np.dot(Sx, u)**2-0.5
    sy = 1.5*np.dot(Sy, u)**2-0.5
    # S = 1/3(2sx+sy)
    S[1:-1] = 2*sx/3+sy/3
    return S



def order(r1, r2=None, noH=True, u=(0, 0, 1)):
    """A function for calculating the order parameter
    params:
        with noH=False:
            r1: ndarray(...,3) of atomic positions C
            r2: ndarray(...,3) of matching coordinates of H
        with noH=True:
            r1: ndarray(m,n,3) of carbon coordinates of m carbons in n chains
            r2: ignored
        u: membrane normal unit vector
    """
    if(not noH):
        w = [None for r in r1]
        for i in range(len(r1)):
            w[i] = 1.5*costheta2(r2[i]-r1[i], u)-0.5
    else:
        w = orderNoH(r1, u)
    return w



def checkBondLengths(r):
    """
    A sanity check, that checks if the bonds are waay too long (for noH calculation), meaning there's a
    problem with the indexing of the group
    """

    # x1 is vec from Cn-1 to Cn
    x1 = r[1:-1]-r[:-2]
    # x2 is vec from Cn to Cn+1
    x2 = r[2:]  -r[1:-1]
    if(np.any(np.linalg.norm(x1, axis=-1)>3) or np.any(np.linalg.norm(x2, axis=-1)>3)):
        raise ValueError("A bond longer than 3 Å found, check your index file / -sel1")


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



def setup_selections(u, sel1, sel2, leaflets, leafletatom, thick, thickatom, noH):
    # Setup defaults to None
    sellower1 = None; sellower2 = None; thickup = None; thickdown=None

    # Get the selections
    print("\nSetting up selections")
    if(isinstance(sel1[0], str)):
        sel1 = [u.select_atoms(s) for s in sel1]
    else:
        sel1 = [MDAnalysis.AtomGroup(s, u) for s in sel1]

    # Thickness selections
    if(thick):
        thickup   = MDAnalysis.AtomGroup([], u)
        for s in sel1: thickup += s

        # Get a selection of the thickatom from each residue
        thickup = thickup.residues.atoms.select_atoms("name %s"%thickatom)
        # And then divide to leaflets
        thickup, thickdown = leafletdiv(thickup, leafletatom)

    # Get carbon names if sepcarbs
    if(len(sel1)>1):
        carbnames = [s.atoms[0].name for s in sel1]
    else:
        carbnames = [""]

    # Divide to leaflets if set
    if(leaflets):
        sellower1 = [None for s in sel1]
        for i in range(len(sel1)):
            sel1[i], sellower1[i] = leafletdiv(sel1[i], leafletatom)

    if(not noH):
        if(sel2==[""]):
            # If second selection isn't give, get the hydrogen
            sel2 = [None for s in sel1]
            for i in range(len(sel1)):
                sel1[i], sel2[i] = get_CtoH_selections(sel1[i])
            if(leaflets):
                sellower2 = [None for s in sellower1]
                for i in range(len(sel1)):
                    sellower1[i], sellower2[i] = get_CtoH_selections(sellower1[i])
        else:
            # Else set it up just like sel1
            if(isinstance(sel2[0], str)):
                sel2 = u.select_atoms(sel2[0])
            else:
                sel2 = MDAnalysis.AtomGroup(sel2[0], u)
            if(leaflets):
                sellower2 = [None]
                sel2[0], sellower2[0] = leafletdiv(sel2[0], leafletatom)


    print()

    # Check the selection sizes
    if(not noH):
        for i in range(len(sel1)):
            if(len(sel1[i].atoms) != len(sel2[i].atoms)):
                raise ValueError("sel1[%d] and sel2[%d] are different sizes: %d and %d" % (i, i, len(sel1[i].atoms), len(sel2[i].atoms)))
            if(leaflets and len(sellower1[i].atoms) != len(sellower2[i].atoms)):
                raise ValueError("sellower1[%d] and sellower2[%d] are different sizes: %d and %d" % (i, i, len(sellower1[i].atoms), len(sellower2[i].atoms)))
    else:
        for i, sel in enumerate(sel1):
            if(len(sel)!=len(sel1[0])):
                raise ValueError("-noH=True, but sel1[0] and sel1[%d] are different sizes: %d and %d" % (i, len(sel1[0].atoms), len(sel.atoms)))
        # Also check the bond lengths
        checkBondLengths(np.array([s.positions for s in sel1]))


    if(leaflets):
        if(noH):
            print("Calculating order parameter for %d selections of %d atoms for upper leaflet"%(len(sel1), len(sel1[0])))
            print("Calculating order parameter for %d selections of %d atoms for lower leaflet"%(len(sel1), len(sel1[0])))
        else:
            print("Calculating order parameter for selections of %d atoms for upper leaflet"%(np.sum([len(s.atoms) for s in sel1])))
            print("Calculating order parameter for selections of %d atoms for lower leaflet"%(np.sum([len(s.atoms) for s in sellower1])))
    else:
        if(noH):
            print("Calculating order parameter for %d selections of %d atoms"%(len(sel1), len(sel1[0])))
        else:
            print("Calculating order parameter for selections of %d atoms"%(np.sum([len(s.atoms) for s in sel1])))

    if(thick):
        print("Calculating thickness with a selection of %d atoms for upper leaflet"%(len(thickup.atoms)))
        print("Calculating thickness with a selection of %d atoms for lower leaflet"%(len(thickdown.atoms)))

    print()

    return sel1, sel2, sellower1, sellower2, thickup, thickdown, carbnames

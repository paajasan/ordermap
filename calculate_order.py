#!/usr/bin/env python3 -O
# -*- coding: UTF-8 -*-
import numpy as np
import MDAnalysis
import sys

import inout as io

def costheta2(vec):
    """
    Calculate the square of the cosine of the angle between the vector and z-axis
     Fortunately the dot product of (x1,x2,x3) and (0,0,1) is simply x3, so the
    calculation is rather trivial
    """
    return (vec[2]**2 / (vec[0]**2+vec[1]**2+vec[2]**2))


def get_CtoH_selections(carbons):
    #print("Starting to find the hydrogen bonded to %d atoms"%len(carbons))
    u = carbons.universe
    sel1 = MDAnalysis.AtomGroup([], u)
    sel2 = MDAnalysis.AtomGroup([], u)
    for C in carbons:
        selH = C.residue.atoms.select_atoms("name H* and bonded group C", C=MDAnalysis.AtomGroup([C]))
        for a in selH:
            sel1 += C
            sel2 += a
    return sel1, sel2





def calculate_order(topol, traj, sel1, sel2="None", dt=2500, ncells=20, center=None, out="order", plot=False, mindat=10):
    """
    A function that does all the magic.
    The parameters are pretty much simply those of the program itself, with the same defaults.
    """

    print("Setting up universe")
    # Start by loading a universe object
    u = MDAnalysis.Universe(topol, traj)
    print(u)
    frames = len(u.trajectory)
    print("Trajectory has %d frames of %g ps"%(frames, u.coord.dt))

    # Get the selections
    print("\nSetting up selections")
    sel1 = u.select_atoms(sel1)
    if(sel2==""):
        sel1, sel2 = get_CtoH_selections(sel1)
    else:
        sel2 = u.select_atoms(sel2)

    # Check the selection sizes
    if(len(sel1.atoms) != len(sel2.atoms)):
        raise ValueError("sel1 and sel2 are different sizes: %d and %d" % (len(sel1.atoms), len(sel2.atoms)))


    print("Calculating order parameter for selections of %d atoms"%(len(sel1.atoms)))


    # Get x and y values from the box vectors (assumes a cubic box)
    dim = u.dimensions
    x = np.linspace(0, dim[0], ncells)
    y = np.linspace(0, dim[1], ncells)

    # If we are centering, also get the centering group
    centering = False
    if(center != None):
        center = u.select_atoms(center)
        print("Using centering:\ncenter %d atoms"%len(center.atoms))
        centering = True
        com = center.center_of_mass()
        x -= com[0]
        y -= com[1]

    # Get the needed variables for calculating grid indices
    xmax = x[-1]
    xmin = x[0]
    ymax = y[-1]
    ymin = y[0]
    cellsizex = (xmax - xmin) / ncells
    cellsizey = (ymax - ymin) / ncells

    # Set up the grids and a variable for storing the satrttime of each dt
    prev = 0
    datagrid = np.zeros((x.shape[0], y.shape[0]))
    ngrid    = np.zeros(datagrid.shape)


    print("\nStarting to iterate trajectory")

    fromatstr = ("Frame %s%d%s"%("%",len(str(frames)), "d/%d"))

    print()
    # And finally start iterating the array
    for ts in u.trajectory:
        if(ts.frame%10==0):
            sys.stdout.write("\033[F\033[K") #Back o prev line and clear it
            print(fromatstr%(ts.frame, frames))

        t=ts.time

        # If this isn't the first frame and the modulo is zero, do stuff
        if(ts.frame % dt == 0 and ts.frame!=0):
            # Silence the warning about zero division, since we won't be using those values anyway
            with np.errstate(invalid='ignore'):
                datagrid /= ngrid


            datagrid *= 3/2
            datagrid -= 1/2
            # Write NaNs to where the ngrid is too small
            datagrid = np.where(ngrid < mindat, np.nan, datagrid)
            io.write_to_file(out, x, y, datagrid, prev, t, plot)
            prev = t

            datagrid = np.zeros(datagrid.shape)
            ngrid    = np.zeros(datagrid.shape)
            print()



        r1 = sel1.positions
        r2 = sel2.positions
        xycoord = r1[:, 0:2]
        vec = r2 - r1
        if(centering):
            xycoord -= center.center_of_mass()[0:2]

        for i, xy in enumerate(xycoord):
            xi = int((xy[0] - xmin) // cellsizex)
            yi = int((xy[1] - ymin) // cellsizey)

            if(xi >= 0 and yi >= 0 and xi < ncells and yi < ncells):
                datagrid[yi, xi] += costheta2(vec[i])
                ngrid[yi, xi]    += 1

    sys.stdout.write("\033[F\033[K") #Back o prev line and clear it
    print(fromatstr%(frames, frames))

    with np.errstate(invalid='ignore'):
        datagrid /= ngrid
    datagrid *= 3/2
    datagrid -= 1/2
    # for the last one we'll allow less data in case theere are less frames
    mindat = (mindat/dt)*(t-prev)/u.coord.dt
    datagrid = np.where(ngrid < mindat, np.nan, datagrid)
    io.write_to_file(out, x, y, datagrid, prev, t, plot)





if __name__ == '__main__':
    options = io.optP()
    calculate_order(
        options.struct_in, options.traj_in, options.sel1, options.sel2,
        options.dt, options.gridn, options.center, options.outfile, options.plot,
        options.mindat
    )

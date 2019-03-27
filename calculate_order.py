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


def processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, plot):
    with np.errstate(invalid='ignore'):
        datagrid /= ngrid

    datagrid *= 3/2
    datagrid -= 1/2
    # Write NaNs to where the ngrid is too small
    datagrid = np.where(ngrid < mindat, np.nan, datagrid)
    io.write_to_file(out, x, y, datagrid, prev, t, plot)





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
    x = np.linspace(0, dim[0], ncells+1)
    y = np.linspace(0, dim[1], ncells+1)

    # If we are centering, also get the centering group
    centering = False
    if(center != None):
        center = u.select_atoms(center)
        print("Using centering:\ncenter %d atoms"%len(center.atoms))
        centering = True
        com = center.center_of_mass()
        x -= com[0]
        y -= com[1]

    xmax = x[-1]
    xmin = x[0]
    ymax = y[-1]
    ymin = y[0]
    # move x and y from bin edges to bin centers
    x = x[:-1]+(x[1]-x[0])/2
    y = y[:-1]+(y[1]-y[0])/2


    # Set up the grids and a variable for storing the satrttime of each dt
    prev = 0
    datagrid = np.zeros((ncells, ncells))
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
            processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, plot)
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

        xcoord, ycoord = xycoord.T

        # A weighted (non normed) histogram is just the sums of the weights in each gridcell
        stat, x_edge, y_edge = np.histogram2d(xcoord, ycoord, weights=costheta2(vec.T),
                                                bins=ncells, range=((xmin, xmax), (ymin, ymax)))
        # And then we calculate the amount of points in each gridcell
        H,    x_edge, y_edge = np.histogram2d(xcoord, ycoord,
                                                bins=ncells, range=((xmin, xmax), (ymin, ymax)))

        datagrid += stat
        ngrid    += H


    # Print the last message
    sys.stdout.write("\033[F\033[K") #Back o prev line and clear it
    print(fromatstr%(frames, frames))

    # for the last one we'll allow less data in case theere are less frames
    mindat = (mindat/dt)*(t-prev)/u.coord.dt
    processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, plot)





if __name__ == '__main__':
    options = io.optP()
    calculate_order(
        options.struct_in, options.traj_in, options.sel1, options.sel2,
        options.dt, options.gridn, options.center, options.outfile, options.plot,
        options.mindat
    )

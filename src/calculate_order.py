#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import numpy as np
import MDAnalysis
import sys

from src.utils import *




def processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, leaflet, plot):
    with np.errstate(invalid='ignore'):
        datagrid /= ngrid

    datagrid *= 3/2
    datagrid -= 1/2
    # Write NaNs to where the ngrid is too small
    datagrid = np.where(ngrid < mindat, np.nan, datagrid)
    io.write_to_file(out, x, y, datagrid, prev, t, leaflet, plot)





def calculate_order(topol, traj, sel1, sel2="None", dt=2500, ncells=20, center=None, out="order", plot=False, mindat=10, leaflets=False, leafletatom="P*"):
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
    sel1leaf = ""
    if(leaflets):
        sel1leaf = "upper"
        sel1, sellower1 = leafletdiv(sel1, leafletatom)

    if(sel2==""):
        sel1, sel2 = get_CtoH_selections(sel1)
        if(leaflets):
            sellower1, sellower2 = get_CtoH_selections(sellower1)
    else:
        sel2 = u.select_atoms(sel2)
        sel2, sellower2 = leafletdiv(sel2, leafletatom)

    # Check the selection sizes
    if(len(sel1.atoms) != len(sel2.atoms)):
        raise ValueError("sel1 and sel2 are different sizes: %d and %d" % (len(sel1.atoms), len(sel2.atoms)))
    if(leaflets and len(sellower1.atoms) != len(sellower2.atoms)):
        raise ValueError("sellower1 and sellower2 are different sizes: %d and %d" % (len(sel1.atoms), len(sel2.atoms)))

    if(leaflets):
        print("Calculating order parameter for selections of %d atoms for upper leaflet"%(len(sel1.atoms)))
        print("Calculating order parameter for selections of %d atoms for lower leaflet"%(len(sellower1.atoms)))
    else:
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
    if(leaflets):
        datagridlow = np.zeros((ncells, ncells))
        ngridlow    = np.zeros(datagrid.shape)


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
            processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, sel1leaf, plot)
            prev = t

            datagrid = np.zeros(datagrid.shape)
            ngrid    = np.zeros(datagrid.shape)

            if(leaflets):
                processAndWrite(datagridlow, ngridlow, mindat, out, x, y, prev, t, "lower", plot)

                datagridlow = np.zeros(datagrid.shape)
                ngridlow    = np.zeros(datagrid.shape)

            print()



        r1 = sel1.positions
        r2 = sel2.positions
        xycoord = r1[:, 0:2]
        vec = r2 - r1
        if(centering):
            xycoord -= center.center_of_mass()[0:2]

        xcoord, ycoord = xycoord.T

        if(leaflets):
            r1low = sellower1.positions
            r2low = sellower2.positions
            xycoordlow = r1low[:, 0:2]
            veclow = r2low - r1low
            if(centering):
                xycoordlow -= center.center_of_mass()[0:2]

            xcoordlow, ycoordlow = xycoordlow.T

        # A weighted (non normed) histogram is just the sums of the weights in each gridcell
        stat, x_edge, y_edge = np.histogram2d(xcoord, ycoord, weights=costheta2(vec.T),
                                                bins=ncells, range=((xmin, xmax), (ymin, ymax)))
        # And then we calculate the amount of points in each gridcell
        H,    x_edge, y_edge = np.histogram2d(xcoord, ycoord,
                                                bins=ncells, range=((xmin, xmax), (ymin, ymax)))

        datagrid += stat
        ngrid    += H

        if(leaflets):
            stat, x_edge, y_edge = np.histogram2d(xcoordlow, ycoordlow, weights=costheta2(veclow.T),
                                                    bins=ncells, range=((xmin, xmax), (ymin, ymax)))
            H,    x_edge, y_edge = np.histogram2d(xcoordlow, ycoordlow,
                                                    bins=ncells, range=((xmin, xmax), (ymin, ymax)))

            datagridlow += stat
            ngridlow    += H


    # Print the last message
    sys.stdout.write("\033[F\033[K") #Back o prev line and clear it
    print(fromatstr%(frames, frames))

    # for the last one we'll allow less data in case theere are less frames
    mindat = (mindat/dt)*(t-prev)/u.coord.dt
    processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, sel1leaf, plot)

    if(leaflets):
        processAndWrite(datagridlow, ngridlow, mindat, out, x, y, prev, t, "lower", plot)






if __name__ == '__main__':

    import src.inout as io
    options = io.optP()
    calculate_order(
        options.struct_in, options.traj_in, options.sel1, options.sel2,
        options.dt, options.gridn, options.center, options.outfile, options.plot,
        options.mindat, options.leaflets, options.leafletatom
    )

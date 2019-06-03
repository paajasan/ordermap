#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import numpy as np
import MDAnalysis
import sys
import warnings

from src.utils import *
import src.inout as io




def processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, leaflet, plot, carbanames):
    with np.errstate(invalid='ignore'):
        datagrid /= ngrid

    # Write NaNs to where the ngrid is too small
    datagrid = np.where(ngrid < mindat, np.nan, datagrid)
    io.write_to_file(out, x, y, datagrid, prev, t, leaflet, plot, carbanames)





def calculate_order(topol, traj, sel1, sel2=[""], davg=2500, b=0, e=-1, dt=1, ncells=20, center=None,
                    out="order", plot=False, mindat=10, leaflets=False, leafletatom="P*", time=True,
                    thick=False, tout="thickness", thickatom="None", **kwargs):
    """
    A function that does all the magic.
    The parameters are pretty much simply those of the program itself, with the same defaults.
     To simplify the input between inout and this function, all keyword arguments not needed
    will be safely ignored through **kwargs, though it does raise a warning.
    """
    if(len(kwargs)>0):
        warnings.warn("Unrecognized keyword arguments: %s"%str(kwargs), Warning, stacklevel=2)

    print("Setting up universe")
    # Start by loading a universe object
    u = MDAnalysis.Universe(topol, traj)
    print(u)
    frames = len(u.trajectory)
    print("Trajectory has %d frames of %g ps"%(frames, u.coord.dt))

    if(thickatom=="None"): thickatom = leafletatom

    # Setup selections
    sel1leaf = "upper" if leaflets else ""

    sel1, sel2, sellower1, sellower2, thickup, thicklow, carbnames = setup_selections(u, sel1, sel2, leaflets, leafletatom, thick, thickatom)


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


    # Set up the grids and a variable for storing the satrttime of each davg
    prev = b*u.coord.dt
    datagrid = np.zeros((len(sel1), ncells, ncells))
    ngrid    = np.zeros(datagrid.shape)
    timedata = np.zeros((len(sel1), frames))
    tx       = np.zeros((len(sel1), frames))

    # dataholders for calculating the maps (see line 177)
    stat = np.empty(datagrid.shape)
    H    = np.empty(ngrid.shape)

    if(leaflets):
        datagridlow = np.zeros((len(sel1), ncells, ncells))
        ngridlow    = np.zeros(datagrid.shape)
        timedatalow = np.zeros((len(sel1), frames))
        statlow     = np.empty(datagridlow.shape)
        Hlow        = np.empty(ngridlow.shape)


    print("\nStarting to iterate trajectory")

    fromatstr = ("Frame %s%d%s"%("%",len(str(frames)), "d/%d"))

    print()


    # And finally start iterating the array
    e+=1
    if(e==0):
        e=None

    ############################################################################
    #     Trajectory iteration start                                           #
    ############################################################################

    for ts in u.trajectory[b:e:dt]:
        frame = ts.frame
        if(ts.frame%(10*dt)==0):
            sys.stdout.write("\033[F\033[K") #Back to prev line and clear it
            print(fromatstr%(frame, frames))

        t=ts.time

        # If this isn't the first frame and the modulo is zero, do stuff
        if(davg>0 and (frame-b) % davg == 0 and frame!=b):
            processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, sel1leaf, plot, carbnames)

            datagrid = np.zeros(datagrid.shape)
            ngrid    = np.zeros(datagrid.shape)

            if(leaflets):
                processAndWrite(datagridlow, ngridlow, mindat, out, x, y, prev, t, "lower", plot, carbnames)

                datagridlow = np.zeros(datagrid.shape)
                ngridlow    = np.zeros(datagrid.shape)

            prev =  t
            print()



        r1 = [s.positions for s in sel1]
        r2 = [s.positions for s in sel2]
        xycoord = [r[:, 0:2] for r in r1]
        vec = [R2-R1 for R1, R2 in zip(r1, r2)]
        if(centering):
            centerCom = center.center_of_mass()[0:2]
            for i in range(len(xycoord)):
                xycoord[i] -= centerCom

        xcoord, ycoord = ([z[:,0] for z in xycoord], [z[:,1] for z in xycoord])

        if(leaflets):
            r1low = [s.positions for s in sellower1]
            r2low = [s.positions for s in sellower2]
            xycoordlow = [r[:, 0:2] for r in r1low]
            veclow = [R2-R1 for R1, R2 in zip(r1low, r2low)]
            if(centering):
                for i in range(len(xycoordlow)):
                    xycoordlow[i] -= centerCom


            xcoordlow, ycoordlow = ([z[:,0] for z in xycoordlow], [z[:,1] for z in xycoordlow])

        w = [np.zeros(xx.shape) for xx in xcoord]

        for i in range(H.shape[0]):
            w[i] = 1.5*costheta2(vec[i].T)-0.5
            # A weighted (non normed) histogram is just the sums of the weights in each gridcell
            stat[i], x_edge, y_edge = np.histogram2d(xcoord[i], ycoord[i], weights=w[i],
                                                    bins=ncells, range=((xmin, xmax), (ymin, ymax)))
            # And then we calculate the amount of points in each gridcell
            H[i],    x_edge, y_edge = np.histogram2d(xcoord[i], ycoord[i],
                                                    bins=ncells, range=((xmin, xmax), (ymin, ymax)))


        if(time):
            timedata[:, ts.frame] = [np.mean(wi) for wi in w]
            tx[:, ts.frame]       = t

        datagrid += stat
        ngrid    += H

        if(leaflets):
            w = [np.zeros(xx.shape) for xx in xcoordlow]

            for i in range(Hlow.shape[0]):
                w[i] = 1.5*costheta2(veclow[i].T) - 0.5
                statlow[i], x_edge, y_edge = np.histogram2d(xcoordlow[i], ycoordlow[i], weights=w[i],
                                                        bins=ncells, range=((xmin, xmax), (ymin, ymax)))
                Hlow[i],    x_edge, y_edge = np.histogram2d(xcoordlow[i], ycoordlow[i],
                                                        bins=ncells, range=((xmin, xmax), (ymin, ymax)))

            if(time):
                timedatalow[:, ts.frame] = [np.mean(wi) for wi in w]

            datagridlow += statlow
            ngridlow    += Hlow


    ############################################################################
    #     Trajectory iteration ends                                            #
    ############################################################################


    # Print the last message
    sys.stdout.write("\033[F\033[K") #Back o prev line and clear it
    print(fromatstr%(frame, frames))


    if(davg>0 and (frame-b) % davg == 0):
        if(time):
            io.write_time_series(out, tx[:, b:e], timedata[:, b:e], sel1leaf, plot, carbnames)
            if(leaflets):
                io.write_time_series(out, tx[:, b:e], timedatalow[:, b:e], "lower", plot, carbnames)
        return

    # for the last one we'll allow less data in case there are less frames
    if(davg>0):
        mindat = (mindat/davg)*(t-prev)/u.coord.dt
    else:
        t=-1
    processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, sel1leaf, plot, carbnames)
    if(time):
        io.write_time_series(out, tx[:, b:e], timedata[:, b:e], sel1leaf, plot, carbnames)

    if(leaflets):
        processAndWrite(datagridlow, ngridlow, mindat, out, x, y, prev, t, "lower", plot, carbnames)
        if(time):
            io.write_time_series(out, tx[:, b:e], timedatalow[:, b:e], "lower", plot, carbnames)

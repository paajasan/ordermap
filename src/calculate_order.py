#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import numpy as np
import MDAnalysis
import sys
import warnings

from src.utils import *
import src.inout as io




def processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, leaflet, plot, carbanames, sepcar=False):
    with np.errstate(invalid='ignore'):
        if(not sepcar and datagrid.shape[0]!=1):
            datagrid = np.sum(datagrid, axis=0).reshape(([1]+list(datagrid.shape[1:])))
            ngrid    = np.sum(ngrid,    axis=0).reshape(([1]+list(datagrid.shape[1:])))
            carbanames = [""]
        datagrid /= ngrid

    # Write NaNs to where the ngrid is too small
    datagrid = np.where(ngrid < mindat, np.nan, datagrid)
    io.write_to_file(out, x, y, datagrid, prev, t, leaflet, plot, carbanames)





def calculate_order(topol, traj, sel1, sel2=[""], noH=False, davg=2500, b=0, e=-1, dt=1, ncells=20, center=None,
                    out="order", plot=False, mindat=10, leaflets=False, leafletatom="P*", time=True, thick=False,
                    tout="thickness", thickatom="None", u=(0,0,1), sepcar=False, unsatInd=None, **kwargs):
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
    univ = MDAnalysis.Universe(topol, traj)
    print(univ)
    frames = len(univ.trajectory)
    print("Trajectory has %d frames of %g ps"%(frames, univ.coord.dt))

    if(thickatom=="None"): thickatom = leafletatom

    # Setup selections
    sel1leaf = "upper" if leaflets else ""

    sel1, sel2, sellower1, sellower2, thickup, thicklow, carbnames = setup_selections(univ, sel1, sel2, leaflets, leafletatom, thick, thickatom, noH)


    # Get x and y values from the box vectors
    dim = univ.dimensions
    xmin,xmax,ymin,ymax = xy_extent(dim)
    x = np.linspace(xmin, xmax, ncells+1)
    y = np.linspace(ymin, ymax, ncells+1)

    # If we are centering, also get the centering group
    centering = False
    if(center != None):
        center = univ.select_atoms(center)
        print("Using centering:\ncenter %d atoms"%len(center.atoms))
        centering = True
        com = center.center_of_mass()
        x -= com[0]
        y -= com[1]

    # move x and y from bin edges to bin centers
    x = x[:-1]+(x[1]-x[0])/2
    y = y[:-1]+(y[1]-y[0])/2


    # Set up the grids and a variable for storing the satrttime of each davg
    prev = b*univ.coord.dt
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

    if(thick):
        # The first dimension is size 1, so that the array works with io.write_to_file() and processAndWrite()
        thickdata = np.zeros((1, ncells, ncells))
        thickn    = np.zeros(thickdata.shape)
        thicktimedata = np.zeros((1, frames))
        thicktx       = np.zeros((1, frames))


    print("\nStarting to iterate trajectory")

    fromatstr = ("Frame %s%d%s"%("%",len(str(frames)), "d/%d"))

    print()


    # And finally start iterating the array
    e+=1 # make e inclusive
    if(e==0):
        e=None # if e is None, all the frames are used

    ############################################################################
    #     Trajectory iteration start                                           #
    ############################################################################

    for ts in univ.trajectory[b:e:dt]:
        frame = ts.frame
        sys.stdout.write("\033[F\033[K") #Back to prev line and clear it
        print(fromatstr%(frame, frames))

        t=ts.time

        # If this isn't the first frame and the modulo is zero, do stuff
        if(davg>0 and (frame-b) % davg == 0 and frame!=b):
            processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, sel1leaf, plot, carbnames, sepcar)

            datagrid = np.zeros(datagrid.shape)
            ngrid    = np.zeros(datagrid.shape)

            if(leaflets):
                processAndWrite(datagridlow, ngridlow, mindat, out, x, y, prev, t, "lower", plot, carbnames)

                datagridlow = np.zeros(datagrid.shape)
                ngridlow    = np.zeros(datagrid.shape)


            if(thick):
                processAndWrite(thickdata, thickn, 0, tout, x, y, prev, t, "", "thicc"*plot, [""])

                thickdata = np.zeros(thickdata.shape)
                thickn    = np.zeros(thickdata.shape)

            prev =  t
            print()



        r1 = [s.positions for s in sel1]
        if(noH):
            r1 = np.asarray(r1)
        r2 = None
        if(not noH):
            r2 = [s.positions for s in sel2]
            vec = [R2-R1 for R1, R2 in zip(r1, r2)]

        xycoord = [r[:, 0:2].copy() for r in r1]
        if(centering):
            centerCom = center.center_of_mass()[0:2]
            for i in range(len(xycoord)):
                xycoord[i] -= centerCom

        xcoord, ycoord = ([z[:,0] for z in xycoord], [z[:,1] for z in xycoord])

        if(leaflets):
            r1low = [s.positions for s in sellower1]
            if(noH):
                r1low = np.asarray(r1low)
            r2low = None
            if(not noH):
                r2low = [s.positions for s in sellower2]
                veclow = [R2-R1 for R1, R2 in zip(r1low, r2low)]

            xycoordlow = [r[:, 0:2].copy() for r in r1low]
            if(centering):
                for i in range(len(xycoordlow)):
                    xycoordlow[i] -= centerCom


            xcoordlow, ycoordlow = ([z[:,0] for z in xycoordlow], [z[:,1] for z in xycoordlow])

        w = order(r1, r2, unsatInd, noH, u)

        for i in range(H.shape[0]):
            # A weighted (non normed) histogram is just the sums of the weights in each gridcell
            stat[i], x_edge, y_edge = np.histogram2d(xcoord[i], ycoord[i], weights=w[i],
                                                    bins=ncells, range=((xmin, xmax), (ymin, ymax)))
            # And then we calculate the amount of points in each gridcell
            H[i],    x_edge, y_edge = np.histogram2d(xcoord[i], ycoord[i],
                                                    bins=ncells, range=((xmin, xmax), (ymin, ymax)))


        if(thick):
            upperpos = thickup.positions
            lowerpos = thicklow.positions
            if(centering):
                upperpos[:, :2] -= centerCom
                lowerpos[:, :2] -= centerCom

            upperheight, x_edge, y_edge  = np.histogram2d(upperpos[:, 0],upperpos[:, 1], weights=upperpos[:, 2], bins=ncells, range=((xmin, xmax), (ymin, ymax)))
            upperdat,    x_edge, y_edge  = np.histogram2d(upperpos[:, 0],upperpos[:, 1], bins=ncells, range=((xmin, xmax), (ymin, ymax)))

            lowerheight, x_edge, y_edge  = np.histogram2d(lowerpos[:, 0],lowerpos[:, 1], weights=lowerpos[:, 2], bins=ncells, range=((xmin, xmax), (ymin, ymax)))
            lowerdat,    x_edge, y_edge  = np.histogram2d(lowerpos[:, 0],lowerpos[:, 1], bins=ncells, range=((xmin, xmax), (ymin, ymax)))

            # Make a "mask" for the data, that is only true when both leaflets have data
            hasData  = (upperdat!=0)*(lowerdat!=0)  # multiplying boolean arrays is the same as elementwise "and"

            # remove zeros (we want the gridcells with no data to have 0 instead of nan, so that averaging works)
            upperdat += upperdat==0
            lowerdat += lowerdat==0

            upperheight /= upperdat
            lowerheight /= lowerdat

            # the data       (upper    -   lower) * (0 where upper is 0 or lower is 0, else 1)      (upp==0 or low==0) is same as (upp!=0 and low!=0)
            thickdata[0] += (upperheight-lowerheight)*hasData
            thickn       += hasData



        if(time):
            timedata[:, ts.frame] = [np.mean(wi) for wi in w]
            tx[:, ts.frame]       = t
            if(thick):
                thicktimedata[0,ts.frame] = np.mean(upperpos[:, 2])-np.mean(lowerpos[:, 2])
                thicktx[0,ts.frame]       = t


        datagrid += stat
        ngrid    += H

        if(leaflets):
            w = order(r1low, r2low, unsatInd, noH, u)

            for i in range(Hlow.shape[0]):
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

    if((not sepcar) and timedata.shape[0]>1):
        timedata  = np.mean(timedata, axis=0).reshape(([1]+list(timedata.shape[1:])))
        carbnames = [""]

    if(davg>0 and (frame-b) % davg == 0):
        if(time):
            io.write_time_series(out, tx[:, b:e], timedata[:, b:e], sel1leaf, plot, carbnames)
            if(thick):
                io.write_time_series(tout, thicktx[:, b:e], thicktimedata[:, b:e], "", plot, [""])
            if(leaflets):
                io.write_time_series(out, tx[:, b:e], timedatalow[:, b:e], "lower", plot, carbnames)
        return

    # for the last one we'll allow less data in case there are less frames
    if(davg>0):
        mindat = (mindat/davg)*(t-prev)/univ.coord.dt
    else:
        t=-1
    processAndWrite(datagrid, ngrid, mindat, out, x, y, prev, t, sel1leaf, plot, carbnames, sepcar)

    if(thick):
        processAndWrite(thickdata, thickn, 0, tout, x, y, prev, t, "", "thicc"*plot, [""], sepcar)


    if(time):
        io.write_time_series(out, tx[:, b:e], timedata[:, b:e], sel1leaf, plot, carbnames)
        if(thick):
            io.write_time_series(tout, thicktx[:, b:e], thicktimedata[:, b:e], "", plot, [""])

    if(leaflets):
        processAndWrite(datagridlow, ngridlow, mindat, out, x, y, prev, t, "lower", plot, carbnames)
        if(time):
            io.write_time_series(out, tx[:, b:e], timedatalow[:, b:e], "lower", plot, carbnames)

#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def getxy(filename):
    with open(filename) as f:
        f.readline(); f.readline()
        x = f.readline().split()[2:]
        y = f.readline().split()[2:]
        for i in range(len(x)):
            x[i] = float(x[i])
        for i in range(len(y)):
            y[i] = float(y[i])
    return np.asarray(x), np.asarray(y)


def plot(filename, name="Order parameter", c=-1, levels=None):
    print("Reading data from "+filename)
    Z = np.loadtxt(filename, comments=["@","#"])
    Z = np.ma.masked_invalid(Z, copy=False)
    x, y = getxy(filename)

    if(np.all(Z.mask)):
        print("data in %s only contains NaN's, moving on"%filename)
        return

    ax = plt.gca()

    im = ax.contourf(x, y, c*Z, corner_mask=True, cmap="plasma", levels=levels)

    # Create colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)

    ax.set_title(name)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    print("Saving figure %s"%(filename[:-4]+".png"))
    plt.gcf().savefig(filename[:-4]+".png")
    plt.clf()


def plot_time_series(fname, name="Order parameter in time"):
    print("Reading data from "+fname)
    x, y = np.loadtxt(fname, unpack=True)

    plt.plot(x/1000, -y)

    plt.xlabel("t (ns)")
    plt.ylabel("order parameter")
    plt.title(name)
    print("Saving figure %s"%(fname[:-4]+"_ts.png"))
    plt.gcf().savefig(fname[:-4]+"_ts.png")
    plt.clf()



def plot_carbs(fname, name="Order parameter per carbon"):
    print("Reading data from "+fname)
    carbs, carbnums, y = ([], [], [])

    with open(fname) as f:
        for line in f:
            parts = line.split()
            if (line.startswith("#") or len(parts)!=2):
                continue
            carbs.append(parts[0])
            carbnums.append(parts[0][1:])
            y.append(-float(parts[1]))


    plt.plot(carbnums, y, "b-")
    plt.gca().set_xticks(carbnums)
    plt.gca().set_xticklabels(carbs)

    plt.xlabel("carbon")
    plt.ylabel("order parameter")
    plt.title(name)
    print("Saving figure %s"%(fname[:-4]+".png"))
    plt.gcf().savefig(fname[:-4]+".png")
    plt.clf()




if __name__ == '__main__':
    import os

    # get current working directory
    if (not os.getcwd().endswith("/analysis")):
        print("Not in a directory named \"analysis\", quitting...")
        os.exit(1)

    for f in os.listdir("order"):
        if(os.path.isfile("order/"+f) and f.endswith(".dat")):
            plot("order/"+f)

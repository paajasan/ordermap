# -*- coding: UTF-8 -*-
import argparse
import numpy as np

from src.create_figs import plot as plotdata, plot_time_series as plot_ts, plot_carbs as plot_c

def getSelection(sel, sepcar):
    if(sepcar):
        sellist = []
    parts = sel.split()
    if(len(parts)!=3):
        return [sel]
    rangeparts = parts[2].split("-")
    if(len(rangeparts)==2):
        try:
            start = int(rangeparts[0])
            end   = int(rangeparts[1])
        except (ValueError):
            return [sel]

        sel = "resname %s and name "%parts[0]
        for i in range(start, end+1):
            if(sepcar):
                sellist.append(sel+"%s%d"%(parts[1], i))
            else:
                sel += "%s%d "%(parts[1], i)
    elif(len(rangeparts)==1):
        try:
            start = int(rangeparts[0])
        except (ValueError):
            return [sel]
        sel = "resname %s and name %s%d "%(parts[0], parts[1], start)
        if(sepcar):
            sellist.append(sel)
    else:
        return [sel]

    if(sepcar):
        return sellist
    return [sel[:-1]]



def optP():

    description = "A program for calculating deuterium order parameter heatmaps"
    optParser = argparse.ArgumentParser(prog="ordermap", description=description)

    optParser.add_argument(
        "-V", '--version',
        action="version",
        version="0.2.0_dev"
    )


    optParser.add_argument(
        '-f', type=str,
        dest='traj_in',
        metavar="inputfilename",
        required=True,
        help="Input trajectory (xtc)"
    )

    optParser.add_argument(
        '-s', type=str,
        dest='struct_in',
        metavar="inputfilename",
        required=True,
        help="Input structure (gro/pdb/tpr)"
    )

    optParser.add_argument(
        '-sel1', type=str,
        dest='sel1',
        metavar="selection_str",
        required=True,
        help="A selection string to find the atoms from which the vector for the order parameter starts from. "\
            "Either a selection string understood by MDAnalysis or to for example get all POPC carbons with atomname C20 to C30 \"POPC C 20-30\"."
    )

    optParser.add_argument(
        '-sel2', type=str,
        dest='sel2',
        metavar="selection_str",
        default="",
        help="A selection string to find the atoms on which the vector for the order parameter ends" \
            " on (must result in an equal sized selection as --sel1). If not given, all the vectors from"\
            " each atom in sel1 to all hydrogen bonded to it, will be used"
    )

    optParser.add_argument(
        '-o',
        dest='outfile',
        metavar="outputdest",
        default="order.dat",
        help="The output filenames without the times. \"-o out.dat\" and \"-o out\" will both result in " \
            "filenaming out0to250.dat (and so on) [Default: %(default)s]"
    )

    optParser.add_argument(
        '-b', type=int,
        dest='b',
        metavar="N",
        default=0,
        help="The first frame to use for analysis [Default: %(default)d]"
    )



    optParser.add_argument(
        '-e', type=int,
        dest='e',
        metavar="N",
        default=-1,
        help="The last frame to use for analysis. Use negative values to count from the end (like with python indices usually, but inclusive, ie with -1 use the last frame too) [Default: %(default)d]"
    )

    optParser.add_argument(
        '-dt', type=int,
        dest='dt',
        metavar="N",
        default=1,
        help="Only use every N frames [Default: %(default)d]"
    )

    optParser.add_argument(
        '-center',
        dest='center',
        metavar="selection_str",
        default=None,
        help="A selection string to define a group, which will be used as a center for the grids (located at (0, 0))"
    )

    optParser.add_argument(
        '-davg', type=int,
        dest='davg',
        metavar="N",
        default=-1,
        help="The amount of frames to average the order parameters angle over [Default: %(default)d]"
    )

    optParser.add_argument(
        '-gridn', type=int,
        dest='gridn',
        metavar="N",
        default=20,
        help="The amount of gridpoints per side (ie. with gridn=N we end up with an NxN grid) [Default: %(default)d]"
    )

    optParser.add_argument(
        '-mindat', type=int,
        dest='mindat',
        metavar="N",
        default=10,
        help="The minimum amount data in a gridcell (cells with Ndata<gridn will have NaN) [Default: %(default)d]"
    )

    optParser.add_argument(
        '-plot', action="store_true",
        dest='plot',
        help="Plots the grids if specified"
    )

    optParser.add_argument(
        '-leafdiv', action="store_true",
        dest='leaflets',
        help="Does the calculation on both leaflets separetely"
    )

    optParser.add_argument(
        '-divatom',
        dest='leafletatom',
        default="P*",
        help="The atomname used for the division of leaflets [Default: %(default)s]"
    )


    optParser.add_argument(
        '-timeseries', action="store_true",
        dest='time',
        help="Also outputs the orderparameter for each timestep"
    )


    optParser.add_argument(
        '-sepcarbs', action="store_true",
        dest='sepcar',
        help="Calculate everything separately for each carbon (only effective if sel2 not given and sel1 is not to be passed straight down to MDAnalysis)"
    )


    optParser.add_argument(
        '-thickness', action="store_true",
        dest='thick',
        help="Calculate thickness heatmaps (-leafdiv has to be specified)"
    )



    optParser.add_argument(
        '-thickatom',
        dest='thickatom',
        default="None",
        help="The atomname used for calculating the thickness (if None use same as -divatom) [Default: %(default)s]"
    )





    options = optParser.parse_args()

    if(options.sel2==""):
        options.sel1 = getSelection(options.sel1, options.sepcar)
        print("Using selection string \"%s\" for sel1"%options.sel1)
    else:
        options.sel2 = [options.sel2]

    if(len(options.traj_in) < 4 or (options.traj_in[-4:] not in (".xtc", ".trr"))):
        raise ValueError(
            "File extension not recognised: %s\nOnly .xtc is allowed" % options.traj_in)

    if(len(options.struct_in) < 4 or (options.struct_in[-4:] not in (".tpr", ".pdb", ".gro"))):
        raise ValueError(
            "File extension not recognised: %s\nLegal extensions are .pdb, .gro and .tpr" % options.traj_in)

    if(options.sel2=="" and not options.struct_in.endswith(".tpr")):
        raise ValueError("If the second selection isn't given, a tpr file is needed to get the bonds between atoms")

    if(options.outfile.endswith(".dat")):
        options.outfile = options.outfile[:-4]

    if(options.thick and not options.leaflets):
        raise ValueError("The thickness cannot be calculated without division to leaflets. Either specify -leafdiv or do not specify -thickness")

    if(options.thickatom==None):
        options.thickatom=options.divatom


    return options



def write_time_series(dest, x, y, leaflet="", plot=False, carbnames=[""]):

    lf = leaflet
    if(lf != ""):
        lf = "_"+lf

    for i in range(len(carbnames)):

        carb = carbnames[i]
        if(carb!=""):
            carb = "_"+carb

        fname = "%s%s%s.xvg"%(dest, carb, lf)

        print("Writing timeseries to file %s"%fname)
        header = "Generated by ordermap\nt   order"

        np.savetxt(fname, np.asarray((x[i], y[i])).T, fmt="%.18g", header=header)

        if plot:
            plot_ts(fname)

    if(len(carbnames)>1):
        fname = "%s%s_carbs.xvg"%(dest, lf)

        print("Writing carbon orderparameters to file %s"%fname)

        with open(fname, "w") as f:
            f.write("# Generated by ordermap\n# aname   order\n")
            for i, c in enumerate(carbnames):
                f.write("%s   %.18g\n"%(c, np.mean(y[i])))


        if plot:
            plot_c(fname)









def write_to_file(dest, x, y, data, starttime, endtime, leaflet="", plot=False, carbnames=[""]):
    """
    As the name suggests writes the data to file. If plot=True, the data will also be plotted
    """
    if(endtime!=-1):
        timestr = "_%04dto%04d"%(starttime/1000, endtime/1000)
    else:
        timestr=""

    lf = leaflet
    if(lf != ""):
        lf = "_"+lf

    for i in range(len(carbnames)):
        carb = carbnames[i]
        if(carb!=""):
            carb = "_"+carb
        fname = "%s%s%s%s.dat"%(dest, timestr, carb, lf)

        print("Opening file %s for writing"%fname)
        with open(fname, "w") as f:
            f.write("# Generated by ordermap\n")
            f.write("# averaged over times between %f and %f ps\n"%(starttime, endtime))

            f.write("# x: ")
            for xi in x:
                f.write(" %f"%(xi/10))

            f.write("\n# y: ")
            for yi in y:
                f.write(" %f"%(yi/10))
            f.write("\n")

            for j in range(data.shape[1]):
                for dat in data[i, j]:
                    f.write(str(dat))
                    f.write(" ")
                f.write("\n")

        if(plot):
            name = "Order parameter between %d to %d ns"%(starttime/1000, endtime/1000)
            if(leaflet != ""):
                name += ", %s leaflet"%leaflet
            if(carbnames[i] != ""):
                name += ", %s"%carbnames[i]
            plotdata(fname, name)

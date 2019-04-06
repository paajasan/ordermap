# -*- coding: UTF-8 -*-
import argparse

from src.create_figs import plot as plotdata

def getSelection(sel):
    parts = sel.split()
    if(len(parts)!=3):
        return sel
    rangeparts = parts[2].split("-")
    if(len(rangeparts)==2):
        try:
            start = int(rangeparts[0])
            end   = int(rangeparts[1])
        except (ValueError):
            return sel

        sel = "resname %s and name "%parts[0]
        for i in range(start, end+1):
            sel += "%s%d "%(parts[1], i)
    elif(len(rangeparts)==1):
        try:
            start = int(rangeparts[0])
        except (ValueError):
            return sel
        sel = "resname %s and name %s%d "%(parts[0], parts[1], start)
    else:
        return sel


    return sel[:-1]



def optP():

    description = "Calculates the order parameter"
    optParser = argparse.ArgumentParser(description=description)

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
        '-center',
        dest='center',
        metavar="selection_str",
        default=None,
        help="A selection string to define a group, which will be used as a center for the grids (located at (0, 0))"
    )

    optParser.add_argument(
        '-dt', type=int,
        dest='dt',
        metavar="N",
        default=2500,
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



    options = optParser.parse_args()

    if(options.sel2==""):
        options.sel1 = getSelection(options.sel1)
        print("Using selection string \"%s\" for sel1"%options.sel1)

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

    return options








def write_to_file(dest, x, y, data, starttime, endtime, leaflet="", plot=False):
    """
    As the name suggests writes the data to file. If plot=True, the data will also be plotted
    """
    fname = "%s%04dto%04d%s.dat"%(dest, starttime/1000, endtime/1000, leaflet)

    print("Opening file %s for writing"%fname)
    with open(fname, "w") as f:
        f.write("# Generated by calculate_order\n")
        f.write("# averaged over times between %f and %f ps\n"%(starttime, endtime))

        f.write("# x: ")
        for xi in x:
            f.write(" %f"%(xi/10))

        f.write("\n# y: ")
        for yi in y:
            f.write(" %f"%(yi/10))
        f.write("\n")

        for i in range(data.shape[0]):
            for dat in data[i]:
                f.write(str(dat))
                f.write(" ")
            f.write("\n")

    if(plot):
        name = "Order parameter between %d to %d ns"%(starttime/1000, endtime/1000)
        if(leaflet != ""):
            name += ", %s leaflet"%leaflet
        plotdata(fname, name)

# Order parameter calculator ordermap
#### v0.2

A program for calculating deuterium order parameter heatmaps. Also the timeseries of the parameter can be made, as well as calculating the average over the whole trajectory.



## Getting started


### Requirements

1. MDAnalysis
2. python3
3. Numpy, (Matplotlib, mpl_toolkits)

Basically install [MDAnalysis](https://www.mdanalysis.org/), and other than that the basic python3 should have everything ready.

**Notice**: I have so far only tested this with `python3.7` (3.7.1, to be exact), and though it should work with other `python3.x`, please keep this in mind. If you do use other versions, please let [me](mailto:santeri.e.paajanen@helsinki.fi) know if it works or not.

### Installing

Since this is in python, you don't of course need to "install" anything (except for the requirements mentioned above). Just copy this project to your computer and you're set.

I would suggest adding the projects root folder to your PATH, if you plan on using the program more than just once. In the next parts I will assume the project is in your path, if not in all commands use `/path/to/folder/ordermap/ordermap` instead of `ordermap` (where `/path/to/folder/ordermap/` is naturally the path to the root folder of this project)



## Basic usage


By running the command
```
ordermap -h
```

you'll get a listing of all the possible arguments.

This could be an example run:
```
ordermap -f trajectories/md_fitted.xtc -s trajectories/start_nopbc.tpr -o order/order -center protein -leafdiv -sel1 "POPC C 29-42" -timeseries -sepcarbs
```


### Options

#### Required arguments

| flag | input | description |
| --- | - | - |
| `-f` | fname | Trajectory [<.xtc/.trr>] |
| `-s` | fname | Topology [<.tpr>] |
| `-n` | fname | Index file * |


\* If `-sel1` is given (see below), index file is not required. If both are given, the index file is used.

##### Index file

If the index file has only as single group, it is assumed to be `-sel1` and the `-sel2` is assumed to not be given. If another one is given this is used as `-sel2`. In case more than two groups are present, they are assumed to be carbons in a chain, in the right order, like with `gmx order`.

In any case, the index file cannot contain any additional groups, like "System" or "Protein".

#### Optional arguments


| flag| input  | description | default |
| --- |:-:| - | :-: |
| `-sel1` | str | Selection string (see syntax below) |
| `-sel2` | str | Selection string [1.1] |
| `-o` | out | Output destination [2] | "order.dat" |
| `-to` | out | Output destination for thickness (basically the same as `-o`) | "thickness.dat" |
| `-b` | int | First frame to read from trajectory | 0 |
| `-e` | int | Last frame to read from trajectory [3] | -1 |
| `-dt` | int | Only use every dt frame | 1 |
| `-center` | str | Selection string [1.2] |
| `-davg` | int | Number of frames to use for each heatmap [4] | -1 |
| `-gridn` | int | Number of gridcells along x and y (NxN grid) | 20 |
| `-mindat` | int | Minimum number of datapoints in each cell (else NaN) | 10 |
| `-plot` | - | If set, makes plots of the data [5] |
| `-leafdiv` | - | If set, the selection is divided to two leaflets [6] |
| `-divatom` | str | The atom name to use for the division to leaflets [6] | "P*" |
| `-timeseries` | - | If set, also plots the timeseries for the selection(s) [7] |
| `-sepcarbs` | - | If set, calculates everything for each carbon in the chain separetely [8] |
| `-thickness` | - | If set, also calculates membrane thickness heatmaps [9] |
| `-thickatom` | str | The atom name to use for calculating the thickness (None means use same as `-divatom`) [9] | None |

1. Selections (both have to be selection strings handled by [`select_atoms`](https://www.mdanalysis.org/docs/documentation_pages/selections.html#selection-keywords) of MDAnalysis)
    1. If `-sel2` is not given, the second selection will be all the hydrogen bonded to the atoms of `-sel1` (and I'd suggest using this always).  
    If for some reason this doesn't work, or you want to use `-sel2` for some other reason, make sure it will be the same size as `-sel1` and when ordered by atom number, the bonded atoms are found in the same index in each list.
    2. If `-center` is given, the center of mass of this selection will be handled as the origon.  
    Could be for example `"protein"` so then the x and y coordinates are simply distances to the protein center of mass along x or y.

2. Output destination with or without the `.dat` extension. Each filename will also include information on the time, leaflet (if `-leafdiv`) and carbon (if `-carbon`). For both `"order"` and `"order.dat"` the final names will be something like `order_0000to0250_C29_upper.dat`.

3. `-b` can also be negative, in which case it will be counted from the end, -1 being the last frame (like with Python list indices, except that the end is inclusive). This way a value of `-1` means use all of the frames.

4. The heatmaps will be saved every `-davg` frames. For example if you have 100 ps frames and `-davg 2500`, your heatmaps will always be averaged over 250 ns and then saved.  
If davg<1, then the whole trajectory will be averaged over.

5. Specifying `-plot` makes it so that the data will also be plotted and saved in `.png`-format.  
The plots are meant to be a rather quick sanity check into the output and I suggest writing your own scripts for actual plotting (though feel free to use the code in `src/create_figs.py`).
To make the plots you will naturally need to have `matplotlib` (and `mpl_toolkits`) installed.

6. With `-leafdiv` set the lipids will be divided to upper and lower leaflet along the membrane normal (for now assumed to be z-axis). The coordinate of `-divatom` will be used to do this.  
The `-divatom` should be a selection string for MDAnalysis that only yields a single atom per lipid residue, for example the default "P*" will use the phosphorus atoms of the lipids.

7. The timeseries will be in files with the same naming convention as the heatmaps' output, but with filename extension `.xvg` (for example `order_C29_upper.xvg`).

8. The possibility to calculate the parameter for each carbon separetely is only available if `-sel2` is not given and `-sel1` follows the syntax specified above.  
If `-timeseries` is also specified, an average per carbon will also be calculated into a file ending with `_carbs.xvg`. *Note:* the first column in this file is the atomname, which is a string, so it can't be read with numpys `loadtxt`-function.



##### Selection string syntax

The selection string should follow a set syntax. At least for now, if it doesn't parse correctly according to said syntax, the string will be passed to MDAnalysis' [`select_atoms`](https://www.mdanalysis.org/docs/documentation_pages/selections.html#selection-keywords) and should then be a valid selection string for that. However, for reliability, I might deprecate that later (or have a separate flag for such selection string) and would suggest simply using the syntax specified here:

The syntax is `rname aname anum[-anum]`  
Here rname is the residue name (for example "POPC"), aname is the "first part" of the atom name (propably C, since we usually just want the carbons) and the anums are the "second part" of the atom name (must be integer) and may specify a range.

**Examples:**  
With `POPC C 29-42` we would choose all carbons between C29 and C42 in all POPC molecules.  
With `POPC C 30` we would only choose the carbon C30 in all POPC molecules.

*Note:*  
This of course won't work if the carbon naming is somehow all weird, or the second part isn't an integer. In this case it's best to use index files.

*Note 2:*  
This syntax is only valid for `-sel1` and not for `-sel2` or `-center`.


## Notes

**Parameter calculation**  
The parameter is the deuterium order parameter defined as  
![3/2<cos2(theta)>-1/2](http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7B3%7D%7B2%7D%5Cleft%5Clangle%5Ccos%5E2%28%5Ctheta%29%5Cright%5Crangle-%5Cfrac%7B1%7D%7B2%7D&bc=Transparent&fc=Black&im=png&fs=12&ff=arev&edit=0)  
where theta is the angle between the membrane normal and the C-H bond. For now it can only be calculated with an all-atomistic simulation and no correction is made for double bonds. Also the membrane normal cannot be specified, and is assumed to be the z-axis.

**Bonded hydrogens**  
The hydrogens' atomnames are assumed to match the string `"H*"` and no other atoms (that are bonded to the atoms of `-sel1`) should match it.

**Thickness**
The thickness is calculated so that the height of the `-thickatom` (along the membrane normal) is averaged into a 2d grid for both leaflets and each cell that has data for both leaflets, will have the difference of these as the thickness. This is repeated over each timestep and the average over time is the final output.

**Arguments**  
Some arguments may not be foolproof, so giving  them a value outside a useful range may yield unexpected results or just crash the program.  
Especially `-sel1` may work unexpectedly when trying to give selection string to MDAnalysis and it happens to comply to the syntax discussed above.

**Saved data**  
The data in for the heatmaps is saved so the y-coordinate runs over rows and x over columns. This way you can read it with numpys `loadtxt` and y will be the 0th axis and x the 1st. Then you can simply give it to matplotlibs `contourf` and it should be oriented the right way.

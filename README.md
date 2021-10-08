# dispersion.py
CASTEP post processing for plotting bandstructures and phonon dispersions

# General Usage


For basic usage, dispersion.py can be called from the commandline to plot bandstructures:
```
dispersion.py <seed>
```
where ```<seed>``` is the CASTEP seed. This will read the ```<seed>.cell``` and ```<seed>.bands``` files to plot a bandstructure. 

If instead a phonon dispersion is required, use the ```--phonon``` flag
```
dispersion.py <seed> --phonon
```
Phonon dispersion and bandstructure figures can both be customised using the same flags, e.g.```--multi``` and ```--fig``` for multicoloured lines and figure captions respectively.

```dispersion.py``` can also be used to format a ```SPECTRAL_KPOINT_PATH``` to be included in your ```<seed>.cell``` file. If you do not know which path to use, one can be automatically generated using the crystal symmetry and Atomic Simulation Environment (ASE): 
```
dispersion.py <seed> --path
```
or if you have a path in mind, you can specify it using the Brillouin zone labels, e.g.:
```
dispersion.py <seed> --path W G X W L G
```

# Advanced Bandstructure Plotting

Here are some other more advanced options for plotting bandstructures. 

Spin Orbit Coupling
-------------------
It is sometimes required that calculations are performed with and without SOC, and it is often useful to compare these bandstructures, using the ```--soc``` flag allows you to overlay another set of bands onto the original. You must have the bands file from a SOC calculation in the same directory, ```<soc-seed>.bands```, and can be plotted by:
```
dispersion.py <seed> --soc <soc-seed>
```
It is important to note that the two calculations must have been performed on the same cell. 

Density of States
-----------------

It is possible to plot a DOS along the righthand side of a bandstructure using the output datafile of an OPTADOS calculation:
```
dispersion.py <seed> --dos <optados file>.dat 
```
Multiple DOS can be plotted by providing multiple data files. As above, it is important the DOS and the bandstructure are from the same unit cell as there is no automatic way of checking.

Projected Density of States
---------------------------

```dispersion.py``` is able to read CASTEP ```<seed>.pdos_bin``` files that contain the projections of the wavefunction onto the atomic orbitals. This can be used to colour the bands in the bandstructure by their orbital character:
```
dispersion.py <seed> --pdos
```

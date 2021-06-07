# ADV_Cor
A Python package for spatially variable advection correction procedures.


# Dependencies 

1. numpy
2. scipy
3. fortran complier

# Description 

ADV_Cor includes spatially variable advection correction procedures that are often used in radar analysis products. These procedures take a gridded scalar field at two input times and finds the pattern translation components and scalar field at analysis times between the two input times. Currently, ADV_Cor includes two-dimensional spatially variable advection correction (Shapiro et al. 2010) and three-dimensional spatially variable advection correction (Gebauer et al. 2021; in prep). There are plans to include radial velocity spatially variable advection correction (Shapiro et al. 2021), but it is not currently included.

# Procdures

ADV_Cor includes the following procedures:

GalChen2D - A function that takes two input two-dimensional scalar fields and finds spatially constant pattern translation components that can be used as a first guess for the the other procedures

GalChen3D - A function similar to GalChen2D, but for three-dimensional scalar fields

ADV2D - The two-dimensional spatially variable advection correction procedure

ADV3D - The three-dimensional spatially variable advection correction procedure


# Installation

Currently, ADV_Cor is not a true python package. Before ADV_Cor will work the relaxation solvers functions will need to be compiled with f2py.

The quickest and easist way to do this is simply:

> python -m numpy.f2py -c relax.f90 -m relax

This will create a .so file and you then should be good to go!
Long term ADV_Cor will be an installable python package and this f2py step will be done on installation. I just have to learn how to do this.....


# How to Use

ADV_Cor should play nicely with PyART! Here is how ADV_Cor could be used:

1. PyART is used to open, edit, and grid radar data

2. (a) The now gridded radar data is passed to either GalChen2D or GalChen3D to obtain arrays with the first guess pattern translation components or (b) The user creates arrays with the first guess pattern translation components of their choosing or (c) this step is ignored and the first guess pattern translation components assumed to be zero (!!NOT RECOMMENDED DUE TO POTENTIAL FOR SOLUTION NONUNIQUENESS!!).

3. The gridded radar field and first guess pattern translation components are passed to ADV2D or ADV3D for the advection correction.

4. Wait. This is a computationally intensive iterative procedure. Depending on the size the gridded radar fields, millions of trajectories may be needed each iteration. If you call the function with the verbose arguement set to true, the procedure will output the maximum change in pattern translation component from each iteration so you know how close the procedure is to completing.

5. Enjoy your advection corrected output!

# References

If you use ADV_Cor please cite the relevant work

Shapiro, A., K. M. Willingham, and C. K. Potvin, 2010: Spatially variable advection correction of radar data. Part I: Theoretical considerations. J. Atmos. Sci.,   67, 3445-3456.

Shapiro, A., K. M. Willingham, and C. K. Potvin, 2010: Spatially variable advection correction of radar data. Part II: Test results. J. Atmos. Sci., 67, 3457-3470.

Shapiro, A., J. G. Gebauer, N. A. Dahl, D. J. Bodine, A. Marhre, and C. K Potvin, 2021: Spatially variable advection correction of Doppler radial velocity data. J. Atmos. Sci., 78, 167-188.

Gebauer, J. G., A. Shapiro, C. K. Potvin, and N. A. Dahl, 2021: Improvements to spatially variable advection correction. In Prep.

# ADV_Cor
A Python package for spatially variable advection correction procedures.


# Dependencies 

1. numpy
2. scipy
3. fortran complier

# Description 

ADV_Cor includes spatially variable advection correction procedures that are often used in radar analysis products. These procedures take a gridded scalar field at two input times and finds the pattern translation components and scalar field at analysis times between the two input times. Currently, ADV_Cor includes two-dimensional spatially variable advection correction (Shapiro et al. 2010) and three-dimensional spatially variable advection correction (Gebauer et al. 2021; in prep). There are plans to include radial velocity spatially variable advection correction (Shapiro et al. 2021), but it is not currently included.

# Procdures

## **GalChen2D**(field1, field2, delta_t, dx, missing=999.)

>A function that takes two input two-dimensional scalar fields and finds spatially constant pattern translation components that can be used as a first guess for   the the other procedures

  Inputs:
  
    field1 - (y,x) array for the first input time
    field2 - (y,x) array for the second input time
    delta_t - Float for the time difference between input fields
    missing - Value for missing data if NaNs are not used in input fields

  Outputs:
  
    U - (y,x) array for spatially constant x-component of pattern translation
    V - (y,x) array for spatially constant y-component of pattern translation

## **GalChen3D**(field1, field2, delta_t, dx, missing=999.)

> A function similar to GalChen2D, but for three-dimensional scalar fields
  
  Inputs:
  
    field1 - (z,y,x) array for the first input time
    field2 - (z,y,x) array for the second input time
    delta_t - Float for the time difference between input fields
    missing - Value for missing data if NaNs are not used in input fields
    
  Outputs:
  
    U - (z,y,x) array for spatially constant x-component of pattern translation
    V - (z,y,x) array for spatially constant y-component of pattern translation
    W - (z,y,x) array for spatially constant z-component of pattern translation

## **ADV2D**(field1, field2, first_U, first_V, dx, dy, bigT, nt ,dt, beta, relax=1, under=1, itermax=20000, intermainmax=100, tol=0.001, tol2=0.01, missing=999, verbose=True)

>The two-dimensional spatially variable advection correction procedure

  Inputs:
    
    field1 - (y,x) array for the first input time
    field2 - (y,x) array for the second input time
    first_U - (y,x) array with the first guess U pattern translation component
    first_V - (y,x) array with the first guess V pattern translation component
    dx - float for the x grid spacing
    dy - float for the y grid spacing
    bigT - float for time difference between two input files
    nt - int for number of timesteps to advection correct data two, including the two input times
    dt - float for the length of the timesteps
    beta - float for smoothing parameter
    relax - float for over-relaxation factor between 1 and 2 for the relaxation solution of the pde's
    under - float for under-relaxation factor beween 0 and 1 for the outer loop solution
    intermax - int for maximum number of iterations for relaxation solution
    intermaxinmax - int for maximum number of iterations for outer loop of procedure
    tol - float for tolerance of relaxation solution
    tol2 - float for tolerance of outer loop of procedure
    missing - Value for missing data in NaNs are not used in the input fields
   
  Outputs:
    
    U - (y,x) array of the x-component of pattern translation
    V - (y,x) array of the y-component of pattern translation
    ref - (y,x,nt) array of the advection corrected scalar field


## **ADV3D**(field1, field2, first_U, first_V, first_W, dx, dy, dz, bigT, nt, dt, beta, gamma, eta, nu, relax=1, under=1, itermax=20000, itermainmax=100, tol=0.001, tol2=0.01,missing=999.,verbose=999.)

>The three-dimensional spatially variable advection correction procedure

  Inputs:
  
    field1 - (z,y,x) array for the first input time
    field2 - (z,y,x) array for the second input time
    first_U - (z,y,x) array with the first guess U pattern translation component
    first_V - (z,y,x) array with the first guess V pattern translation component
    first_W - (z,y,x) array with the first guess W pattern translation component
    dx - float for the x grid spacing
    dy - float for the y grid spacing
    dz - float for the z grid spacing
    bigT - float for time difference between two input files
    nt - int for number of timesteps to advection correct data two, including the two input times
    dt - float for the length of the timesteps
    beta - float for horizontal smoothing parameter of the horizontal pattern translation components
    gamma - float for the vertical smoothing parameter of the horizontal pattern translation components
    eta - float for the horizontal smoothing parameter of the vertical pattern translation components
    nu - float for the vertical smoothing parameter of the vertical pattern translation components
    relax - float for over-relaxation factor between 1 and 2 for the relaxation solution of the pde's
    under - float for under-relaxation factor beween 0 and 1 for the outer loop solution
    intermax - int for maximum number of iterations for relaxation solution
    intermaxinmax - int for maximum number of iterations for outer loop of procedure
    tol - float for tolerance of relaxation solution
    tol2 - float for tolerance of outer loop of procedure
    missing - Value for missing data in NaNs are not used in the input fields
    
  Outputs:
    
    U - (z,y,x) array for the x-component of the pattern translation component
    V - (z,y,x) array for the y-component of the pattern translation component
    W - (z,y,x) array for the z-component of the pattern translation component
    ref - (z,y,x,t) array of the advection corrected field
    
## **precomputed_ADV2D**(field1, field2, u, v, dx, dy, dt, nt)

> Advect a two-dimensional field from computed pattern translation components

  Inputs:
      
    field1 - (y,x) array for the first input time
    field2 - (y,x) array for the second input time
    u - (y,x) array for the x-component of the pattern translation
    v - (y,x) array for the y-component of the pattern translation
    dx - float for the x grid spacing
    dy - float for the y grid spacing
    dt - float for the length of the timesteps
    nt - int for number of timesteps to advection correct data two, including the two input times
  
  Outputs:
    
    ref - (y,x,nt) array of the advection corrected field

## **precomputed_ADV3D**(field1, field2, u, v, w, dx, dy, dz, dt, nt)

> Advect a two-dimensional field from computed pattern translation components

  Inputs:
      
    field1 - (z,y,x) array for the first input time
    field2 - (z,y,x) array for the second input time
    u - (z,y,x) array for the x-component of the pattern translation
    v - (z,y,x) array for the y-component of the pattern translation
    w - (z,y,x) array for the z-component of the pattern translation
    dx - float for the x grid spacing
    dy - float for the y grid spacing
    dz - float for the z grid spacing
    dt - float for the length of the timesteps
    nt - int for number of timesteps to advection correct data two, including the two input times
  
  Outputs:
    
    ref - (z,y,x,nt) array of the advection corrected field

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

5. If an additional field needs to be advection corrected using the retrieved pattern translation components this can be done with precomputed_ADV2D or precomputed_ADV3D.

6. Enjoy your advection corrected output!

# References

If you use ADV_Cor please cite the relevant work

Shapiro, A., K. M. Willingham, and C. K. Potvin, 2010: Spatially variable advection correction of radar data. Part I: Theoretical considerations. J. Atmos. Sci.,   67, 3445-3456.

Shapiro, A., K. M. Willingham, and C. K. Potvin, 2010: Spatially variable advection correction of radar data. Part II: Test results. J. Atmos. Sci., 67, 3457-3470.

Shapiro, A., J. G. Gebauer, N. A. Dahl, D. J. Bodine, A. Marhre, and C. K Potvin, 2021: Spatially variable advection correction of Doppler radial velocity data. J. Atmos. Sci., 78, 167-188.

Gebauer, J. G., A. Shapiro, C. K. Potvin, and N. A. Dahl, 2021: Improvements to spatially variable advection correction. In Prep.

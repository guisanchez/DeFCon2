2D-Debris
=========

Debris flow simulation code. Current status: ... in progress...

Shallow water equations are taken to solve numerically the dynamic of debris flow and avalanches. The numerical method is upwind and first order.


Compile 2D-Debris
-----------------

Requisites for linux users: Make and gfortran compiler. 

* Get into the 2D-Debris folder. Run makefile (Type make in terminal).
* Check that an executable file called main2D has been created in your directory. 

Windows users should wait. 

Input files
-----------

There are 3 required input files: a mesh file, an initial conditions file and a code inputs/options file.

1.- Mesh

The code reads [ASCII grid](http://en.wikipedia.org/wiki/Esri_grid) files. The format is:

    nvert      	     4
    ncols      	     120
    nrows      	     130
    xllcorner  	     12.0
    yllcorner 	     -48.6
    cellsize	     4.0
    NODATA_value     -9999
    Z(1,1:ncols) 
    ...
    Z(nrows,1:ncols)

Where Z(i,j) is the bed height at row i and column j

2.- Initial

Init file is [ASCII grid](http://en.wikipedia.org/wiki/Esri_grid) as well. 

    nvert      	     4
    ncols      	     120
    nrows      	     130
    xllcorner  	     12.0
    yllcorner 	     -48.6
    cellsize	     4.0
    NODATA_value     -9999
    h(1,1:ncols) 
    ...
    h(nrows,1:ncols)

Where h(i,j) is the debris depth at row i and column j

You may find some examples in the [2d-debris project page ](http://guisanchez.github.com/2D-Debris/).

3.- Inputs/options file.

The code reads the [JSON](http://www.json.org/) file `input.js`. Special thanks to [@josephalevin](https://github.com/josephalevin/fson), who developed the fortran library we are working with. Several data are read from `input.js`:

    {
      "time" :
        {
          "cfl"      : 0.9,
          "end"      : 300.0,
          "vtklag"   : 5.0,
          "stdprint" : 50
        },
      
      "params" :
        {
          "tanphi" : 0.226295,
          "Xi"     : 0.004888,
          "k"      : 1.0
        },
        
      "preproc_mesh"   : 1,
      "preproc_err"    : 1,
      "standard_out"   : 1,
      "generate_vtk"   : 1,
      "compare_result" : 0,
      
      "mesh"           : "mesh.txt",
      "h_initial"      : "h_initial.txt"
    }
    	
Further details on this file in soon available in the [2d-debris project page ](http://guisanchez.github.com/2D-Debris/).

    
Run 2D_Debis
------------

Run the file generated after make. Type ./main_2D

Watch your result!
------------------

Requisites: [paraview](http://www.paraview.org)

If you generated vtk files in your run, you may watch them with paraview
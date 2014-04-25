DeFCon2
=======

Debris flow simulation code. Current status: ... in progress...

Shallow water equations are taken to solve numerically the dynamic of debris flow and avalanches. The numerical method is upwind and first order.


Compile DeFCon2
---------------

Requisites for linux users: Make and gfortran compiler. 

* Get into the DeFCon2 folder. Run makefile (Type make in terminal).
* Check that an executable file called main2D has been created in your directory. 

Windows users should wait. 

Input files
-----------

There are 3 required input files: a mesh file, an initial conditions file and a code inputs/options file.

1.- Mesh

The code reads ASCII mesh files. The format is:

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

Init file is ASCII as well. 

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

You may find some examples in the [DeFCon2 project page ](http://guisanchez.github.com/2D-Debris/).

3.- Inputs/options file.

The code reads a json format file called input.json. Special thanks to [@josephalevin](https://github.com/josephalevin/fson), who developed the fortran library we are working with. Several data are read from input.json:

    {
    	"time" : 
    	{
		"cfl": 0.9,
		"end": 300.0,
		"vtklag": 5.0,
		"stdprint": 50
    	},
    	"params" : 
    	{
		"tanphi": 0.226295,
		"Xi": 0.004888,
		"k": 1.0
    	},
    	"preproc_mesh": 1,
    	"preproc_err": 1,
    	"standard_out": 1,
    	"compare_result": 0,
    	"generate_vtk": 1,
    	"mesh": "malla0.txt",
    	"h_initial": "init0.txt",
		     "stop" :
		     {
		     "stime": 100.0,
		     "Ekstop": 10000.0
		     }
	}
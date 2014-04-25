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
    	"time": 
    	{
		"cfl": 0.9,
		"end": 300.0,
		"vtklag": 5.0,
		"stdprint": 50
    	},
    	"friclaw": "voellmy",
    	"params": 
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
	"stop":
    	{
		"stime" : 100.0,
		"Ekstop" : 10000.0
    	}

    }

In "time" you can set the cfl number, "end", the simulation lenght of time (in seconds), the time lapse (in seconds) and the period for standard output (in cycles)
In "friclaw" you can choose the rheology law appropriated for your simulation. Please type "voellmy", "coulomb" or "bingham".
Next, you need to introduce the rheology parameters: "params". These include the earth pressure coefficient ("k"), and, in voellmy rheology, the equilibrium slope ("tanphi") and the turbulent coefficient ("Xi"). If you work with Bingham or Coulomb rheologies, you may need to input the yield stress ("yieldstress") and the consistency index ("mu"). Please note that Coulomb rheology needs also an apparent equilibrium slope.
Next, there are several binary orders (0 / 1):
If you are expecting to run the code several times for the same mesh, I suggest to you to set "preproc_mesh": 1. Then, the code will create a binary file so it can access to mesh data much faster, in the following runs. If, in addition you are evaluating the simulation error, set "preproc_err" to 1 as well. 
If you like to see the standard output set "standard_out" to 1. The code types the step number, time, last time step, "Kinetic Energy" and mass error.
If you have some expected deposit for your debris flow, you can evaluate the error of your simulation. Set "compare_result": 1 and the code will provide you a Nash-Shutcliffe statistic. 
If you set "generate_vtk" to 1, the code will write vtk files so that you can watch your simulation in paraview.
Next, three input mesh files: The bed ("mesh"), the initial depth ("h_initial") and observed/expected deposit ("h_observ"). The former is optional (needed for error estimation) and its format is the same as first two.

Sometimes, it is interseting to stop the simulation, specially if the sliding mass is almost halt, i.e. if its kinetic energy is lower than some threshold you may introduce. This order will be affective only for time > stime.


Compile and run
---------------

After downloading the code, go to the DeFCon2 folder and type: make
Several warnings may appear, but an executable file is created.
Make sure you have the mesh files and the input.js file.
To run the code, type : ./main2D input.js

Watch your result
-----------------

vtk files can be open by paraview (and others). Enjoy!
2D_Debris
=========

Debris flow simulation code. Current status: ... in progress...

Shallow water equations are taken to solve numerically the dynamic of debris flow and avalanches. The numerical method is upwind and first order.


Compile 2D_Debris
_________________

Requisites for linux users: Make and gfortran compiler. 

1.- Get into the 2D-Debris folder. Run makefile (Type make in terminal).
2.- Check that an executable file called main2D has been created in your directory. 

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

3.- Inputs/options file.
The code reads a json format file called input.json. Special thanks to @josephalevin, who developed the fortran library we are using. Several data are read from input.json:
    time data:
    	 cfl (real number)  
	     Courant number for the time step
	 end (real number)
	     Simulation (max) lenght of time
	 vtklag (real number)
	     Simulated time between vtk outputs
	 stdprint (integer)
	     Number of time steps between standard out
    parameters data:
	 tanphi (real number)
	 	Friction expression: equilibrium slope
	 Xi (real number)
                Friction expression: dynamic component
	 k (real number)
                Debris flow dynamics: Earth pressure term.
    preproc_mesh (integer)
	 1 if mesh has already been processed and there exist a binary mesh file
    preproc_err (integer)
    	 1 if expected values file has been processed and written to a binary file
    standard_out (integer)
         1 if you want the code to print standard outputs
    compare_result (integer)
         1 if you have (expected) data to compare with your final output
    generate_vtk (integer) 
         1 if you want the code to write vtk files (eg for data visualization)
    mesh (string)
    	 Name of the mesh file
    h_initial (string)
    	 Name of the inital h file

Run 2D_Debis
------------

Run the file generated after make. Type ./main_2D

Watch your result!
------------------

Requisites: paraview

If you generated vtk files in your run, you may watch them with paraview
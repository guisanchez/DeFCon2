{
    "time" : 
    {
	"cfl": 0.9,
	"end": 100.0,
	"vtklag": 2.0,
	"stdprint": 50
    },
    "friclaw": "voellmy",
    "params" : 
    {
	"rho" : 1.0,
	"yieldstress" : 0.1,
	"mu" : 35.0,
	"tanphi": 0.15,
	"Xi": 0.01,
	"k": 1.0
    },
    "preproc_mesh": 1,
    "preproc_err": 1,
    "standard_out": 1,
    "compare_result": 1,
    "generate_vtk": 1,
    "mesh": "meshfile.asc",
    "h_initial": "hfile.asc",
    "h_observ": "deposit.asc",
    "stop" :
    {
	"stime" : 100.0,
	"Ekstop" : 10000.0
    }
}
	
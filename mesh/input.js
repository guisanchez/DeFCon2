{
    "time" : 
    {
	"cfl": 0.9,
	"end": 30.0,
	"vtklag": 1.0,
	"stdprint": 10
    },
    "friclaw": "Voellmy",
    "params" : 
    {
	"tanphi": 0.1,
	"Xi": 0.02,
	"k": 1.0
    },
    "preproc_mesh": 0,
    "preproc_err": 1,
    "standard_out": 1,
    "compare_result": 0,
    "generate_vtk": 1,
    "mesh": "examplemesh1.txt",
    "h_initial": "exampleinitial1.txt",
    "stop" :
    {
	"stime" : 100.0,
	"Ekstop" : 10000.0
    }
}

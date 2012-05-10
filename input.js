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
    
  "mesh"           : "mesh.txt",
  "h_initial"      : "h_initial.txt",
  "preproc_mesh"   : 1,
  "preproc_err"    : 1,
  "standard_out"   : 1,
  "generate_vtk"   : 1,
  "compare_result" : 0
}
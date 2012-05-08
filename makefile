src = src/
fson = $(src)fson/
cflags = -c -O3 -Wall -I$(fson)
cc = gfortran
obj_lib = $(fson)fson_string_m.o $(fson)fson_value_m.o $(fson)fson_path_m.o \
	$(fson)fson.o
obj_main = $(src)main2D.o $(src)calc_paredes2D_Sf_fix3.o $(src)method2Drd.o \
	$(src)matprod.o $(src)compute_Sf2D.o $(src)read_regmesh2DP.o \
	$(src)rutina_comparadorP.o $(src)out_vtk.o

main2D: $(fson)libfson.a $(obj_main) makefile
	$(cc) -o main2D $(obj_main) $(fson)libfson.a

$(fson)fson_string_m.o: $(fson)fson_string_m.f95 makefile
	$(cc) $(cflags) $(fson)fson_string_m.f95 -o $(fson)fson_string_m.o

$(fson)fson_value_m.o: $(fson)fson_value_m.f95 makefile
	$(cc) $(cflags) $(fson)fson_value_m.f95 -o $(fson)fson_value_m.o

$(fson)fson_path_m.o: $(fson)fson_path_m.f95 makefile
	$(cc) $(cflags) $(fson)fson_path_m.f95 -o $(fson)fson_path_m.o

$(fson)fson.o: $(fson)fson.f95 makefile
	$(cc) $(cflags) $(fson)fson.f95 -o $(fson)fson.o

$(fson)libfson.a: $(obj_lib) makefile
	ar -rv $(fson)libfson.a $(obj_lib)

$(src)main2D.o: $(src)main2D.f makefile
	$(cc) $(cflags) $(src)main2D.f -o $(src)main2D.o

$(src)calc_paredes2D_Sf_fix3.o: $(src)calc_paredes2D_Sf_fix3.f makefile
	$(cc) $(cflags) $(src)calc_paredes2D_Sf_fix3.f \
		-o $(src)calc_paredes2D_Sf_fix3.o

$(src)method2Drd.o: $(src)method2Drd.f makefile
	$(cc) $(cflags) $(src)method2Drd.f -o $(src)method2Drd.o

$(src)matprod.o: $(src)matprod.f makefile
	$(cc) $(cflags) $(src)matprod.f -o $(src)matprod.o

$(src)compute_Sf2D.o: $(src)compute_Sf2D.f makefile
	$(cc) $(cflags) $(src)compute_Sf2D.f -o $(src)compute_Sf2D.o

$(src)read_regmesh2DP.o: $(src)read_regmesh2DP.f makefile
	$(cc) $(cflags) $(src)read_regmesh2DP.f -o $(src)read_regmesh2DP.o

$(src)rutina_comparadorP.o: $(src)rutina_comparadorP.f makefile
	$(cc) $(cflags) $(src)rutina_comparadorP.f \
		-o $(src)rutina_comparadorP.o

$(src)out_vtk.o: $(src)out_vtk.f makefile
	$(cc) $(cflags) $(src)out_vtk.f -o $(src)out_vtk.o

$(src)latex/refman.pdf: main2D $(src)Doxyfile makefile
	cd $(src); doxygen; cd latex; make;

clean:
	rm $(obj_main) $(obj_lib)

build:
	make clean
	make

doc:
	make $(src)latex/refman.pdf

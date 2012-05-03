#
#
# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -O3 -Wall
PROGRAM = main2D
main2D: read_regmesh2DP.o calc_paredes2D_Sf_fix3.o method2Drd.o out_vtk.o rutina_comparadorP.o matprod.o compute_Sf2D.o
all: $(PROGRAM)
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
%.o: %.f
	$(FC) $(FCFLAGS) -c $<


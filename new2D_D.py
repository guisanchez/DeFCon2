import os
n1 = 'j'
noes = 0
yes = ['Y','y']
no = ['N','n']
ans = [ yes[0], yes[1], no[0], no[1] ]
while (n1 not in ans):
    n1 = raw_input("Need standard out? (y/n): ")

if (n1 in no):
   f1 = open('madre.f','r')
   f2 = open('temp0','w')
   key = 'write(*,*)'
   for line in f1:
       if (key not in line):
           f2.write(line)
       else:
           line2 = '!' + line
           f2.write(line2)
   f2.close()
   f1.close()
   noes = noes + 1
else:
    os.system('cp madre.f temp0')
try:
    f1.close()
    f2.close()
except:
    pass    


n1 = 'j'
while (n1 not in ans):
    n1 = raw_input("Need to compare your output with observed result? (y/n): ")

if (n1 in no):
    f1 = open('temp0','r')
    f2 = open('temp1','w')
    key = 'call comparador'
    for line in f1:
        if (key not in line):
            f2.write(line)
        else:
            line2 = '!' + line
            f2.write(line2)
    f2.close()
    f1.close()
    noes = noes + 1
else:
    os.system('cp temp0 temp1')

n1 = 'j'
while (n1 not in ans):
    n1 = raw_input("Need vtk (paraview) output files? (y/n): ")

if (n1 in no):
    f1 = open('temp1','r')
    f2 = open('main2D.f','w')
    key = 'call write_plt'
    for line in f1:
        if (key not in line):
            f2.write(line)
        else:
            line2 = '!' + line
            f2.write(line2)
    f2.close()
    f1.close()
    noes = noes + 1
else:
    os.system('cp temp1 main2D.f')

if (noes == 3): 
    print 'Does this make sense?\n'
os.system('rm temp1 temp0')
print 'Created main2D.f file\n'
print 'Please compile : gfortran main2D.f calc_paredes2D_Sf_fix3.f matprod.f read_regmesh2DP.f compute_Sf2D.f method2Drd.f out_vtk.f rutina_comparadorP.f -O3 -Wall -o main2D'


    


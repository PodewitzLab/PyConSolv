%Chk=LIG_small_opt.chk
%Mem=3000MB
%NProcShared=2
# B3LYP/6-31G* Geom=PrintInputOrient Integral=(Grid=UltraFine) Opt
 
CLR
 
0  1
Pt         -2.880    0.137    0.180
N          -0.903    0.643    0.527
H          -0.793    1.559    0.991
H          -0.336    0.670   -0.335
H          -0.479   -0.067    1.149
N          -3.179    1.800   -1.017
H          -3.030    2.690   -0.515
H          -4.156    1.800   -1.356
H          -2.567    1.815   -1.848
Cl         -2.495   -1.736    1.546
Cl         -5.129   -0.399   -0.247
 
 

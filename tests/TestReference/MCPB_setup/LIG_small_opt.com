%Chk=LIG_small_opt.chk
%Mem=3000MB
%NProcShared=2
# B3LYP/6-31G* Geom=PrintInputOrient Integral=(Grid=UltraFine) Opt
 
CLR
 
0  1
Pt         -2.872    0.147    0.175
N          -0.890    0.636    0.536
H          -0.770    1.545    1.010
H          -0.312    0.658   -0.318
H          -0.492   -0.090    1.158
N          -3.192    1.807   -1.026
H          -3.053    2.702   -0.532
H          -4.175    1.781   -1.348
H          -2.593    1.825   -1.866
Cl         -2.481   -1.719    1.540
Cl         -5.116   -0.382   -0.253
 
 

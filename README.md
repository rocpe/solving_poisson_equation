# Program which calculates electric potential from density of charge (Poisson eqn) in C++
It uses armadillo library. Program calculates potential u(x,y) from given equation
of density rho(x,y). It also calculates rho_prime(x,y) which is density of electric charge
calculated from calculated previously u(x,y) and delta(x,y) = rho_prime(x, y) - rho(x,y)
# Compilation
Just type "make".
# Files
Program creates 3 files: u.csv, rho_prime.csv, delta.csv which are matrices. You can use your
favorite plotting program to plot them (in R it is just "filled.contour(as.matrix(u.csv)) ! ).

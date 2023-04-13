/*
 * Program which calculates electric potential when density of electric charge is
 * given (Poisson eqn). It uses over-relaxation method. Results of calculation
 * are saved to the text files. You can (and should) tweak in this file based on
 * your needs.
 *
 */
#include <iostream>
#include <armadillo>

#define N 31
#define REPS 500
#define w 1.9

using namespace arma;

/* USED FUNCTIONS */
double neighborValue(const arma::mat& u, const int& i, const int& j);//sum of neighbor next to  values of given point
double rho(const double& xgrid, const double& ygrid, const double& x0 , const double& d );//density of charge
double rho(const double& xgrid, const double& ygrid);//same as above(overload)

/* MAIN */
int main() {

    /* grid in x and y direction is the same */
    const double dx = 1;
    const vec grid = regspace<vec>(-N,dx,N);//vector which holds xes for rho function

    /* vector which holds values of next iterations */
    vec a_vals(REPS, fill::zeros);

    /* matrix 63x63 filled with 0 */
    mat u(2*N+1, 2*N+1, fill::zeros);

    /* initial condition */
    u.submat(N-15, N-15, N+15, N+15) = mat(N,N, fill::ones);


    /* calculating potential u(x,y) */
    for(int rep = 0; rep < REPS; rep++)
    {
        for(long long unsigned i = 1; i < u.n_rows-2; i++)
            for(long long unsigned j = 1; j < u.n_cols-2; j++)
                /* over-relaxation method, check out eqn (6) in https://home.agh.edu.pl/~bszafran/mofit/poi.pdf for more info */
                u(i,j) = (1-w)*u(i,j) + w*(neighborValue(u, i, j) + rho(grid(i), -grid(j)))*(pow(dx,2))/4;
    }

    /* saving matrix to the file */
    u.save("u.csv", csv_ascii);


    /* rho_prime is calculated rho of calculated previously u(x,y) */
    mat rho_prime(2*N+1, 2*N+1, fill::zeros);
    for(long long unsigned i = 1; i < u.n_rows-2; i++)
        for(long long unsigned j = 1; j < u.n_cols-2; j++)
            rho_prime(i,j) = -(neighborValue(u,i,j) - 4*u(i,j))/pow(dx,2);
                                                                         
    /* saving matrix to the file */
    rho_prime.save("rho_prime.csv", csv_ascii);


    
    /* delta is difference between given rho and rho_prime */
    mat delta(2*N+1, 2*N+1, fill::zeros);
    for(long long unsigned i = 1; i < u.n_rows-2; i++)
        for(long long unsigned j = 1; j < u.n_cols-2; j++)
            delta(i,j) = rho_prime(i,j) - rho(grid(i), -grid(j));

    /* saving matrix to the file */
    delta.save("delta.csv", csv_ascii);
}

/* FUNCTIONS */
double neighborValue(const arma::mat& u, const int& i, const int& j)
{
    return( u(i+1,j) + u(i-1,j) + u(i,j-1) + u(i,j+1) );
}

double rho(const double& xgrid, const double& ygrid, const double& x0 , const double& d )
{   /* density of  electric charge, you can change it */
    return( exp( -(pow(xgrid-x0,2) + pow(ygrid,2) )/pow(d,2) ) -
            exp( -(pow(xgrid+x0,2) + pow(ygrid,2) )/pow(d,2) ) );
}
double rho(const double& xgrid, const double& ygrid)
{
    const double x0 = 4; const double d = 4;//default values, you can change them!
    return rho(xgrid, ygrid, x0, d) ;
}

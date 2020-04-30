#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>



double T(double a, double b, int n, double (*f)(double));
double x_square(double x);

const unsigned double Pi = 3.14159265359;
const unsigned double Epsilon0 = 8.8541878128e-12;

int main(int argc, char **argv){

    printf("%f\n", T(1., 10., 5, sin));

    return 0;
}

double x_square(double x){
    double res = x*x;
    return res;
}

double T(double a, double b, int n, double (*f)(double)){
    double res = 1.;
    double sum = 0;
    // double h = (b-a)/n;


    // Speicher allocalisieren ? zumindest so ähnlich
    double* h = (double*)malloc(sizeof(double)* (n+1) );

    double* lines = (double*)malloc(sizeof(double) * 1000000 );
    assert(lines != NULL);
    double **T = (double**)malloc(sizeof(double*) * 1000000);
    assert(T != NULL);
    for (int i = 0; i <= n; i++){
        T[i] = lines + n*i;
    }

    h[0] = (b-a);
    for (int i = 1; i <= n; i++){
        h[i] = h[i-1]/2;
        // printf("%f\n", h[i]);
    }

    // Trapezsumme zu den Schrittweiten h0 - hn :'( (T_i,0)
    for (int i = 0; i <= n; i++){
        sum = 0;
        for (int j = 1; j < (b-a)/h[i]; j++){
            sum += f(a+j*h[i]);
        }
        T[i][0] = h[i]/2*(f(a)+2*sum+f(b));
        // printf("%f\n", T[i][0]);
    }

    // printf("dummer Test: %f\n", T[n][n]);

    // Berechnen der T_i,k-Werte
    for (int i = 1; i <= n; i++){
        for (int k = 1; k <= i; k++){
            T[i][k] = T[i][k-1] * (-pow(h[i-k], 2) / (pow(h[i], 2) - pow(h[i-k], 2) ) ) - T[i-1][k-1] * (-pow(h[i], 2) / ( pow(h[i], 2) - pow(h[i-k], 2) ) );
        }
    }
    
    res = T[n][n];
    return res;
}

double T_inf(double epsilon, int n, double (*f)(double)){
    double res = 1.;
    double sum = 0;
    // double h = (b-a)/n;


    // Speicher allocalisieren ? zumindest so ähnlich
    double* h = (double*)malloc(sizeof(double)* (n+1) );

    double* lines = (double*)malloc(sizeof(double) * 1000000 );
    assert(lines != NULL);
    double **T = (double**)malloc(sizeof(double*) * 1000000);
    assert(T != NULL);
    for (int i = 0; i <= n; i++){
        T[i] = lines + n*i;
    }

    h[0] = (b-a);
    for (int i = 1; i <= n; i++){
        h[i] = h[i-1]/2;
        // printf("%f\n", h[i]);
    }

    // Trapezsumme zu den Schrittweiten h0 - hn :'( (T_i,0)
    for (int i = 0; i <= n; i++){
        sum = 0;
        for (int j = 1; j < (b-a)/h[i]; j++){
            sum += f(a+j*h[i]);
        }
        T[i][0] = h[i]/2*(f(a)+2*sum+f(b));
        // printf("%f\n", T[i][0]);
    }

    // printf("dummer Test: %f\n", T[n][n]);

    // Berechnen der T_i,k-Werte
    for (int i = 1; i <= n; i++){
        for (int k = 1; k <= i; k++){
            T[i][k] = T[i][k-1] * (-pow(h[i-k], 2) / (pow(h[i], 2) - pow(h[i-k], 2) ) ) - T[i-1][k-1] * (-pow(h[i], 2) / ( pow(h[i], 2) - pow(h[i-k], 2) ) );
        }
    }
    
    res = T[n][n];
    return res;
}

double V_integrate(int x){
    res = exp(-pow(x,2)/pow(a,2)/sqrt(pow(x,2)+pow(z,2)));
    return res;
}
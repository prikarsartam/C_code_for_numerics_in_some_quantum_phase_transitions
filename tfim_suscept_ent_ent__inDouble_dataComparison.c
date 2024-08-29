#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

// everything is in Double


double finite_size_sigma_exp(double h, int n) {
    const double PI = 2.0 * acos(0.0);
    double sum = 0.0;

    for (int m_init = - n/2; m_init <= n/2 - 1 ; m_init++) {
        double m = (double)m_init + 0.5;
        double k = 2.0 * PI * (double)m / (double)n;
        double term = (h - cos(k)) / sqrt(pow(h - cos(k), 2) + pow(sin(k), 2));
        sum += term;
    }
    return sum / (double)n;
}


double Derivative_of_finite_size_sigma_exp(double h, int n) {
    const double PI = 2.0 * acos(0.0);
    double sum_n = 0.0;
    for (int m_init = - n/2; m_init <= n/2 - 1 ; m_init++) {
        double m = (double)m_init + 0.5;
        double kval = 2.0 * PI * m / (double)n;
        double term_n = (sin(kval) * sin(kval)) / pow((1.0 + h * h - 2.0 * h * cos(kval)), 1.5);
        sum_n += term_n;
    }
    return sum_n / n;
}




/*
double Derivative_of_finite_size_sigma_exp(double h, int n) {
 double delta = 1e-8;
 double delta_mz = (finite_size_sigma_exp(h + delta, n) - finite_size_sigma_exp(h, n));
 return delta_mz / delta;
 
}
*/



double FiniteSizeSusceptibility(double h, int n) {
    double mz = finite_size_sigma_exp(h, n);
    double dhmz = Derivative_of_finite_size_sigma_exp(h, n);
    
    double numerator = dhmz * dhmz;
    double denominator = 2.0 * (1.0 - mz * mz);
    
    return numerator / denominator;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <samplesize> <n>\n", argv[0]);
        return 1;
    }
    

    int n = atoi(argv[1]);


    char filename[100];
    snprintf(filename, sizeof(filename), "data/tfim_susceptibility_dataComparison_inDoub_%d.csv", n);

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Unable to create file %s\n", filename);
        return 1;
    }

    fprintf(file, "h,susceptibility\n");


    
    double h_vals_array[] = {0.990000000000000000, 0.991000000000000000, 0.992500000000000000,0.995000000000000000,0.997500000000000000,0.999000000000000000, 1.000000000000000000, 1.000500000000000000, 1.001000000000000000, 1.002500000000000000, 1.005000000000000000, 1.007500000000000000, 1.009000000000000000, 1.010000000000000};


const int h_vals_size = sizeof(h_vals_array) / sizeof(h_vals_array[0]);


       printf("\n\n10 Susceptibility values for equally spaced h values:\n\n");
     printf("------------------------------------------------------------\n\n");

    for (int i = 0; i < h_vals_size; i++) {
        double h_val = h_vals_array[i];
        double suscept_val = FiniteSizeSusceptibility(h_val, n);
        fprintf(file, "%2.17f,%2.17f\n", h_val, suscept_val);
        
        printf("h = %2.17f, suscept = %2.17f\n", h_val, suscept_val);
    }

    printf("\n------------------------------------------------------------\n");

    fclose(file);
    printf("\nData written to %s\n\n\n", filename);

    return 0;
}

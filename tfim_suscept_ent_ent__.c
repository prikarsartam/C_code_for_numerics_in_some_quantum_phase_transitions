#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>



double finite_size_sigma_exp(double h, int n) {
    const double PI = 2.0 * acos(0.0);
    double sum = 0.0;

    for (int m_init = - n/2; m_init <= n/2 - 1 ; m_init++) {
        double m = m_init + 0.5;
        double k = 2.0 * PI * (double)m / (double)n;
        double term = (h - cos(k)) / sqrt(pow(h - cos(k), 2) + pow(sin(k), 2));
        sum += term;
    }
    return sum / (double)n;
}

/*
double Derivative_of_finite_size_sigma_exp(double h, int n) {
    const double PI = 2.0 * acos(0.0);
    double sum_n = 0.0;
    for (int k = -(n-1)/2; k <= (n-1)/2; k++) {
        double kval = 2.0 * PI * k / n;
        double term_n = (sin(kval) * sin(kval)) / pow((1.0 + h * h - 2.0 * h * cos(kval)), 1.5);
        sum_n += term_n;
    }
    return sum_n / n;
}
*/


double Derivative_of_finite_size_sigma_exp(double h, int n) {
 double delta = 1e-8;
 double delta_mz = (finite_size_sigma_exp(h + delta, n) - finite_size_sigma_exp(h, n));
 return delta_mz / delta;
 
}


double FiniteSizeSusceptibility(double h, int n) {
    double mz = finite_size_sigma_exp(h, n);
    double dhmz = Derivative_of_finite_size_sigma_exp(h, n);
    
    double numerator = dhmz * dhmz;
    double denominator = 2.0 * (1.0 - mz * mz);
    
    return numerator / denominator;
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <samplesize> <n>\n", argv[0]);
        return 1;
    }
    
    
    double h_start = atof(argv[1]);
    double h_end = atof(argv[2]);
    
    int samplesize = atoi(argv[3]);
    int n = atoi(argv[4]);

    if (n <= 0 || samplesize <= 0) {
        fprintf(stderr, "Error: n and samplesize must be positive integers\n");
        return 1;
    }

    char filename[100];
    snprintf(filename, sizeof(filename), "data/tfim_susceptibility_data_%d_%d.csv", n, samplesize);

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Unable to create file %s\n", filename);
        return 1;
    }

    fprintf(file, "h,susceptibility\n");


    double h_step = (h_end - h_start) / (samplesize - 1);

    printf("10 Susceptibility values for equally spaced h values:\n");
    printf("h,susceptibility\n");

    for (int i = 0; i < samplesize; i++) {
        double h = h_start + i * h_step;
        double suscept = FiniteSizeSusceptibility(h, n);
        fprintf(file, "%.15e,%.15e\n", h, suscept);
        
        // Print 10 equally spaced values
        if (i % (samplesize / 10) == 0 || i == samplesize - 1) {
            printf("%.15e,%.15e\n", h, suscept);
        }
    }

    fclose(file);
    printf("\nData written to %s\n", filename);

    return 0;
}

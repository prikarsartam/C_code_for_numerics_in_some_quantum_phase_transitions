#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>



long double finite_size_sigma_exp(long double h, int n) {
    const long double PI = 2.0 * acos(0.0);
    long double sum = 0.0;

    for (int m_init = - n/2; m_init <= n/2 - 1 ; m_init++) {
        long double m = (long double)m_init + 0.5;
        long double k = 2.0 * PI * (long double)m / (long double)n;
        long double term = (h - cos(k)) / sqrt(pow(h - cos(k), 2) + pow(sin(k), 2));
        sum += term;
    }
    return sum / (long double)n;
}


long double Derivative_of_finite_size_sigma_exp(long double h, int n) {
    const long double PI = 2.0 * acos(0.0);
    long double sum_n = 0.0;
    for (int m_init = - n/2; m_init <= n/2 - 1 ; m_init++) {
        long double m = (long double)m_init + 0.5;
        long double kval = 2.0 * PI * m / (long double)n;
        long double term_n = (sin(kval) * sin(kval)) / pow((1.0 + h * h - 2.0 * h * cos(kval)), 1.5);
        sum_n += term_n;
    }
    return sum_n / n;
}




/*
long double Derivative_of_finite_size_sigma_exp(long double h, int n) {
 long double delta = 1e-8;
 long double delta_mz = (finite_size_sigma_exp(h + delta, n) - finite_size_sigma_exp(h, n));
 return delta_mz / delta;
 
}
*/



long double FiniteSizeSusceptibility(double h, int n) {
    long double mz = finite_size_sigma_exp(h, n);
    long double dhmz = Derivative_of_finite_size_sigma_exp(h, n);
    
    long double numerator = dhmz * dhmz;
    long double denominator = 2.0 * (1.0 - mz * mz);
    
    return numerator / denominator;
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <samplesize> <n>\n", argv[0]);
        return 1;
    }
    
    
    double h_initial = atof(argv[1]);
    double h_final = atof(argv[2]);
    
    long double h_start = (long double)h_initial;
    long double h_end = (long double)h_final;
    
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


    long double h_step = (h_end - h_start) / ((long double)samplesize - 1);

    printf("10 Susceptibility values for equally spaced h values:\n");
    printf("h,susceptibility\n");

    for (int i = 0; i < samplesize; i++) {
        long double h = h_start + i * h_step;
        long double suscept = FiniteSizeSusceptibility(h, n);
        fprintf(file, "%.30Lf,%.30Lf\n", h, suscept);
        
        // Print 10 equally spaced values
        if (i % (samplesize / 10) == 0 || i == samplesize - 1) {
            printf("%.30Lf,%.30Lf\n", h, suscept);
        }
    }

    fclose(file);
    printf("\nData written to %s\n", filename);

    return 0;
}

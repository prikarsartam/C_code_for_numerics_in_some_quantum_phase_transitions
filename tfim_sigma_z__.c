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
    snprintf(filename, sizeof(filename), "data/tfim_sigma_data_%d_%d.csv", n, samplesize);

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Unable to create file %s\n", filename);
        return 1;
    }

    fprintf(file, "h,sigma\n");


    double h_step = (h_end - h_start) / (double)(samplesize - 1);

    printf("10 Sigma_z values for equally spaced h values:\n");
    printf("h,sigma_z\n");

    for (int i = 0; i < samplesize; i++) {
        double h = h_start + (double)i * h_step;
        double sigma = finite_size_sigma_exp(h, n);
        fprintf(file, "%.15e,%.15e\n", h, sigma);
        
        if (i % (samplesize / 10) == 0) {
            printf("h = %.15e,   sigma_z(h) = %.15e\n", h, sigma);
        }
    }

    fclose(file);
    printf("\nData written to %s\n", filename);

    return 0;
}

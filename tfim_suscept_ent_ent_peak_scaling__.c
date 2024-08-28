#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


#define GOLDEN_RATIO 1.6180339887498948482045868
#define TOLERANCE 1e-20
#define MAX_ITERATIONS 100000





long double finite_size_sigma_exp(long double h, long double gamma, int n) {
    const long double PI = 2.0 * acos(0.0);
    long double sum = 0.0;

    for (int k_init = - n/2; k_init <= n/2 - 1 ; k_init++) {
        long double k = (long double)k_init + 0.5;
        double kval = 2.0 * PI * k / (long double)n;
        long double term = (h - cos(kval)) / sqrt(pow(h - cos(kval), 2) + pow(gamma * sin(kval), 2));
        sum += term;
    }
    return sum / (long double)n;
}



double Derivative_of_finite_size_sigma_exp(long double h, long double gamma, int n) {
    const long double PI = 2.0 * acos(0.0);
    long double sum_n = 0.0;
    for (int k_init = - n/2; k_init <= n/2 - 1 ; k_init++) {
        long double k = k_init + 0.5;
        long double kval = 2.0 * PI * k / (long double)n;
        long double term_n = (pow(gamma * sin(kval) , 2)) / pow( ( pow(h - cos(kval), 2) + pow(gamma * sin(kval), 2) ), 1.5);
        sum_n += term_n;
    }
    return sum_n / (long double)n;
}



/*

double Derivative_of_finite_size_sigma_exp(double h, int n) {
 double delta = 1e-8;
 double delta_mz = (finite_size_sigma_exp(h + delta, n) - finite_size_sigma_exp(h, n));
 return delta_mz / delta;
 
}
*/


long double calculate_g11(long double h, long double gamma, int n) {

    
    long double mz = finite_size_sigma_exp(h, gamma, n);
    long double dhmz = Derivative_of_finite_size_sigma_exp(h, gamma, n);
    
    long double numerator =  pow(dhmz, 2) ;
    long double denominator =  2 * (1 - pow(mz, 2));
    
    return numerator/denominator ;
}


long double parabolic_interpolation(long double (*f)(long double, long double, int), long double x0, long double x1, long double x2, long double gamma, int n) {
    long double f0 = f(x0, gamma, n);
    long double f1 = f(x1, gamma, n);
    long double f2 = f(x2, gamma, n);
    
    long double x = x1 - 0.5 * ((x1 - x0) * (x1 - x0) * (f1 - f2) - (x1 - x2) * (x1 - x2) * (f1 - f0)) /
                          ((x1 - x0) * (f1 - f2) - (x1 - x2) * (f1 - f0));
    
    return x;
}

long double optimize_g11(int n) {
    long double gamma = 1.0;
    long double a = 1- 1/((long double) n), b = 1 + 1/((long double) n); // Assuming the peak is near h=1
    long double x0 = a, x3 = b;

/*
    double x1 = x0 + (x3 - x0) / GOLDEN_RATIO;
    double x2 = x3 - (x3 - x0) / GOLDEN_RATIO;
*/
    
    long double x1 = x0 + (x3 - x0)/5.0;
    long double x2 = x3 - (x3 - x0)/5.0;
    
    
    for (int i = 0; i < MAX_ITERATIONS; i++) {
        long double x_new = parabolic_interpolation(calculate_g11, x0, x1, x3, gamma, n);
        
        if (fabs(x_new - x1) < TOLERANCE) {
            return x_new;
        }
        
        if (x_new > x1) {
            if (calculate_g11(x_new, gamma, n) > calculate_g11(x1, gamma, n)) {
                x0 = x1;
                x1 = x_new;
            } else {
                x3 = x_new;
            }
        } else {
            if (calculate_g11(x_new, gamma, n) > calculate_g11(x1, gamma, n)) {
                x3 = x1;
                x1 = x_new;
            } else {
                x0 = x_new;
            }
        }
        
    }
    
    return x1;
}


long double tfim_size_vs_peak_distance_from_criticality(int n) {
    long double h_critical = optimize_g11(n);
    return fabs(1.0 - h_critical);
}

void calculate_and_store_data(const int* size_list, int size_count, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }
    
    fprintf(file, "System_Size,Distance_from_Criticality\n");
    
    printf("\n\n Calculating the distances\n\n");
    
    for (int i = 0; i < size_count; i++) {
        int size = size_list[i];
        long double distance = tfim_size_vs_peak_distance_from_criticality(size);
        fprintf(file, "%.30Lf,%.30Lf\n", 1/(long double)size, distance); 

        printf("| h_c(%d) - h_c(thermo) |  = %.30Lf\n", size, distance);
    }
    
    fclose(file);
    printf("\n\nData stored in %s\n\n", filename);
}

int main() {

/*
    int size_list[] = { 4, 8, 12, 16, 24, 36, 48, 56, 62, 76, 90, 96, 102, 150, 286, 378, 500, 650, 748, 800, 1000, 2000, 2750, 3500, 5000, 7000, 8000, 9000, 10000, 11500, 12000, 13500, 15000, 17500, 20000, 25000, 50000, 65000, 77800, 85000, 94200, 100000};
    
*/

int size_list[] = {
    4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 
    54, 56, 58, 60, 62, 64, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 
    110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 286, 300, 350, 378, 400, 450, 500, 550, 650, 
    748, 800, 850, 900, 950, 1000, 1250, 1400, 1650, 1800, 2000, 3000, 4000, 5000, 7500, 10000, 20000, 36000, 50000, 75000, 100000,150000, 200000, 250000, 500000, 750000,1000000
};


    
    int size_count = sizeof(size_list) / sizeof(size_list[0]);
    
    calculate_and_store_data(size_list, size_count, "tfim_suscept_peak_data__28072024_1650__inLongDouble.csv");
    
    return 0;
}

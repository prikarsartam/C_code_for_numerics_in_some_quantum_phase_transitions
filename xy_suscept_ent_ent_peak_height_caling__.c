#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


#define TOLERANCE 1e-24
#define MAX_ITERATIONS 10000





long double finite_size_sigma_exp(long double h, long double gamma, int n) {
    const long double PI = 2.0 * acos(0.0);
    long double sum = 0.0;

    for (int k = - (n-1)/2; k <= (n-1)/2 ; k++) {
        double kval = 2.0 * PI * k / (long double)n;
        long double term = (h - cos(kval)) / sqrt(pow(h - cos(kval), 2) + pow(gamma * sin(kval), 2));
        sum += term;
    }
    return sum / (long double)n;
}



double Derivative_of_finite_size_sigma_exp(long double h, long double gamma, int n) {    // In XY --> derivative is w.r.t. gamma
    const long double PI = 2.0 * acos(0.0);
    long double sum_n = 0.0;
    for (int k = - (n-1)/2; k <= (n-1)/2 - 1 ; k++) {
        long double kval = 2.0 * PI * k / (long double)n;
        long double term_n =  ((cos(kval)-h) * gamma * sin(kval) * sin(kval)) / pow((pow((h - cos(kval)), 2) + pow( gamma * sin(kval), 2)), 1.5000000000);
        sum_n += term_n;
    }
    return sum_n / (long double)n;
}



long double calculate_g22(long double h, long double gamma, int n) {

    
    long double mz = finite_size_sigma_exp(h, gamma, n);
    long double dhmz = Derivative_of_finite_size_sigma_exp(h, gamma, n);
    
    long double suscept =  pow(dhmz, 2)/ ( 2 * (1 - pow(mz, 2)));
    
    return suscept ;
}


long double parabolic_interpolation(long double (*f)(long double, long double, int), long double x0, long double x1, long double x2, long double global_h, int n) {
    long double f0 = f(global_h, x0, n);
    long double f1 = f(global_h, x1, n);
    long double f2 = f(global_h, x2, n);
    
    long double x = x1 - 0.5 * ((x1 - x0) * (x1 - x0) * (f1 - f2) - (x1 - x2) * (x1 - x2) * (f1 - f0)) /
                          ((x1 - x0) * (f1 - f2) - (x1 - x2) * (f1 - f0));
    
    return x;
}

long double optimize_g22(int n) {    

    long double global_h = 0.0;
    
    long double a = 0.0, b = 0.5; // Assuming the peak is near gamma=0
    long double x0 = a, x3 = b;

    
    long double x1 = x0 + (x3 - x0)/1.5;
    long double x2 = x3 - (x3 - x0)/1.5;
    
    
    for (int i = 0; i < MAX_ITERATIONS; i++) {
        long double x_new = parabolic_interpolation(calculate_g22, x0, x1, x3, global_h, n);
        
        if (fabs(x_new - x1) < TOLERANCE) {
            return x_new;
        }
        
        if (x_new > x1) {
            if (calculate_g22(global_h, x_new, n) > calculate_g22(global_h, x1, n)) {
                x0 = x1;
                x1 = x_new;
            } else {
                x3 = x_new;
            }
        } else {
            if (calculate_g22(global_h, x_new, n) > calculate_g22(global_h, x1 , n)) {
                x3 = x1;
                x1 = x_new;
            } else {
                x0 = x_new;
            }
        }
        
    }
    
    return x1;
}


long double xy_size_vs_peak_height(int n) {
    long double gamma_critical = optimize_g22(n);
    long double peak_value = calculate_g22(0.0, gamma_critical, n);
    return peak_value;
}

void calculate_and_store_data(const int* size_list, int size_count, const char* filename) {

    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }
    
    fprintf(file, "System_Size,Peak_Height\n");
    
    printf("\n\n Calculating the peak heights\n\n");
    
    for (int i = 0; i < size_count; i++) {
        int size = size_list[i]+1;   // to make all sizes odd
        long double peak_height = xy_size_vs_peak_height(size);
        fprintf(file, "%d,%.48Lf\n", size, peak_height); 

        printf("suscept_peak(%d)  = %.48Lf\n", size, peak_height);
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
    748, 800, 850, 900, 950, 1000, 1250, 1400, 1650, 1800, 2000, 3000, 4000, 5000
};


    
    int size_count = sizeof(size_list) / sizeof(size_list[0]);
    
    calculate_and_store_data(size_list, size_count, "data/xy_suscept_peak_height_data__30072024_0150__inLongDouble.csv");
    
    return 0;
}

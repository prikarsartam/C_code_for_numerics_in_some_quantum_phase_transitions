#include<stdio.h>
#include<math.h>

int main() {
    printf("%ld\n", sizeof(double));      // some compilers print 8
    printf("%ld\n", sizeof(long double)); // some compilers print 16
    return 0;
}

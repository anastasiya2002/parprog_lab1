#include <stdio.h>

double t_max = 50;
double x_max = 60;
double t_step = 0.1;
double x_step = 0.1;

double func(double t, double x){
    return x*t;
}

double fi(double x){
    return x*x*x / 12;
}

double ksi(double t){
    return t*t*t / 12;
}

void shema(int k, int m, double** u) {
    if (m != x_max / x_step - 1) {
        // Явная центральная трехточечная схема
        u[k + 1][m] = func(k, m) * t_step + 1/2.*(u[k][m+1] + u[k][m-1]) - t_step/2./x_step * (u[k][m+1] - u[k][m-1]);
    } else {
        // Явный левый уголок
        u[k + 1][m] = func(k, m) * t_step + u[k][m] - t_step/x_step * (u[k][m] - u[k][m-1]);
    }
}


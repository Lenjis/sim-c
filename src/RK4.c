#include "Sim.h"

void rk4(void (*fun)(), double t, double x[], double u[], int dim, double T,
         double x2[], double *t2) {
    int idx;
    double xd[30], q[30], xx[30], x_a[30], tt;

    (*fun)(t, x, u, xd, dim);
    for (idx = 0; idx < dim; idx++) x_a[idx] = xd[idx] * T;  // K1
    for (idx = 0; idx < dim; idx++)
        xx[idx] = x[idx] + x_a[idx] / 2;  // xx = x(t) + K1/2
    tt = t + T / 2;

    (*fun)(tt, xx, u, xd, dim);  // xd = f(x+K1/2, t+T/2)
    for (idx = 0; idx < dim; idx++) q[idx] = xd[idx] * T;  // q=K2
    for (idx = 0; idx < dim; idx++)
        xx[idx] = x[idx] + q[idx] / 2;  // xx=x(t)+K2/2
    for (idx = 0; idx < dim; idx++)
        x_a[idx] = x_a[idx] + 2 * q[idx];  // K1+2*K2

    (*fun)(tt, xx, u, xd, dim);  // xd = f(x+K2/2, t+T/2)
    for (idx = 0; idx < dim; idx++) q[idx] = xd[idx] * T;       // q=K3
    for (idx = 0; idx < dim; idx++) xx[idx] = x[idx] + q[idx];  // xx=x(t)+K3
    for (idx = 0; idx < dim; idx++)
        x_a[idx] = x_a[idx] + 2 * q[idx];  // K1+2K2+2K3

    tt = t + T;
    (*fun)(tt, xx, u, xd, dim);  // xd = f(x+K3/2, t+T)
    for (idx = 0; idx < dim; idx++)
        x2[idx] =
            x[idx] + (x_a[idx] + xd[idx] * T) / 6;  // x(t) +[K1+2K2+2K3+K4]/6
    *t2 = t + T;
}

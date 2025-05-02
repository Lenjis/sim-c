#include "model.h"  // Library functions declared here
#include "matrix.h"
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

// [Aircraft parameters]
static double SA = 1.3536,               // Wing area [m^2]
    b = 3.2,                             // Wingspan [m]
    cbar = 0.423,                        // Mean aerodynamic chord [m]
    Ixx = 1.71, Iyy = 5.54, Izz = 4.15,  // Moments of inertia [kg*m^2]
    Ixy = 0, Ixz = 0, Iyz = 0,           // Products of inertia [kg*m^2]
    g = 9.81, mass = 17;                 // Gravitational acceleration and mass

static void uav_density(double H, double VT, double *ru, double *mach);
static double uav_CD(double alpha_deg);
static double uav_CL0(double alpha_deg);
static double uav_CM(double alpha_deg);
static double uav_interp1(double *A, double idx[], int Len1, double xi);
static double uav_interp2(double *A, double idx[], int Len1, double xi, double Dim2[],
                          int Len2, double X2);
static double uav_interp3(double *A, double idx[], int Len1, double xi, double Dim2[],
                          int Len2, double X2, double Dim3[], int Len3, double X3);
static int uav_find(double A[], double X, int len);

double ac_dH;

void model6dof(double t, double x[], double u[], double dx[], int dim) {
    double Vt, alpha, beta, phi, theta, psi, P, Q, R, PN, PE, H;
    double alpha_deg, beta_deg;
    double dVt, dalpha, dbeta, dphi, dtheta, dpsi, dP, dQ, dR, dPN, dPE, dH;
    double salpha, sbeta, sphi, stheta, spsi, calpha, cbeta, cphi, ctheta, cpsi;
    double ele, ail, rud, eng, Pow;
    double U, V, W;
    double dU, dV, dW;
    double ru, mach, qs;
    double CD, CL0, CM0, D, L, Y;
    double Lbar, M, N;

    const double CL_ele = 0.00636;
    const double CY_beta = -0.00909;
    const double CR_beta = -0.00600, CR_ail = -0.003618, CR_rud = 0.000144,
                 CR_P = -0.52568, CR_R = 0.01832;
    const double CM_ele = -0.02052, CM_Q = -9.3136 / 3.0, CM_dalpha = -4.0258;
    const double CN_beta = 0.00235, CN_ail = 0.000132, CN_rud = -0.00111,
                 CN_P = 0.01792, CN_R = -0.15844;

    Vt = x[0];
    alpha = x[1];
    beta = x[2];
    phi = x[3];
    theta = x[4];
    psi = x[5];
    P = x[6];
    Q = x[7];
    R = x[8];
    PN = x[9];
    PE = x[10];
    H = x[11];

    // control input
    ele = u[0];  // elevator deflection angle [deg]
    ail = u[1];  // aileron  deflection angle [deg]
    rud = u[2];  // aileron  deflection angle [deg]
    eng = u[3];  // engine input

    uav_density(H, Vt, &ru, &mach);  // [air density] [mach number]
    qs = SA * (ru * Vt * Vt / 2);    // [Dynamic pressure](kg/m^2)

    alpha_deg = alpha * 180.0 / M_PI;
    beta_deg = beta * 180.0 / M_PI;

    salpha = sin(alpha);
    calpha = cos(alpha);
    sbeta = sin(beta);
    cbeta = cos(beta);
    sphi = sin(phi);
    cphi = cos(phi);
    stheta = sin(theta);
    ctheta = cos(theta);
    spsi = sin(psi);
    cpsi = cos(psi);

    U = Vt * calpha * cbeta;
    V = Vt * sbeta;
    W = Vt * salpha * cbeta;
    double uvw[3] = {U, V, W};  // [u,v,w]  ----- body axis velocity

    Pow = eng / 100 * (mass * g / 4.0);

    double J[3][3] = {
        {Ixx, 0, -Ixz}, {0, Iyy, 0}, {-Ixz, 0, Izz}};  // [inertia matrix]
    double J_inv[3][3] = {
        {-Izz / (Ixz * Ixz - Ixx * Izz), 0, -Ixz / (Ixz * Ixz - Ixx * Izz)},
        {0, 1 / Iyy, 0},
        {-Ixz / (Ixz * Ixz - Ixx * Izz), 0, -Ixx / (Ixz * Ixz - Ixx * Izz)}};
    double S[3][3] = {
        {calpha * cbeta, -calpha * sbeta, -salpha},
        {sbeta, cbeta, 0},
        {salpha * cbeta, -salpha * sbeta, calpha}};  // [rotation matrix]
    double B[3][3] = {
        {ctheta * cpsi, ctheta * spsi, -stheta},
        {sphi * stheta * cpsi - cphi * spsi, sphi * stheta * spsi + cphi * cpsi,
         sphi * ctheta},
        {cphi * stheta * cpsi + sphi * spsi, cphi * stheta * spsi - sphi * cpsi,
         cphi * ctheta}};           // [rotation matrix]
    double pqr[3] = {P, Q, R};      // [roll rate] [yaw rate] [pitch rate]
    double Pow_v[3] = {Pow, 0, 0};  // [thrust vector]

    qs = SA * ru * Vt * Vt / 2;  // [Dynamic pressure](kg/m^2)

    CD = uav_CD(alpha_deg);  // [drag coefficient]
    CL0 = uav_CL0(alpha_deg);
    CM0 = uav_CM(alpha_deg);  // [moment coefficient]

    D = qs * CD;
    Y = qs * CY_beta * beta_deg;
    L = qs * (CL0 + CL_ele * ele);

    double F0_v[3] = {-D, Y, -L};  // [force vector]
    double F_v[3];
    double Fxyz[3];
    multiply(&S[0][0], F0_v, F_v, 3, 3, 3, 1);  // F = S * [-D; Y; -L]
    add(Pow_v, F_v, Fxyz, 3, 1);  // Fxyz = [Pow; 0; 0] + S * [-D; Y; -L];

    // duvw = Fxyz / mass - cross(pqr, uvw) + g * [-stheta; sphi * ctheta;
    // cphi * ctheta];
    double duvw[3] = {Fxyz[0] / mass, Fxyz[1] / mass, Fxyz[2] / mass};
    vector3 pqr_v = {P, Q, R};
    vector3 uvw_v = {U, V, W};
    vector3 cv1 = cross(pqr_v, uvw_v);  // cross(pqr, uvw)
    double v1[3] = {cv1.x, cv1.y, cv1.z};
    subtract(duvw, v1, duvw, 3, 1);
    double c2[3] = {-stheta * g, sphi * ctheta * g, cphi * ctheta * g};
    add(duvw, c2, duvw, 3, 1);
    dU = duvw[0];
    dV = duvw[1];
    dW = duvw[2];

    dVt = (U * dU + V * dV + W * dW) / Vt;
    dbeta = (dV * Vt - V * dVt) / (Vt * Vt * cbeta);
    dalpha = (U * dW - W * dU) / (U * U + W * W);

    Lbar = qs * b *
           (CR_beta * beta_deg + CR_ail * ail + CR_rud * rud +
            (CR_P * P + CR_R * R) * b / Vt / 2.0);
    M = qs * cbar * (CM0 + CM_ele * ele + CM_Q * Q * cbar / Vt / 2.0);
    N = qs * b *
        (CN_beta * beta_deg + CN_ail * ail + CN_rud * rud +
         (CN_P * P + CN_R * R) * b / Vt / 2.0);

    // dpqr = inv(J) * (-cross(pqr, (J * pqr)) + [Lbar; M; N]);
    double dpqr[3];
    double v3[3];
    multiply(&J[0][0], pqr, v3, 3, 3, 3, 1);  // v3 = J*pqr
    vector3 v3_v = {v3[0], v3[1], v3[2]};
    vector3 c4 = cross(pqr_v, v3_v);
    double v4[3] = {-c4.x, -c4.y, -c4.z};  // v4 = -pqr x (J*pqr)
    multiply(&J_inv[0][0], v4, dpqr, 3, 3, 3,
             1);  // J_inv * -(pqr x (J*pqr))
    double v5[3] = {Lbar, M, N};
    add(dpqr, v5, dpqr, 3, 1);  // dpqr = J_inv*(pqr x (J*pqr)) + [Lbar,M,N]
    dP = dpqr[0];
    dQ = dpqr[1];
    dR = dpqr[2];

    dphi = P + (stheta / ctheta) * (Q * sphi + R * cphi);
    dtheta = Q * cphi - R * sphi;
    dpsi = (Q * sphi + R * cphi) / ctheta;

    // dVe = B' * uvw;
    double B_t[3][3];
    transpose(&B[0][0], &B_t[0][0], 3, 3);  // B_t = B^T
    double dVe[3];
    multiply(&B_t[0][0], uvw, dVe, 3, 3, 3, 1);  // dVe = B^T * uvw
    dPN = dVe[0];
    dPE = dVe[1];
    dH = -dVe[2];  // [dVe] = [dPN,dPE,dH]
    ac_dH = dH;

    dx[0] = dVt;
    dx[1] = dalpha;
    dx[2] = dbeta;
    dx[3] = dphi;
    dx[4] = dtheta;
    dx[5] = dpsi;
    dx[6] = dP;
    dx[7] = dQ;
    dx[8] = dR;
    dx[9] = dPN;
    dx[10] = dPE;
    dx[11] = dH;
}

// Return air density and mach number
//===================================
static void uav_density(double H, double VT, double *ru, double *mach) {
    double T, PP, RR, TS, sonic;
    const double K = 34.163195;
    const double RL = 1.225;
    const double C1 = 0.001;
    const double AL = 340.294;
    const double BT = 1.458E-06;

    H = C1 * H / (1 + C1 * H / 6356.766);

    if (H < 11) {
        T = 288.15 - 6.5 * H;
        PP = pow(288.15 / T, -K / 6.5);
    } else if (H < 20) {
        T = 216.65;
        PP = 0.22336 * exp(-K * (H - 11) / 216.65);
    } else if (H < 32) {
        T = 216.65 + (H - 20);
        PP = 0.054032 * pow(216.65 / T, K);
    } else if (H < 47) {
        T = 228.65 + 2.8 * (H - 32);
        PP = 0.0085666 * pow(228.65 / T, K / 2.8);
    } else if (H < 51) {
        T = 270.65;
        PP = 0.0010945 * exp(-K * (H - 47) / 270.65);
    } else if (H < 71) {
        T = 270.65 - 2.8 * (H - 51);
        PP = 0.00066063 * pow(270.65 / T, -K / 2.8);
    } else {
        T = 214.65 - 2 * (H - 71);
        PP = 3.9046e-05 * pow(214.65 / T, -K / 2);
    }

    RR = PP / (T / 288.15);
    TS = T / 288.15;
    sonic = AL * sqrt(TS);
    *ru = RL * RR;
    *mach = VT / sonic;
}

// return drag coefficient
//========================
static double uav_CD(double alpha_deg) {
    static double IDX_alpha[9] = {-4.0, -2.0, 0.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20.0};
    static double TBL_CD[9] = {0.026, 0.024, 0.024, 0.028, 0.036,
                               0.061, 0.102, 0.141, 0.173};
    return uav_interp1(TBL_CD, IDX_alpha, 9, alpha_deg);
}

// return lift coefficient
//========================
static double uav_CL0(double alpha_deg) {
    static double IDX_alpha[9] = {-4.0, -2.0, 0.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20.0};
    static double TBL_CL0[9] = {-0.219, -0.04, 0.139, 0.299, 0.455,
                                0.766,  1.083, 1.409, 1.743};
    return uav_interp1(TBL_CL0, IDX_alpha, 9, alpha_deg);
}

// return moment coefficient
//==========================
static double uav_CM(double alpha_deg) {
    static double IDX_alpha[9] = {-4.0, -2.0, 0.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20.0};
    static double TBL_CM0[9] = {0.1161,  0.0777,  0.0393,  0.0009, -0.0375,
                                -0.0759, -0.1527, -0.2295, -0.3063};
    return uav_interp1(TBL_CM0, IDX_alpha, 9, alpha_deg);
}

// 1D interpolation
//=================
static double uav_interp1(double *A, double idx[], int Len1, double xi) {
    int r;
    double DA, Y;

    r = uav_find(idx, xi, Len1);

    DA = (xi - idx[r]) / (idx[r + 1] - idx[r]);
    Y = A[r] + (A[r + 1] - A[r]) * DA;
    return (Y);
}

// find the index of the nearest point
//====================================
static int uav_find(double A[], double X, int len) {
    int result = 0;
    int i, P_Start, P_End;

    if (A[0] < A[len - 1]) { /*[数组顺序排列]*/
        P_Start = 0;
        P_End = len;
        if (X < A[0]) {
            result = 0;
        } else if (X < A[len - 1]) {
            if (X > A[len / 2])
                P_Start = len / 2;
            else
                P_End = len / 2 + 1;

            for (i = P_Start; i < P_End; i++) {
                if (X < A[i]) {
                    result = i - 1;
                    break;
                }
            }
        } else {
            result = len - 2;
        }
    } else { /*[数组逆序排列]*/
        P_Start = len - 1;
        P_End = 0;
        if (X >= A[0])
            result = 0;
        else if (X > A[len - 1]) {
            if (X > A[len / 2])
                P_Start = len / 2 + 1;
            else
                P_End = len / 2;
            for (i = P_Start; i >= P_End; i--) {
                if (X < A[i]) {
                    result = i;
                    break;
                }
            }
        } else
            result = len - 2;
    }
    return result;
}

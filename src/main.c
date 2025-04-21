#include "Sim.h"
#include "model.h"
#include "matrix.h"
#include <windows.h>
#include <mmsystem.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#pragma comment(lib, "winmm.lib")

#define KP_H 1.2
#define KI_H 0.1
#define KD_H 0

#define KP_THETA 0.7
#define KI_THETA 0.3
#define KD_THETA 0.1

#define KA_PHI 0.7
#define KA_P 0.1

void (*aircraft)();
double T, x[13], u[3], t = 0;
int sim_step, sim_status;
short flag_Stop = 1;
short ctrl_state = 0;
double theta_cmd, theta_var, phi_cmd, phi_var, H_cmd, H_out = 0;

void ctrl_alt(void) {
    const double Kp_H = KP_H, Ki_H = KI_H, Kd_H = KD_H, dt = 0.01;
    static double H_i = 0, H_e = 0, H_prev = 0, H_d;

    H_e = H_cmd - ac_H;  // 高度误差

    H_d = (ac_H - H_prev) / dt;  // 高度导数
    H_prev = ac_H;

    if (H_e > 20)
        H_i += 20 * dt;
    else if (H_e < -20)
        H_i += -20 * dt;
    else
        H_i += H_e * dt;  // 积分项, dt=0.01s

    // if (H_i > 20)
    //     H_i = 20;
    // else if (H_i < -20)
    //     H_i = -20;

    H_out = Kp_H * H_e + Ki_H * H_i + Kd_H * H_d;  // 高度控制
}

// void ctrl_long(void) {
//     static double Kp_theta = 1, Ki_theta = 0.1, K_Q = 0.3;
//     static double theta_i = 0, theta_e = 0;

//     theta_e = theta_cmd - ac_theta * Rad2Deg;  //[deg]

//     if (theta_e > 10)
//         theta_i += 10 * 0.01;
//     else if (theta_e < -10)
//         theta_i += -10 * 0.01;
//     else
//         theta_i += theta_e * 0.01;  // 积分项, dt=0.01s

//     if (theta_i > 10)
//         theta_i = 10;
//     else if (theta_i < -10)
//         theta_i = -10;

//     ac_ele = K_Q * ac_Q * Rad2Deg - Kp_theta * (theta_e + H_out) - Ki_theta * theta_i;  // K_Q极性相反
// }

void ctrl_long(void) {  // incremental PID
    const double Kp_theta = KP_THETA, Ki_theta = KI_THETA, Kd_theta = KD_THETA, dt = 0.01;
    static double theta_e = 0, theta_e1 = 0, theta_e2 = 0;  // 当前、上一次、上上次误差
    static double du = 0;
    static double ac_ele_last = 0;

    theta_e = theta_cmd + H_out - ac_theta * Rad2Deg;

    du =
        Kp_theta * (theta_e - theta_e1) + Ki_theta * theta_e * dt + Kd_theta * (theta_e - 2 * theta_e1 + theta_e2) / dt;

    ac_ele -= du;  // 舵量输入相反

    theta_e2 = theta_e1;
    theta_e1 = theta_e;

    // 加上角速率反馈项
    // ac_ele -= 0.3 * ac_Q * Rad2Deg;
}

// /*飞行器纵向控制*/
// void ctrl_long(void) {
//     static double Ke_theta = 0.45, Ke_Q = 0.25, Ke_phi = 0.05;

//     ac_ele = Ke_theta * (ac_theta * Rad2Deg - theta_cmd) + Ke_Q * ac_Q * Rad2Deg;
// }

/*飞行器横侧向控制*/
void ctrl_late(void) {
    static double Ka_phi = KA_PHI, Ka_P = KA_P;

    ac_ail = Ka_phi * (ac_phi * Rad2Deg - phi_cmd) + Ka_P * ac_P * Rad2Deg;
}

/*飞行器控制模块*/
void ctrl_level(void) {
    if (t > 40) flag_Stop = 0;
    phi_cmd = 0;
    theta_cmd = 1.1190407073;
    H_cmd = 500;
    ctrl_alt();
    ctrl_long();
    // ctrl_late();
}

void ctrl_rectangular(void) {
    // during:dpsi=psi_cmd-psi_deg;
    // while(dpsi>180) dpsi=dpsi-360;end
    // while(dpsi<-180) dpsi=dpsi+360;end
    // phi_cmd=1.0*dpsi;
    // if(phi_cmd>45) phi_cmd=45;end
    // if(phi_cmd<-45) phi_cmd=-45;end
    static double dpsi = 0, psi_cmd;

    switch (ctrl_state) {
        case 0:
            psi_cmd = 0;
            if (ac_PN >= 500) ctrl_state++;
            break;
        case 1:
            psi_cmd = 90;
            if (ac_PE >= 500) ctrl_state++;
            break;
        case 2:
            psi_cmd = 180;
            if (ac_PN <= 0) ctrl_state++;
            break;
        case 3:
            psi_cmd = 270;
            if (ac_PE <= 0) ctrl_state = 0;
            break;
        default:
            break;
    }
    if (t > 200) flag_Stop = 0;
    dpsi = psi_cmd - ac_psi * Rad2Deg;
    while (dpsi > 180) dpsi = dpsi - 360;
    while (dpsi < -180) dpsi = dpsi + 360;
    phi_cmd = 1.0 * dpsi;
    if (phi_cmd > 45) phi_cmd = 45;
    if (phi_cmd < -45) phi_cmd = -45;
    theta_cmd = 1.1190407073;
    H_cmd = 500;
    ctrl_alt();
    ctrl_long();
    ctrl_late();
}

void ctrl_approach(void) {}

/*飞行器模型解算模块，无需看懂*/
void simu_run(void) {
    // static double g = 9.81;
    // static double stheta, ctheta, sphi, cphi, spsi, cpsi;

    u[0] = ac_ele;
    u[1] = ac_ail;
    u[2] = ac_rud;
    u[3] = ac_eng;

    x[0] = ac_Vt;
    x[1] = ac_alpha;
    x[2] = ac_beta;
    x[3] = ac_phi;
    x[4] = ac_theta;
    x[5] = ac_psi;
    x[6] = ac_P;
    x[7] = ac_Q;
    x[8] = ac_R;
    x[9] = ac_PN;
    x[10] = ac_PE;
    x[11] = ac_H;

    rk4(aircraft, t, x, u, 12, T, x, &t);

    ac_Vt = x[0];
    ac_alpha = x[1];
    ac_beta = x[2];
    ac_phi = x[3];
    ac_theta = x[4];
    ac_psi = x[5];
    ac_P = x[6];
    ac_Q = x[7];
    ac_R = x[8];
    ac_PN = x[9];
    ac_PE = x[10];
    ac_H = x[11];

    ac_track = (ac_psi)*Rad2Deg;
    if (ac_track < 0.0) ac_track = ac_track + 360.0;
    if (ac_track >= 360.0) ac_track = ac_track - 360.0;
    ac_track = ac_track / Rad2Deg;
}

/*飞行器模型解算初始化，无需看懂*/
void simu_init(void) {
    ac_Vt = 30.0;
    ac_alpha = 1.1190407073 / Rad2Deg;
    ac_beta = 0 / Rad2Deg;
    ac_phi = 0.0 / Rad2Deg;
    ac_theta = 1.1190407073 / Rad2Deg;
    ac_psi = 0.0 / Rad2Deg;
    ac_P = 0 / Rad2Deg;
    ac_Q = 0 / Rad2Deg;
    ac_R = 0 / Rad2Deg;
    ac_PN = 0.0 / Rad2Deg;
    ac_PE = 0.0 / Rad2Deg;
    ac_H = 500;

    sim_step = 5;  // 5ms
    sim_status = 0;
    T = sim_step / 1000.0f;
    t = 0;
    aircraft = model6dof;

    ac_ele = 0.868149045790946;
    ac_ail = 0.0;
    ac_rud = 0.0;
    ac_eng = 44.753409573328163;
}

void CALLBACK Timerdefine(UINT uTimerID, UINT uMsg, DWORD_PTR dwUser, DWORD_PTR dw1, DWORD_PTR dw2) {
    static short cnt = 0;
    cnt++;
    cnt %= 100;

    simu_run(); /*无人机模型解算 周期5ms*/

    if (cnt % 2 == 1) {
        ctrl_rectangular(); /*简单的飞行控制*/
    }
}

void main(void) {
    MMRESULT idtimer;
    FILE *fp;
    static short count = 0;

    simu_init();

    idtimer = timeSetEvent(5, 5, Timerdefine, 0, TIME_PERIODIC);  // 5ms中断一次

    fp = fopen("simu.txt", "w");  // openfile
    fprintf(fp, "t\tVt\tphi\ttheta\tpsi\tPN\tPE\tH\tele\tail\trud\teng\n");

    printf("Running simulation! \n");
    while (flag_Stop) {
        Sleep(5);
        count++;
        count %= 200;
        if (count == 0) {
            printf("t: %6.2lf Vt: %6.2lf phi: %6.2lf theta: %6.2lf psi: %6.2lf PN: %8.2lf PE: %8.2lf H: %6.2lf", t,
                   ac_Vt, ac_phi * Rad2Deg, ac_theta * Rad2Deg, ac_psi * Rad2Deg, ac_PN, ac_PE, ac_H);
            printf(" | ele: %6.2lf ail: %6.2lf rud: %6.2lf eng: %6.2lf\n", ac_ele, ac_ail, ac_rud, ac_eng);
        }
        fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, ac_Vt, ac_phi * Rad2Deg,
                ac_theta * Rad2Deg, ac_psi * Rad2Deg, ac_PN, ac_PE, ac_H, ac_ele, ac_ail, ac_rud, ac_eng);
    };
    fclose(fp);
}

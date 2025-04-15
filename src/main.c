#include "Sim.h"
#include "model.h"
#include "matrix.h"
#include "UAVCtrl.h"
#include <windows.h>
#include <mmsystem.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#pragma comment(lib, "winmm.lib")

void (*aircraft)();
double T, x[13], u[3], t = 0;
int sim_step, sim_status;
short flag_Stop = 1;
short ctrl_state = 0;
double theta_cmd, theta_var, phi_cmd, phi_var;
double H_cmd, H_out, H_i = 0, H_e, H_et, H_d;

void ctrl_alt(void) {
    static double Ke_theta = 1, Ke_Q = 0.3, Ke_phi = 0.05;
    static double Kp_H = 3, Ki_H = 1, Kd_H = 1.5;

    H_e = H_cmd - ac_H;  // 高度误差

    if (H_e > 0.5)
        H_i += 0.5 * T;
    else if (H_e < -0.5)
        H_i += -0.5 * T;
    else
        H_i += H_e * T;  // 积分项, dt=0.01s

    H_d = (H_e - H_et) / T;  // 微分项

    H_out = Kp_H * (H_e) + Ki_H * H_i + Kd_H * H_d;  // 高度控制
    ac_ele = Ke_Q * ac_Q * Rad2Deg - Ke_theta * (theta_cmd - ac_theta * Rad2Deg) - H_out;
}

/*飞行器纵向控制*/
void ctrl_long(void) {
    static double Ke_theta = 0.45, Ke_Q = 0.25, Ke_phi = 0.05;

    ac_ele = Ke_theta * (ac_theta * Rad2Deg - theta_cmd) + Ke_Q * ac_Q * Rad2Deg;
}

/*飞行器横侧向控制*/
void ctrl_late(void) {
    static double Ka_phi = 0.5, Ka_P = 0.1;

    ac_ail = Ka_phi * (ac_phi * Rad2Deg - phi_cmd) + Ka_P * ac_P * Rad2Deg;
}

/*飞行器控制模块*/
void ctrl_task(void) {
    // switch (ctrl_state) {
    //     case 0:
    //         theta_cmd = 2;
    //         phi_cmd = 0;
    //         if (t >= 10) ctrl_state++;
    //         break;
    //     case 1:
    //         theta_cmd = 2;  // t*3.14/5- 周期T=10s
    //         phi_cmd = 30;
    //         if (t >= 20) ctrl_state++;
    //         break;
    //     case 2:
    //         theta_cmd = 2;  // t*3.14/5- 周期T=10s
    //         phi_cmd = 0;
    //         if (t >= 30) ctrl_state++;
    //         break;
    //     case 3:
    //         theta_cmd = 2;  // t*3.14/5- 周期T=10s
    //         phi_cmd = -30;
    //         if (t >= 40) ctrl_state++;
    //         break;
    //     case 4:
    //         flag_Stop = 0;
    //         break;
    //     default:
    //         break;
    // }
    // if (t > 20) flag_Stop = 0;
    phi_cmd = 0;
    theta_cmd = 0.9997;
    H_cmd = 500;
    // ctrl_long();
    // ctrl_late();
    // ctrl_alt();
}

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
    ac_alpha = 1.119040707310974 / Rad2Deg;
    ac_beta = 0 / Rad2Deg;
    ac_phi = 0.0 / Rad2Deg;
    ac_theta = 1.119040707309775 / Rad2Deg;
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
    ac_eng = 44.7534;
}

void CALLBACK Timerdefine(UINT uDelay, UINT uMsg, DWORD dwUser, DWORD dw1, DWORD dw2) {
    static short cnt = 0;
    cnt++;
    cnt %= 100;

    simu_run(); /*无人机模型解算 周期5ms*/

    if (cnt % 2 == 1) {
        // ctrl_task(); /*简单的飞行控制*/
    }
}

void main(void) {
    MMRESULT idtimer;
    FILE *fp;
    static short count = 0;

    simu_init();

    idtimer = timeSetEvent(5, 5, Timerdefine, 0, TIME_PERIODIC);  // 5ms中断一次

    fp = fopen("simu.txt", "w");  // openfile

    printf("Running simulation! \n");
    while (flag_Stop) {
        Sleep(100);
        count++;
        count %= 10;
        if (count == 1) {
            printf(
                "t: %6.2lf |Vt: %8.4lf |phi: %8.4lf "
                "|theta: %8.4lf "
                "|psi: %8.4lf |PN: %8.4lf |H: %8.4lf\n",
                t, ac_Vt, ac_phi * Rad2Deg, ac_theta * Rad2Deg, ac_psi * Rad2Deg, ac_PN, ac_H);
        }
        fprintf(fp,
                "t: %8.4lf |Vt: %8.4lf |phi: %8.4lf "
                "|theta: %8.4lf "
                "|psi: %8.4lf |PN: %8.4lf |H: %8.4lf\n",
                t, ac_Vt, ac_phi * Rad2Deg, ac_theta * Rad2Deg, ac_psi * Rad2Deg, ac_PN, ac_H);
    };
    fclose(fp);
}

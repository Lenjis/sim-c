#include "Sim.h"
#include "model.h"
#include "matrix.h"
#include <windows.h>
#include <mmsystem.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#pragma comment(lib, "winmm.lib")

#define DT 0.01

#define KP_H 1.2
#define KI_H 0.1
#define KD_H 0

#define KP_THETA 0.7
#define KI_THETA 0.3
#define KD_THETA 0.1

#define KP_PHI 3
#define KI_PHI 0.5
#define KD_PHI 0.1

#define KP_SPEED 1000
#define KI_SPEED 1500
#define KD_SPEED 0

void (*aircraft)();
double T, x[13], u[3], t = 0;
int sim_step, sim_status;
short flag_Stop = 1;
short ctrl_state = 0;
double theta_cmd = 0, phi_cmd = 0, H_cmd = 0, speed_cmd = 0, H_out = 0;

void ctrl_alt(void) {
    static double H_i = 0, H_e = 0, H_prev = 0, H_d;

    H_e = H_cmd - ac_H;  // 高度误差

    H_d = (ac_H - H_prev) / DT;  // 高度导数
    H_prev = ac_H;

    if (H_e > 20)
        H_i += 20 * DT;
    else if (H_e < -20)
        H_i += -20 * DT;
    else
        H_i += H_e * DT;  // 积分项, dt=0.01s

    // if (H_i > 20)
    //     H_i = 20;
    // else if (H_i < -20)
    //     H_i = -20;

    H_out = KP_H * H_e + KI_H * H_i + KD_H * H_d;  // 高度控制
}

/*飞行器纵向控制*/
void ctrl_long(void) {                                      // incremental PID
    static double theta_e = 0, theta_e1 = 0, theta_e2 = 0;  // 当前、上一次、上上次误差
    static double du = 0;

    theta_e = theta_cmd + H_out - ac_theta * Rad2Deg;

    du =
        KP_THETA * (theta_e - theta_e1) + KI_THETA * theta_e * DT + KD_THETA * (theta_e - 2 * theta_e1 + theta_e2) / DT;

    ac_ele -= du;  // 舵量输入相反

    theta_e2 = theta_e1;
    theta_e1 = theta_e;
}

void ctrl_speed() {
    static double speed = 0, speed_e1 = 0, speed_e2 = 0;  // 当前、上一次、上上次误差
    static double du = 0;

    speed = speed_cmd - ac_Vt;

    du = KP_SPEED * (speed - speed_e1) + KI_SPEED * speed * DT + KD_SPEED * (speed - 2 * speed_e1 + speed_e2) / DT;

    ac_eng += du;  // 舵量输入相反

    if (ac_eng > 100) ac_eng = 100;
    if (ac_eng < 0) ac_eng = 0;

    speed_e2 = speed_e1;
    speed_e1 = speed;
}

/*飞行器横侧向控制*/
void ctrl_late(void) {
    const double Kp_phi = KP_PHI, Ki_phi = KI_PHI, Kd_phi = KD_PHI;
    static double phi_e = 0, phi_e1 = 0, phi_e2 = 0;  // 当前、上一次、上上次误差
    static double du = 0;

    // ac_ail = Kp_phi * (ac_phi * Rad2Deg - phi_cmd) + Kd_phi * ac_P * Rad2Deg;

    phi_e = phi_cmd - ac_phi * Rad2Deg;

    du = Kp_phi * (phi_e - phi_e1) + Ki_phi * phi_e * DT + Kd_phi * (phi_e - 2 * phi_e1 + phi_e2) / DT;

    ac_ail -= du;  // 舵量输入相反

    phi_e2 = phi_e1;
    phi_e1 = phi_e;
}

/*水平巡航*/
void ctrl_level(void) {
    if (t > 20) flag_Stop = 0;
    speed_cmd = 31;
    theta_cmd = 1.1190407073;
    H_cmd = 500;
    ctrl_alt();
    ctrl_long();
    ctrl_late();
    ctrl_speed();
}

/*矩形轨迹巡航*/
void ctrl_rectangular(void) {
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

/*进场着陆*/
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
        ctrl_level(); /*简单的飞行控制*/
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
        Sleep(10);
        count++;
        if (count % 100 == 0) {
            printf("t: %6.2lf Vt: %6.2lf phi: %6.2lf theta: %6.2lf psi: %6.2lf PN: %8.2lf PE: %8.2lf H: %6.2lf", t,
                   ac_Vt, ac_phi * Rad2Deg, ac_theta * Rad2Deg, ac_psi * Rad2Deg, ac_PN, ac_PE, ac_H);
            printf(" | ele: %6.2lf ail: %6.2lf rud: %6.2lf eng: %6.2lf\n", ac_ele, ac_ail, ac_rud, ac_eng);
        }
        fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, ac_Vt, ac_phi * Rad2Deg,
                ac_theta * Rad2Deg, ac_psi * Rad2Deg, ac_PN, ac_PE, ac_H, ac_ele, ac_ail, ac_rud, ac_eng);
    };
    fclose(fp);
}

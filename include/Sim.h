#ifndef _SIM_H_
#define _SIM_H_

typedef unsigned char BYTE;
typedef unsigned short WORD;
typedef unsigned long DWORD;
typedef unsigned int UINT;

typedef struct {
    int port;          //[串口号]
    short Parity;      //[奇偶校验]
    short StopBit;     //[停止位]
    short WordLength;  //[字长]
    int BAUD;          //[波特率]
} SensorStruc;

SensorStruc ADIS16455;

#define Rad2Deg 57.295779513082323

// 飞机模型仿真参数定义
double ac_Vt,                  //[m/s]
    ac_alpha, ac_beta,         //[rad]
    ac_P, ac_Q, ac_R,          //[rad/s]
    ac_phi, ac_theta, ac_psi,  //[rad]
    ac_heading, ac_track,      //[rad]
    ac_PN, ac_PE, ac_H,        //[m]
    ac_lon,  //[rad]
    ac_lat;  //[rad]

double ac_ele,  //[deg]
    ac_ail,     //[deg]
    ac_rud,     //[deg]
    ac_eng;     //[0-100%]

void rk4(void (*fun)(), double t, double x[], double u[], int dim, double T, double x2[], double *t2);

#endif
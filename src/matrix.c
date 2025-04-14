#include "matrix.h"
#include <stdio.h>

void add(double *a, double *b, double *c, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            c[i * col + j] = a[i * col + j] + b[i * col + j];
        }
    }
}

void subtract(double *a, double *b, double *c, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            c[i * col + j] = a[i * col + j] - b[i * col + j];
        }
    }
}

void multiply(double *a, double *b, double *c, int row1, int col1, int row2, int col2) {
    if (col1 != row2) {
        printf("Matrix multiplication error: incompatible dimensions\n");
        return;
    }
    for (int i = 0; i < row1; i++) {
        for (int j = 0; j < col2; j++) {
            c[i * col2 + j] = 0;  // c: row1*col2
            for (int k = 0; k < col1; k++) {
                c[i * col2 + j] += a[i * col1 + k] * b[k * col2 + j];
            }
        }
    }
}

void transpose(double *a, double *b, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            b[j * row + i] = a[i * col + j];
        }
    }
}

vector3 cross(vector3 a, vector3 b) {
    vector3 result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

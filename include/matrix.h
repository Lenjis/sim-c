
typedef struct {
    double x;
    double y;
    double z;
} vector3;

void add(double *a, double *b, double *c, int row, int col);
void subtract(double *a, double *b, double *c, int row, int col);
void multiply(double *a, double *b, double *c, int row1, int col1, int row2, int col2);
void transpose(double *a, double *b, int row, int col);
vector3 cross(vector3 a, vector3 b);

#include "cn.h"

double ymw16_ne(double x_pc, double y_pc, double z_pc) {
    //printf("For ncrd=2, input gl, gb in deg, Dist in pc\n");
    int ncrd = 1;
    int vbs = 0;
    char dirname[256] = "data/";
    char text[64] = "";
    double gl, gb, dd;
    return ne_crd(&x_pc, &y_pc, &z_pc, &gl, &gb, &dd, ncrd, vbs, dirname, text);
}

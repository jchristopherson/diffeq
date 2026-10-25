#ifndef C_PLOT_HELPER_H_
#define C_PLOT_HELPER_H_

#ifdef __cplusplus
extern "C" {
#endif

void c_plot_xy(int n, const double *x, const double *y, const char *xlbl, const char *ylbl);
void c_plot_xy2(int n1, int n2, const double *x1, const double *y1, const char *name1,
    const double *x2, const double *y2, const char *name2, const char *xlbl, const char *ylbl);

#ifdef __cplusplus
}
#endif
#endif

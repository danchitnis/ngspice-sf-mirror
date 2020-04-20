#ifndef hicumL2_H
#define hicumL2_H
#include "hicum2defs.h"

#ifdef __cplusplus
extern "C" {
#endif
    void hicum_diode(double T, double IS, double UM1, double U, double *Iz, double *Gz, double *Tz);
    void hicum_qjmodf(double T, double c_0, double u_d, double z, double a_j, double U_cap, double *C, double *C_dU, double *C_dvt, double *Qz, double *Qz_dU, double *Qz_dvt);
    static double HICUMlimitlog( double deltemp, double deltemp_old, double LIM_TOL, int *check);
    int hicum_thermal_update(HICUMmodel *, HICUMinstance *);
    int HICUMload(GENmodel *inModel, CKTcircuit *ckt);
#ifdef __cplusplus
}
#endif

#endif /* hicumL2_H */
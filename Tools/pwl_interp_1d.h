#ifndef PWL_INTERP_1D_H
#define PWL_INTERP_1D_H

double *pwl_basis_1d ( int nd, double xd[], int ni, double xi[] );
double *pwl_value_1d ( int nd, double xd[], double yd[], int ni, double xi[] );
double pwl_value_1d_scalar ( int nd, double xd[], double yd[], double xi );

#endif
#ifndef FMM_SPLINE_H
#define FMM_SPLINE_H

/* Splines a la Forsythe Malcolm and Moler (from splines.c in the stats package).
 *
 * 'n' is the number of points, 'x' is the vector of x-coordinates, 'y' is the 
 * vector of y-coordinates (unchanged and equal to the constant for the interpolating
 * cubic spline upon output), 'b' is the coefficient degree 1, 'c' is the coefficient 
 * of degree '2', 'd' is the coefficient degree 3. We need  to access the coefficients 
 * directly, which makes evaluation from R undesirable.
 */

extern "C" {
void fmm_spline(int n, double *x, double *y, double *b, double *c, double *d);
}

#endif

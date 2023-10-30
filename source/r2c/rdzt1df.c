#include <math.h>
#include <complex.h>

/*
    A simple recursive implementation of 1-dimentional recursive Z-transform for 32-bit floats
*/

void rdzt1df(const float* x, const float* y, const int n, const int k, const float df, complex float *out, const float fmin){
/*
    x - time of measurements
    y - value of measurements
    n - number of measurements

    fmin - start frequency of transform (for pruned transforms, default = 0)
    df - frequency step size

    k - number of output frequencies
    out - complex output array
*/

    float *cosdx = malloc(n * sizeof(float)),
          *sindx = malloc(n * sizeof(float)),
          *cosx = malloc(n * sizeof(float)),
          *sinx = malloc(n * sizeof(float));

    float tmp, fk;

    //prepare trigonometric recurrences
    for (int i=0; i<n; i++) {
        cosdx[i] = cos(2 * M_PI * df * x[i]);
        sindx[i] = sin(2 * M_PI * df * x[i]);
    }

    for (int j=0; j<k; j++) {
        fk = fmin + (df * float(j));
        out[j] = CMPLX(0, 0); //CMPLX(real, imag);

        for (int i=0; i<n; i++) {

            if (j % 256 == 0) { /* refresh recurrences to stop error propagation */
                cosx[i] = cos(2 * M_PI * fk * x[i]);
                sinx[i] = sin(2 * M_PI * fk * x[i]);
                }

            out[j] += CMPLX(y[i] * cosx[i], y[i] * sinx[i]);

            tmp = (cosx[i] * cosdx[i]) - (sinx[i] * sindx[i]);
            sinx[i] = (cosx[i] * sindx[i]) + (sinx[i] * cosdx[i]);
            cosx[i] = tmp;
        }
    }
free(cosx); free(cosdx); free(sinx); free(sindx);
}

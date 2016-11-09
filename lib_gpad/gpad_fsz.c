/* 
 * File: gpad_fsz.c
 * Created on: 4 July, 2013
 * Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
 * Institute: IMT Lucca, Lucca Italy
 *    __  __      __  
 *   / _ |__) /\ |  \ 
 *   \__)|   /--\|__/  ~ CODE GENERATOR
 *
 */

#include "gpad_data.h"

#ifdef GPAD_DEBUG
#include <stdio.h>
#endif

#ifndef NULL
#define NULL 0
#endif

real_t square_root(creal_t x, creal_t epsilon, cus_t max_iter) {
    real_t a, a_p, diff;
    us_t i = 1, mi;
    if (max_iter == 0) {
        mi = 12;
    } else {
        mi = max_iter;
    }
    a = (real_t) x;
    a_p = a;
    do {
        a = 0.5 * (a_p + x / a_p);
        diff = a - a_p;
        if (diff < 0) diff = -diff;
        if (i > mi || diff / a < epsilon) break;
        a_p = a;
        i++;
    } while (1);
    return a;
}

/* y = y + a*x */
void zvadd(cus_t n, creal_t alpha, creal_t *x, real_t *y) {
    us_t i = 0;
    for (i = 0; i < n; i++) y[i] += alpha * x[i];
}

/* y = y + alpha*a*x + b */
void zmvmult(creal_t a[], creal_t x[], creal_t alpha, creal_t b[], cus_t nRows, cus_t nCols, real_t y[]) {
    us_t i, j;
    for (i = 0; i < nRows; i++) {
        for (j = 0; j < nCols; j++) y[i] += alpha * a[i + j * nRows] * x[j];
        if (b != NULL) y[i] += b[i];
    }
}

/* y = y + alpha*a'*x + b */
void zmvmult_trans(creal_t a[], creal_t x[], creal_t alpha, creal_t b[], cus_t nRows, cus_t nCols, real_t y[]) {
    us_t i, j;
    for (j = 0; j < nCols; j++) {
        for (i = 0; i < nRows; i++) y[j] += alpha * a[i + j * nRows] * x[i];
        if (b != 0) y[j] += b[j];
    }
}

/* returns: x'*q*x */
real_t zvdmquad(creal_t x[], cus_t nx, creal_t q[]) {
    real_t result;
    us_t i;
    for (i = 0; i < nx; i++) result += x[i] * x[i] * q[i];
    return result;
}

/* return alpha*x'*x */
real_t zvnrmsq(creal_t x[], cus_t n, creal_t alpha) {
    real_t norm = 0.0;
    us_t i;
    for (i = 0; i < n; i++) norm += alpha * x[i] * x[i];
    return norm;
}

real_t zcdot(creal_t x[], creal_t y[], cus_t nx, creal_t alpha) {
    real_t result = 0.0;
    us_t i;
    for (i = 0; i < nx; i++) result += alpha * x[i] * y[i];
    return result;
}

real_t zsym_quad(creal_t x[], creal_t s[], cus_t n) {
    us_t i, j;
    real_t result = 0.0;
    for (i = 0; i < n; i++) result += x[i] * x[i] * s[i * (n + 1)];
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) result += 2.0 * x[i] * x[j] * s[i + j * n];
    }
    return result;
}

/* z = x'*Q'*y */
real_t zquad(creal_t x[], creal_t q[], creal_t y[], cus_t nx, cus_t ny, creal_t alpha) {
    us_t i, j;
    real_t result = 0.0;
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) result += alpha * x[i] * y[j] * q[j + i * nx];
    }
    return result;
}

void zlowslv(creal_t low_triang[], real_t x[], cus_t n) {
    us_t i, j, offset;
    real_t sum = 0.0;
    x[0] = x[0] / low_triang[0];
    for (i = 1; i < n; i++) {
        offset = ((i + 1) * i) / 2;
        sum = x[i];
        for (j = 0; j < i; j++) {
            sum -= (x[j] * low_triang[offset + j]);
        }
        x[i] = sum / low_triang[offset + i];
    }
}

void zlowtrslv(creal_t low_triang[], real_t x[], cus_t n) {
    us_t i, j, s;
    cus_t size_low = (n * (n + 1)) / 2;
    real_t sum = 0.0;
    x[n - 1] = x[n - 1] / low_triang[size_low - 1];
    for (i = 1; i < n; i++) {
        sum = x[n - i - 1];
        s = size_low - i - 1;
        for (j = 0; j < i; j++) {
            sum -= (low_triang[s] * x[n - 1 - j]);
            s -= (n - j - 1);
        }
        x[n - i - 1] = sum / low_triang[s];
    }
}

void zcholslv(creal_t low_triang[], real_t x[], cus_t n) {
    zlowslv(low_triang, x, n);
    zlowtrslv(low_triang, x, n);
}


void control(creal_t x0[]) {
    us_t i, doloop = TRUE, iter = 0;
    gBETA = gTHETA * (1 / gTHETA_P - 1);
    while (doloop == TRUE && iter <= MAX_ITER_) {
        for (i = 0; i < NP_N_; i++) gW[i] = gY[i] + gBETA * (gY[i] - gY_PREV[i]);
        for (i = 0; i < NF_; i++) gWN[i] = gYN[i] + gBETA * (gYN[i] - gYN_PREV[i]);
        for (i = 1; i < NP_N_; i++) gY_PREV[i] = gY[i];
        for (i = 1; i < NF_; i++) gYN_PREV[i] = gYN[i];
    }
}

void calculate_e(void) {
    us_t i, k;
    for (i = 0; i < NX_; i++) gE[i] = 0.0;
    zmvmult_trans(FN_, gWN, 1.0, NULL, NF_, NX_, gE + NX_ * (N_ - 1));
#if (defined(IS_EMPTY_qN) && !IS_EMPTY_qN)
    for (i = 0; i < NX_; i++) gE[i] += Q_SMALL_N_[i];
#endif
    for (k = N_; k > 1; k--) { /* THE FOR LOOP */
        zmvmult(L_, gE + NX_*k, 1.0, NULL, NX_, NX_, gE + NX_ * (k - 1));
        i = (k - 1) * NP_;
#if !IS_EMPTY_CMAX
        zmvmult(C_, gW + i, 1.0, NULL, NX_, NC_, gE + NX_ * (k - 1));
        i += NC_;
#endif
#if !IS_EMPTY_CMIN
        zmvmult(C_, gW + i, -1.0, NULL, NX_, NC_, gE + NX_ * (k - 1));
        i += NC_;
#endif
#if !IS_EMPTY_XMAX
        zvadd(NX_, 1.0, gW + i, gE + NX_ * (k - 1));
        i += NX_;
#endif
#if !IS_EMPTY_XMIN
        zvadd(NX_, -1.0, gW + i, gE + NX_ * (k - 1));
        i += NX_;
#endif 
#if !IS_EMPTY_UMAX
        zmvmult_trans(K_, gW + i, 1.0, NULL, NU_, NX_, gE + NX_ * (k - 1));
        i += NU_;
#endif
#if !IS_EMPTY_UMIN
        zmvmult_trans(K_, gW + i, -1.0, NULL, NU_, NX_, gE + NX_ * (k - 1));
#endif
    }
}

void primal_feasibility(creal_t x[], creal_t u[]) {
    /* Declarations */
    us_t i = 0, j = 0;

#if ((!IS_EMPTY_CMAX) || (!IS_EMPTY_CMIN))
    zmvmult(F_, x, 1.0, NULL, NC_, NX_, gFXGU);
    zmvmult(G_, u, 1.0, NULL, NC_, NU_, gFXGU);
#if(!IS_EMPTY_CMAX)
    for (j = 0; j < NC_; j++) {
        gPR_FEAS[j] = gFXGU[j] - CMAX_[j];
    }
    i += NC_;
#endif /* #if(!IS_EMPTY_CMAX) */
#if (!IS_EMPTY_CMIN)    
    for (j = 0; j < NC_; j++) {
        gPR_FEAS[j + i] = CMIN_[j] - gFXGU[j];
    }
    i += NC_;
#endif /* #if (!IS_EMPTY_CMIN) */
#endif /* #if ((!IS_EMPTY_CMAX) || (!IS_EMPTY_CMIN)) */

#if (!IS_EMPTY_XMAX)
    for (j = 0; j < NX_; j++) {
        gPR_FEAS[j + i] = x[j] - XMAX_[j];
    }
    i += NX_;
#endif
#if (!IS_EMPTY_XMIN)
    for (j = 0; j < NX_; j++) {
        gPR_FEAS[j + i] = XMIN_[j] - x[j];
    }
    i += NX_;
#endif    
#if (!IS_EMPTY_UMAX)
    for (j = 0; j < NU_; j++) {
        gPR_FEAS[j + i] = u[j] - UMAX_[j];
    }
    i += NU_;
#endif
#if (!IS_EMPTY_UMIN)
    for (j = 0; j < NU_; j++) {
        gPR_FEAS[j + i] = UMIN_[j] - u[j];
    }
#endif

}

int main(void) {

    us_t i;
    for (i = 0; i < NF_; i++) gWN[i] = (real_t) i / 10.0;
    calculate_e();
#if defined(_STDIO_H_) && defined(GPAD_DEBUG)
    for (i = 0; i < NX_N_; i++) printf("e[%d]=%3.8f\n", i, gE[i]);
#endif
    return 2;
}


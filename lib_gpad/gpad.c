/* 
 * File: gpad.c
 * Created on: 15 June, 2013
 * Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
 * Institute: IMT Lucca, Lucca Italy
 *    __  __      __  
 *   / _ |__) /\ |  \ 
 *   \__)|   /--\|__/ 
 */

#include "gpad.h"

/* max */
__INLINE__ real_t max(creal_t *array, cus_t n) {
    us_t i;
    real_t max_val = -VERY_BIG_DBL;
    for (i = 0; i < n; i++) {
        if (array[i] >= max_val) max_val = array[i];
    }
    return max_val;
}

/* min */
__INLINE__ real_t min(creal_t *array, cus_t n) {
    us_t i;
    real_t min_val = VERY_BIG_DBL;
    for (i = 0; i < n; i++) {
        if (array[i] <= min_val) min_val = array[i];
    }
    return min_val;
}

/* handle_error */
static __INLINE__ error_t handle_error(const error_t err, const char *funct_name) {
#if (defined(GPAD_DEBUG))
    if (err != SUCCESS_OK) {
        char error_message[70];
        if (err != SUCCESS_OK) {
            sprintf(error_message, "ERROR %s - error %d\n", funct_name, err);
            printf(error_message);
        }
    }
#endif
    return err;
}

/* square_root */
__INLINE__ real_t square_root(creal_t x, creal_t epsilon, cus_t *max_iter) {
    real_t a, a_p, diff;
    us_t i = 1, mi;
    if (max_iter == NULL) {
        mi = 12;
    } else {
        mi = *max_iter;
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

/* dynsys_factory */
dynsys_t *dynsys_factory(creal_t *A, creal_t *B, creal_t *f, creal_t *F, creal_t *G,
        creal_t *cmin, creal_t *cmax, creal_t *xmin, creal_t *xmax, creal_t *umin,
        creal_t *umax, cus_t nx, cus_t nu, cus_t nc) {
    dynsys_t *p = NULL;
    if ((p = (dynsys_t*) malloc(sizeof (dynsys_t))) != NULL) {
        p->A = (real_t*) A;
        p->B = (real_t*) B;
        p->f = (real_t*) f;
        p->F = (real_t*) F;
        p->G = (real_t*) G;
        p->cmin = (real_t*) cmin;
        p->cmax = (real_t*) cmax;
        p->xmin = (real_t*) xmin;
        p->xmax = (real_t*) xmax;
        p->umin = (real_t*) umin;
        p->umax = (real_t*) umax;
        p->nx = (us_t) nx;
        p->nu = (us_t) nu;
        p->nc = (us_t) nc;
    } else {
        return NULL;
    }
    return p;
}

/* dynsys_factory_meagre */
dynsys_t *dynsys_factory_meagre(cus_t nx, cus_t nu, cus_t nc) {
    dynsys_t *p = NULL;
    if ((p = (dynsys_t*) malloc(sizeof (dynsys_t))) != NULL) {
        p->nx = (us_t) nx;
        p->nu = (us_t) nu;
        p->nc = (us_t) nc;
    } else {
        return NULL;
    }
    return p;
}

/* mpc_factory */
mpc_t *mpc_factory(creal_t *Q, creal_t *R, creal_t *S, creal_t *q, creal_t *r, creal_t *QN,
        creal_t *qN, creal_t *K, creal_t *M, creal_t *D, creal_t *L, creal_t *C, creal_t *s, creal_t *FN,
        creal_t *gN, cus_t classQ, cus_t classR, cus_t N, cus_t nf, const dynsys_t *sys) {
    mpc_t *mpc = NULL;
    if ((mpc = (mpc_t*) malloc(sizeof (mpc_t))) != NULL) {
        mpc->Q = (real_t*) Q;
        mpc->R = (real_t*) R;
        mpc->S = (real_t*) S;
        mpc->q = (real_t*) q;
        mpc->r = (real_t*) r;
        mpc->QN = (real_t*) QN;
        mpc->qN = (real_t*) qN;
        mpc->K = (real_t*) K;
        mpc->M = (real_t*) M;
        mpc->D = (real_t*) D;
        mpc->L = (real_t*) L;
        mpc->C = (real_t*) C;
        mpc->s = (real_t*) s;
        mpc->FN = (real_t*) FN;
        mpc->gN = (real_t*) gN;
        mpc->classQ = (us_t) classQ;
        mpc->classR = (us_t) classR;
        mpc->N = (us_t) N;
        mpc->nf = (us_t) nf;
        mpc->sys = (dynsys_t*) sys;
        return mpc;
    } else {
        return NULL;
    }
}

/* mpc_factory_simple */
mpc_t *mpc_factory_simple(creal_t *Q, creal_t *R, creal_t *S, creal_t *q, creal_t *r, creal_t *QN,
        creal_t *qN, cus_t classQ, cus_t classR, cus_t N, cus_t nf, const dynsys_t *sys) {
    mpc_t *mpc = NULL;
    if ((mpc = (mpc_t*) malloc(sizeof (mpc_t))) != NULL) {
        mpc->Q = (real_t*) Q;
        mpc->R = (real_t*) R;
        mpc->S = (real_t*) S;
        mpc->q = (real_t*) q;
        mpc->r = (real_t*) r;
        mpc->QN = (real_t*) QN;
        mpc->qN = (real_t*) qN;
        mpc->classQ = (us_t) classQ;
        mpc->classR = (us_t) classR;
        mpc->N = (us_t) N;
        mpc->nf = (us_t) nf;
        mpc->sys = (dynsys_t*) sys;
        return mpc;
    } else {
        return NULL;
    }
}

/* state_update */
error_t state_update(creal_t *A, creal_t *B, creal_t *f, creal_t *x,
        creal_t *u, cus_t nx, cus_t nu, real_t *x_new) {
    error_t error_zmvmult;
    us_t doAddOver = DO_NOT_ADD_OVER;
    if (A == NULL || B == NULL || x == NULL || u == NULL) return NULL_INPUT;
    if (x_new == NULL) return NULL_OUTPUT;
    if (f != NULL) {
        memcpy(x_new, f, nx * DBL_SIZE);
        doAddOver = DO_ADD_OVER;
    }
    error_zmvmult = zmvmult(A, x, NULL, NULL, nx, nx, doAddOver, x_new);
    if (error_zmvmult != SUCCESS_OK) return error_zmvmult;
    error_zmvmult = zmvmult(B, u, NULL, NULL, nx, nu, DO_ADD_OVER, x_new);
    if (error_zmvmult != SUCCESS_OK) return error_zmvmult;
    return SUCCESS_OK;
}

/* state_update */
error_t state_update_i(dynsys_t *sys, creal_t *x, creal_t *u, real_t *x_new) {
    error_t error_zmvmult;
    us_t doAddOver = DO_NOT_ADD_OVER;
    if (sys == NULL || x == NULL || u == NULL) return NULL_INPUT;
    if (sys->A == NULL || sys->B == NULL) return NULL_INPUT;
    if (x_new == NULL) return NULL_OUTPUT;
    if (sys->f != NULL) {
        memcpy(x_new, sys->f, (sys->nx) * DBL_SIZE);
        doAddOver = DO_ADD_OVER;
    }
    error_zmvmult = zmvmult(sys->A, x, NULL, NULL, sys->nx, sys->nx, doAddOver, x_new);
    if (error_zmvmult != SUCCESS_OK) return error_zmvmult;
    error_zmvmult = zmvmult(sys->B, u, NULL, NULL, sys->nx, sys->nu, DO_ADD_OVER, x_new);
    if (error_zmvmult != SUCCESS_OK) return error_zmvmult;
    return SUCCESS_OK;
}

/* sym_quad */
error_t zsym_quad(creal_t *x, creal_t *S, cus_t n, cus_t doAddOver, real_t *L) {
    us_t i, j;
    if (L == NULL) return NULL_OUTPUT;
    if (x == NULL || S == NULL) return NULL_INPUT;
    if (doAddOver == DO_NOT_ADD_OVER) *L = 0.0;
    for (i = 0; i < n; i++) (*L) += x[i] * x[i] * S[i * (n + 1)];
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) (*L) += 2.0 * x[i] * x[j] * S[i + j * n];
    }
    return SUCCESS_OK;
}

/* zvnrmsq */
error_t zvnrmsq(creal_t *x, cus_t n, creal_t *alpha,
        cus_t doAddOver, real_t *norm) {
    us_t i;
    real_t *ytemp;
    if (doAddOver == DO_NOT_ADD_OVER) *norm = 0.0;
    if (norm == NULL) return NULL_OUTPUT;
    ytemp = calloc(1, DBL_SIZE);
    if (ytemp == NULL) return ALLOCATION_ERROR;
    for (i = 0; i < n; i++) *ytemp += x[i] * x[i];
    if (alpha != 0) *ytemp *= *alpha;
    *norm += *ytemp;
    if (ytemp != NULL) free(ytemp);
    return SUCCESS_OK;
}

/* zvdmquad: y += alpha*x'*Q*x */
error_t zvdmquad(creal_t *x, cus_t nx, creal_t *Q,
        creal_t *alpha, cus_t doAddOver, real_t *y) {
    real_t *ytemp;
    us_t i;
    if (y == NULL) return NULL_OUTPUT;
    if (x == NULL || Q == NULL) return NULL_INPUT;
    if (doAddOver == DO_NOT_ADD_OVER) *y = 0.0;
    ytemp = calloc(1, DBL_SIZE);
    if (ytemp == NULL) return ALLOCATION_ERROR;
    for (i = 0; i < nx; i++) *ytemp += x[i] * x[i] * Q[i];
    if (alpha != 0) *ytemp *= *alpha;
    *y += *ytemp;
    if (ytemp != NULL) free(ytemp);
    return SUCCESS_OK;
}

/* zquad :  z = x'*Q'*y */
error_t zquad(creal_t *x, creal_t *Q, creal_t *y, cus_t nx,
        cus_t ny, creal_t *alpha, cus_t doAddOver, real_t *z) {
    /* z=alpha*x'*Q'*y */
    us_t i, j;
    real_t *ztemp = NULL;
    if (x == NULL) return handle_error(NULL_INPUT, "x==NULL/zquad");
    if (Q == NULL) return handle_error(NULL_INPUT, "Q==NULL/zquad");
    if (z == 0) return handle_error(NULL_OUTPUT, "z==NULL/zquad");
    if (doAddOver == DO_NOT_ADD_OVER) *z = 0.0;
    ztemp = calloc(1, DBL_SIZE);
    if (ztemp == NULL) return ALLOCATION_ERROR;
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) *ztemp += x[i] * y[j] * Q[j + i * nx];
    }
    if (alpha != 0) *ztemp *= *alpha;
    *z += *ztemp;
    if (ztemp != NULL) free(ztemp);
    return SUCCESS_OK;
}

/* zmvmult:  y = alpha*A*x+b */
error_t zmvmult(creal_t *A, creal_t *x, creal_t *alpha, creal_t *b, cus_t nRows,
        cus_t nCols, cus_t doAddOver, real_t *y) {
#if defined(BLAS)
    double a;
#else
    us_t i, j; /* i: Row Index,   j; Column Index */
    real_t *ytemp = NULL;
#endif   
    if (A == NULL || x == NULL) return handle_error(NULL_INPUT, "zmvmult");
    if (y == NULL) return NULL_OUTPUT;
#if defined(BLAS) /* Matrix-Vector Multiplication with BLAS */
    a = 1.0;
    if (alpha != NULL) a = *alpha;
    cblas_dgemv(CblasColMajor, CblasNoTrans, nRows, nCols, a, A, nRows, x, 1, (double) doAddOver, y, 1);
    if (b != NULL) cblas_daxpy(nRows, 1.0, b, 1, y, 1);
#else /* Do not use BLAS */    
    ytemp = calloc(nRows, DBL_SIZE);
    if (ytemp == NULL) return ALLOCATION_ERROR;
    for (i = 0; i < nRows; i++) {
        for (j = 0; j < nCols; j++) ytemp[i] += A[i + j * nRows] * x[j];
        if (alpha != NULL) ytemp[i] *= *alpha;
        if (b != NULL) ytemp[i] += b[i];
        if (doAddOver == DO_ADD_OVER) ytemp[i] += y[i];
    }
    memcpy(y, ytemp, nRows * DBL_SIZE);
    if (ytemp != NULL) free(ytemp);
#endif /* end defined(BLAS) */
    return SUCCESS_OK;
}

/* zmvmult_trans:  y = alpha*A'*x+b */
error_t zmvmult_trans(creal_t *A, creal_t *x, creal_t *alpha, creal_t *b, cus_t nRows,
        cus_t nCols, cus_t doAddOver, real_t *y) {
#if defined(BLAS)
    double a;
#else
    us_t i, j; /* i: Row Index,   j; Column Index */
    real_t *ytemp;
#endif
    if (A == NULL || x == NULL) return NULL_INPUT;
    if (y == NULL) return NULL_OUTPUT;
#if defined(BLAS) /* Matrix-Vector Multiplication with BLAS */
    a = 1.0;
    if (alpha != NULL) a = *alpha;
    cblas_dgemv(CblasColMajor, CblasTrans, nRows, nCols, a, A, nRows, x, 1, (double) doAddOver, y, 1);
    if (b != NULL) cblas_daxpy(nCols, 1.0, b, 1, y, 1);
#else                
    ytemp = calloc(nCols, DBL_SIZE);
    if (ytemp == NULL) return ALLOCATION_ERROR;
    for (j = 0; j < nCols; j++) {
        for (i = 0; i < nRows; i++) ytemp[j] += A[i + j * nRows] * x[i];
        if (alpha != 0) ytemp[j] *= *alpha;
        if (b != 0) ytemp[j] += b[j];
        if (doAddOver == DO_ADD_OVER) ytemp[j] += y[j];
    }
    memcpy(y, ytemp, nCols * DBL_SIZE);
    if (ytemp != NULL) free(ytemp);
#endif
    return SUCCESS_OK;
}

/* zcdot: z = alpha*x'*y */
error_t zcdot(creal_t *x, creal_t *y, cus_t nx, creal_t *alpha,
        cus_t doAddOver, real_t *z) {
    /* Declarations */
#if defined(BLAS) /* Dot product with BLAS*/
    double dot_prod;
#else
    us_t i;
    real_t *ztemp;
#endif
    if (x == NULL || y == NULL) return handle_error(NULL_INPUT, "zcdot");
    if (z == NULL) return handle_error(NULL_OUTPUT, "z==NULL/zcdot");
#if defined(BLAS) /* Dot product with BLAS*/
    dot_prod = cblas_ddot(nx, x, 1, y, 1);
    if (doAddOver == DO_ADD_OVER) *z += dot_prod;
    if (doAddOver == DO_NOT_ADD_OVER) *z = dot_prod;
#else /* Dot product without BLAS */       
    if (doAddOver == DO_NOT_ADD_OVER) *z = 0.0;
    ztemp = calloc(1, DBL_SIZE);
    if (ztemp == NULL) return ALLOCATION_ERROR;
    for (i = 0; i < nx; i++) *ztemp += x[i] * y[i];
    if (alpha != 0) *ztemp *= *alpha;
    *z += *ztemp;
#endif /* end defined(BLAS) */
    return SUCCESS_OK;
}

/* zvadd: Do y = y + ax */
error_t zvadd(cus_t N, creal_t alpha, creal_t *X, real_t *Y) {
#if !defined(BLAS)
    us_t i;
#endif
    if (X == NULL) return NULL_INPUT;
    if (Y == NULL) return NULL_OUTPUT;
#if defined(BLAS)
    cblas_daxpy(N, alpha, X, 1, Y, 1);
#else
    for (i = 0; i < N; i++) Y[i] = Y[i] + alpha * X[i];
#endif
    return SUCCESS_OK;
}

/* zlowslv */
error_t zlowslv(creal_t *L, real_t *x, cus_t n) {
#if !defined(BLAS)
    us_t i, j;
    real_t sum;
#endif
    if (L == NULL) return NULL_OUTPUT;
    if (x == NULL) return NULL_INPUT;
#if defined(BLAS)    
    cblas_dtrsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, L, n, x, 1);
#else
    if (*L == 0) return handle_error(DEVISION_BY_ZERO, "zlowslv/dbz$781");
    *x = *x / *L;
    for (i = 1; i < n; i++) {
        sum = x[i];
        for (j = 0; j < i; j++) {
            sum -= (x[j] * L[j * n + i]);
        }
        if (L[i * (n + 1)] == 0) return handle_error(DEVISION_BY_ZERO, "zlowslv/dbz$783");
        x[i] = sum / L[i * (n + 1)];
    }
#endif
    return SUCCESS_OK;
}

#undef BLAS

/* zlowtrslv */
error_t zlowtrslv(creal_t *L, real_t *x, cus_t n) {
#if !defined(BLAS)
    us_t i, j;
    real_t sum;
#endif
#if defined(BLAS)    
    cblas_dtrsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit, n, L, n, x, 1);
#else    
    if (L[n * n - 1] == 0) return handle_error(DEVISION_BY_ZERO, "zlowtrslv/dbz$791");
    x[n - 1] = x[n - 1] / L[n * n - 1];
    for (i = 1; i < n; i++) {
        sum = x[n - 1 - i];
        for (j = 0; j < i; j++) {
            sum -= (L[n * (n - 1 - i) + (n - i + j)] * x[n - i + j]);
        }
        x[n - 1 - i] = sum / L[(n - 1 - i)*(1 + n)];
    }

#endif
    return SUCCESS_OK;
}

/* zchlslv: solve b<--x, x<--(LL')\b*/
error_t zchlslv(
        creal_t *L,
        real_t *x,
        cus_t n
        ) {
    error_t err;
    err = zlowslv(L, x, n);
    if (err != SUCCESS_OK) return handle_error(err, "zchlslv/zlowslv");
    err = zlowtrslv(L, x, n);
    if (err != SUCCESS_OK) return handle_error(err, "zchlslv/zlowtrslv");
    return SUCCESS_OK;
}

/* zchlslv_cpy: y = (LL')\x */
error_t zchlslv_cpy(
        creal_t *L,
        creal_t *x,
        real_t *y,
        cus_t n
        ) {
    error_t err;
    if (L == NULL || x == NULL) return handle_error(NULL_INPUT, "zchlslv_cpy$1");
    if (y == NULL) return handle_error(NULL_OUTPUT, "zchlslv_cpy$2");
    memcpy(y, x, n * DBL_SIZE);
    err = zchlslv(L, y, n);
    if (err != SUCCESS_OK) return handle_error(err, "zchlslv_cpy/zchlslv");
    return SUCCESS_OK;
}

/* pcost_i */
error_t pcost_i(
        const mpc_t *mpc,
        creal_t *X,
        creal_t *U,
        real_t *J
        ) {
    if (mpc == NULL || X == NULL || U == NULL) return NULL_INPUT;
    if (J == NULL) return NULL_OUTPUT;
    return pcost(mpc->Q, mpc->R, mpc->S, mpc->q, mpc->r, mpc->QN, mpc->qN,
            mpc->N, X, U, mpc->sys->nx, mpc->sys->nu, mpc->classQ, mpc->classR, J);
}

/* pcost */
error_t pcost(
        creal_t *Q,
        creal_t *R,
        creal_t *S,
        creal_t *q,
        creal_t *r,
        creal_t *QN,
        creal_t *qN,
        cus_t N,
        creal_t *X,
        creal_t *U,
        cus_t nx,
        cus_t nu,
        cus_t classQ,
        cus_t classR,
        real_t *J
        ) {
    creal_t *currentState = NULL, *currentInput = NULL;
    us_t k;
    real_t *cost_temp = NULL;
    if (Q == NULL || R == NULL || X == NULL || U == NULL) return NULL_INPUT;
    if (N == 0) return handle_error(BAD_PARAM_VALUE, "pcost (N==0)");
    if (J == NULL) return NULL_OUTPUT;
    currentState = X;
    currentInput = U;

    cost_temp = malloc(DBL_SIZE);
    if (cost_temp == NULL) return handle_error(ALLOCATION_ERROR, "pcost(cost_temp)");
    *cost_temp = 0.0;


    /* For the cost of all stages: */
    /* __builtin_prefetch(J); */
    for (k = 0; k < N; k++) { /* Terms x'*Q*x ... */
        if (classQ == MAT_aI) {
            zvnrmsq(currentState, nx, Q, DO_ADD_OVER, cost_temp);
        } else if (classQ == MAT_DIAG) {
            zvdmquad(currentState, nx, Q, NULL, DO_ADD_OVER, cost_temp);
        } else if (classQ == MAT_DENSE) {
            zquad(currentState, Q, currentState, nx, nx, NULL, DO_ADD_OVER, cost_temp);
        } else {
            return NOT_IMPLEMENTED_YET;
        }
        *J += (*cost_temp / 2.0);
        *cost_temp = 0.0;
        /* Terms u'*R*u ... */
        if (classR == MAT_aI) {
            zvnrmsq(currentInput, nu, R, DO_ADD_OVER, cost_temp);
        } else if (classR == MAT_DIAG) {
            zvdmquad(currentInput, nu, R, NULL, DO_ADD_OVER, cost_temp);
        } else if (classR == MAT_DENSE) {
            zquad(currentInput, R, currentInput, nu, nu, NULL, DO_ADD_OVER, cost_temp);
        } else {
            return NOT_IMPLEMENTED_YET;
        }
        *J += (*cost_temp / 2.0);
        *cost_temp = 0.0;
        if (S != NULL) zquad(currentState, S, currentInput, nu, nx, NULL, DO_ADD_OVER, J);
        if (q != NULL) zcdot(q, currentState, nx, NULL, DO_ADD_OVER, J);
        if (r != NULL) zcdot(r, currentInput, nu, NULL, DO_ADD_OVER, J);
        currentState += nx;
        currentInput += nu;
    } /* end: stage cost k=0,...,N-1 */
    /* For the terminal cost: */
    if (QN != NULL) zsym_quad(currentState, QN, nx, DO_ADD_OVER, cost_temp); /* QN: is symmetric */
    *J += (*cost_temp / 2.0);
    *cost_temp = 0.0;
    if (qN != NULL) zcdot(currentState, qN, nx, NULL, DO_ADD_OVER, J);
    if (cost_temp != NULL) free(cost_temp);
    return SUCCESS_OK;
}

/* dcost */
error_t dcost(
        creal_t *A,
        creal_t *B,
        creal_t *f,
        creal_t *Q,
        creal_t *R,
        creal_t *S,
        creal_t *q,
        creal_t *r,
        creal_t *QN,
        creal_t *qN,
        creal_t *FN,
        creal_t *gN,
        creal_t *F,
        creal_t *G,
        creal_t *cmin,
        creal_t *cmax,
        creal_t *xmin,
        creal_t *xmax,
        creal_t *umin,
        creal_t *umax,
        creal_t *D,
        creal_t *L,
        creal_t *M,
        creal_t *K,
        creal_t *C,
        creal_t *s,
        creal_t *d,
        creal_t *Rchol,
        creal_t *y,
        creal_t *yN,
        creal_t *x0,
        cus_t N,
        cus_t nx,
        cus_t nu,
        cus_t nc,
        cus_t nf,
        cus_t classQ,
        cus_t classR,
        real_t *Xstar,
        real_t *Ustar,
        real_t *dual_cost
        ) {

    /* Declarations */
    error_t err;
    us_t k, i, np = 0;
    real_t *e = NULL, *minus_one = NULL, *ytemp = NULL, *slack = NULL, *slackN = NULL;
    real_t *cost_temp = NULL;
    cus_t isempty_cmin = (cmin == NULL), isempty_cmax = (cmax == NULL), isempty_xmin = (xmin == NULL),
            isempty_xmax = (xmax == NULL), isempty_umin = (umin == NULL), isempty_umax = (umax == NULL);
    /* End of Declarations */

    /* Check input variables: */
    if (dual_cost == NULL) return handle_error(NULL_OUTPUT, "dual_cost==NULL/dcost");
    if (A == NULL) return handle_error(NULL_INPUT, "A==NULL/dcost");
    if (B == NULL) return handle_error(NULL_INPUT, "B==NULL/dcost");
    if (K == NULL) return handle_error(NULL_INPUT, "K==NULL/dcost");
    if (L == NULL) return handle_error(NULL_INPUT, "L==NULL/dcost");
    if (M == NULL) return handle_error(NULL_INPUT, "M==NULL/dcost");
    if (FN == NULL) return handle_error(NULL_INPUT, "FN==NULL/dcost");
    if (gN == NULL) return handle_error(NULL_INPUT, "gN==NULL/dcost");
    if (Q == NULL) return handle_error(NULL_INPUT, "Q==NULL/dcost");
    if (R == NULL) return handle_error(NULL_INPUT, "R==NULL/dcost");
    if (y == NULL) return handle_error(NULL_INPUT, "y==NULL/dcost");
    if (yN == NULL) return handle_error(NULL_INPUT, "yN==NULL/dcost");
    if (x0 == NULL) return handle_error(NULL_INPUT, "x0==NULL/dcost");
    if (Rchol == NULL) return handle_error(NULL_INPUT, "Rchol==NULL/dcost");
    if (Xstar == NULL) return handle_error(NULL_INPUT, "Xstar==NULL/dcost");
    if (Ustar == NULL) return handle_error(NULL_INPUT, "Ustar==NULL/dcost");

    cost_temp = malloc(DBL_SIZE);
    if (cost_temp == NULL) return handle_error(ALLOCATION_ERROR, "pcost(cost_temp)");


    np = nc * (2 - isempty_cmax - isempty_cmin) + nx * (2 - isempty_xmax - isempty_xmin) +
            nu * (2 - isempty_umax - isempty_umin);
    if (np == 0) return handle_error(BAD_PARAM_VALUE, "np==0 (No constraints given)");

    slack = malloc(np * DBL_SIZE);
    if (slack == NULL) return handle_error(ALLOCATION_ERROR, "dcost/malloc:slack");
    slackN = malloc(nf * DBL_SIZE);
    if (slackN == NULL) return handle_error(ALLOCATION_ERROR, "dcost/malloc:slackN");

    /* CALCULATE THE MATRIX e : */
    e = malloc(nx * (N + 1) * DBL_SIZE);
    err = calculate_e(FN, qN, K, L, C, s, y, yN, isempty_cmin, isempty_cmax,
            isempty_xmin, isempty_xmax, isempty_umin, isempty_umax, nx, nu, nc, nf, N, e);
    if (err != SUCCESS_OK) return handle_error(err, "calculate_e/dcost");

    /* Solution Step: */
    memcpy(Xstar, x0, nx * DBL_SIZE);
    *dual_cost = 0.0;
    minus_one = malloc(DBL_SIZE);
    if (minus_one == NULL) return handle_error(ALLOCATION_ERROR, "dcost/minus_one");
    *minus_one = -1.0;

    ytemp = malloc(nu * DBL_SIZE);
    if (ytemp == NULL) return handle_error(ALLOCATION_ERROR, "dcost/malloc ytemp");


    for (k = 0; k < N; k++) {
        *cost_temp = 0.0;
        err = zmvmult(K, Xstar + k*nx, NULL, NULL, nu, nx, DO_NOT_ADD_OVER, Ustar + k * nu); /* u = Kx */
        if (err != SUCCESS_OK) return handle_error(err, "dcost/zmvmult");
        err = zmvmult(M, e + (k + 1) * nx, NULL, NULL, nu, nx, DO_ADD_OVER, Ustar + k * nu); /* u += M*e(k+1)*/
        if (err != SUCCESS_OK) return handle_error(err, "dcost/zmvmult");
        i = 0;
        if (isempty_cmax == FALSE) {
            err = zmvmult(D, y + k * np + i, NULL, NULL, nu, nc, DO_ADD_OVER, Ustar + k * nu);
            if (err != SUCCESS_OK) return handle_error(err, "dcost/zmvmult");
            i += nc;
        }
        if (isempty_cmin == FALSE) {
            err = zmvmult(D, y + k * np + i, minus_one, NULL, nu, nc, DO_ADD_OVER, Ustar + k * nu);
            if (err != SUCCESS_OK) return handle_error(err, "dcost/zmvmult");
            i += nc;
        }
        i += nx * (2 - isempty_xmax - isempty_xmax);
        if (isempty_umax == FALSE) {
            err = zchlslv_cpy(Rchol, y + k * np + i, ytemp, nu);
            if (err != SUCCESS_OK) return handle_error(err, "dcost/zchlslv_cpy");
            if (zvadd(nu, -1.0, ytemp, Ustar + k * nu) != SUCCESS_OK) return handle_error(err, "dcost/zvadd");
            i += nu;
        }
        if (isempty_umin == FALSE) {
            zchlslv_cpy(Rchol, y + k * np + i, ytemp, nu);
            if (zvadd(nu, 1.0, ytemp, Ustar + k * nu) != SUCCESS_OK) return handle_error(err, "dcost/zvadd");
        }
        /* Term x'Qx/2*/
        if (classQ == MAT_aI) {
            zvnrmsq(Xstar + k*nx, nx, Q, DO_ADD_OVER, cost_temp);
        } else if (classQ == MAT_DIAG) {
            zvdmquad(Xstar + k*nx, nx, Q, NULL, DO_ADD_OVER, cost_temp);
        } else if (classQ == MAT_DENSE) {
            zquad(Xstar + k*nx, Q, Xstar + k*nx, nx, nx, NULL, DO_ADD_OVER, cost_temp);
        } else {
            return handle_error(NOT_IMPLEMENTED_YET, "Sparse Matrix - NIY");
        }
        *dual_cost += (*cost_temp / 2.0);
        *cost_temp = 0.0;
        /* Terms u'*R*u/2 ... */
        if (classR == MAT_aI) {
            zvnrmsq(Ustar + k*nu, nu, R, DO_ADD_OVER, cost_temp);
        } else if (classR == MAT_DIAG) {
            zvdmquad(Ustar + k*nu, nu, R, NULL, DO_ADD_OVER, cost_temp);
        } else if (classR == MAT_DENSE) {
            zquad(Ustar + k*nu, R, Ustar + k*nu, nu, nu, NULL, DO_ADD_OVER, cost_temp);
        } else {
            return handle_error(NOT_IMPLEMENTED_YET, "Sparse Matrix - NIY");
        }
        *dual_cost += (*cost_temp / 2.0);
        *cost_temp = 0.0;

        if (S != NULL) zquad(Xstar + k * nx, S, Ustar + k * nu, nu, nx, NULL, DO_ADD_OVER, dual_cost);
        if (q != NULL) zcdot(q, Xstar + k * nx, nx, NULL, DO_ADD_OVER, dual_cost);
        if (r != NULL) zcdot(r, Ustar + k * nu, nu, NULL, DO_ADD_OVER, dual_cost);
        primal_feasibility(F, G, cmin, cmax, umin, umax, xmin, xmax, Xstar + k*nx, Ustar + k*nu, nx, nu, nc, np, slack);
        zcdot(slack, y + k*np, np, NULL, DO_ADD_OVER, dual_cost);
        state_update(A, B, f, Xstar + k*nx, Ustar + k*nu, nx, nu, Xstar + (k + 1) * nx);
    }
    zquad(Xstar + N*nx, QN, Xstar + N*nx, nx, nx, NULL, DO_ADD_OVER, cost_temp);
    *dual_cost += (*cost_temp / 2.0);

    zmvmult(FN, Xstar + N*nx, NULL, NULL, nf, nx, DO_NOT_ADD_OVER, slackN);
    for (i = 0; i < nf; i++) {
        slackN[i] -= gN[i];
    }
    zcdot(slackN, yN, nf, NULL, DO_ADD_OVER, dual_cost);

    if (slack != NULL) free(slack);
    if (slackN != NULL) free(slackN);
    if (cost_temp != NULL) free(cost_temp);
    if (ytemp != NULL) free(ytemp);
    if (minus_one != NULL) free(minus_one);
    return SUCCESS_OK;
}

error_t update_theta(real_t *theta) {
    double *theta2 = malloc(DBL_SIZE);
    if (theta2 == NULL) return ALLOCATION_ERROR;
    *theta2 = *theta * *theta;
    *theta = 0.5 * (square_root(*theta2 * *theta2 + 4.0 * *theta2, 0.001, NULL) - *theta2);
    if (theta2 != NULL) free(theta2);
    return SUCCESS_OK;
}

/* calculate_e */
error_t calculate_e(
        creal_t *FN,
        creal_t *qN,
        creal_t *K,
        creal_t *L,
        creal_t *C,
        creal_t *s,
        creal_t *w,
        creal_t *wN,
        cus_t isempty_cmin,
        cus_t isempty_cmax,
        cus_t isempty_xmin,
        cus_t isempty_xmax,
        cus_t isempty_umin,
        cus_t isempty_umax,
        cus_t nx,
        cus_t nu,
        cus_t nc,
        cus_t nf,
        cus_t N,
        real_t *e
        ) {
    /* Declarations: */
    us_t k, i, do_add_CtildeY = FALSE, nw = 0;

    error_t err;
    double minus_one = -1;
    /* End of Declarations*/

    /* e(first_column) = 0.0; */
    memset(e, 0, nx * DBL_SIZE);

    /* Check input-output */
    if (FN == NULL || L == NULL || w == NULL || wN == NULL) {
        return handle_error(NULL_INPUT, "Null input/calculate_e");
    }
    if (e == NULL) handle_error(NULL_OUTPUT, "Null output (e)/calculate_e");
    if (isempty_cmax == FALSE) {
        do_add_CtildeY = TRUE;
        nw += nc;
    }
    if (isempty_cmin == FALSE) {
        do_add_CtildeY = TRUE;
        nw += nc;
    }
    if (isempty_xmax == FALSE) {
        do_add_CtildeY = TRUE;
        nw += nx;
    }
    if (isempty_xmin == FALSE) {
        do_add_CtildeY = TRUE;
        nw += nx;
    }
    if (isempty_umax == FALSE) {
        do_add_CtildeY = TRUE;
        nw += nu;
    }
    if (isempty_umin == FALSE) {
        do_add_CtildeY = TRUE;
        nw += nu;
    }
    if (nw == 0) return handle_error(BAD_PARAM_VALUE, "nw==0 (No constraints provided)");

    err = zmvmult_trans(FN, wN, NULL, NULL, nf, nx, DO_NOT_ADD_OVER, e + N * nx);
    if (err != SUCCESS_OK) handle_error(err, "calculate_e/zmvmult_trans");
    /* if ~isempty(qN), e(:,N+1) = e(:,N+1) + qN; end */
    if (qN != NULL) zvadd(nx, 1.0, qN, e + N * nx);
    /* e(:,N+1) = FN' * yN [+ qN]; */
    for (k = N; k > 1; k--) { /* THE FOR LOOP */
        /* e(:, k) = L * e(:, k+1); */
        zmvmult(L, e + nx*k, NULL, NULL, nx, nx, DO_NOT_ADD_OVER, e + nx * (k - 1));
        if (s != NULL) zvadd(nx, 1.0, s, e + nx * (k - 1));
        if (do_add_CtildeY == TRUE) {
            /* w(:,k) <--> w + (k-1)*nw  */
            i = (k - 1) * nw;
            if (isempty_cmax == FALSE) {
                err = zmvmult(C, w + i, NULL, NULL, nx, nc, DO_ADD_OVER, e + nx * (k - 1));
                if (err != SUCCESS_OK) return handle_error(err, "calculate_e/zmvmult");
                i += nc;
            }
            if (isempty_cmin == FALSE) {
                err = zmvmult(C, w + i, &minus_one, NULL, nx, nc, DO_ADD_OVER, e + nx * (k - 1));
                if (err != SUCCESS_OK) return handle_error(err, "calculate_e/zmvmult");
                i += nc;
            }
            if (isempty_xmax == FALSE) {
                err = zvadd(nx, 1.0, w + i, e + nx * (k - 1));
                if (err != SUCCESS_OK) return handle_error(err, "calculate_e/zvadd");
                i += nx;
            }
            if (isempty_xmin == FALSE) {
                err = zvadd(nx, -1.0, w + i, e + nx * (k - 1));
                if (err != SUCCESS_OK) return handle_error(err, "calculate_e/zvadd");
                i += nx;
            }
            if (isempty_umax == FALSE) {
                err = zmvmult_trans(K, w + i, NULL, NULL, nu, nx, DO_ADD_OVER, e + nx * (k - 1));
                if (err != SUCCESS_OK) return handle_error(err, "calculate_e/zmvmult_trans");
                i += nu;
            }
            if (isempty_umin == FALSE) {
                err = zmvmult_trans(K, w + i, &minus_one, NULL, nu, nx, DO_ADD_OVER, e + nx * (k - 1));
                if (err != SUCCESS_OK) return handle_error(err, "calculate_e/zmvmult_trans");
            }
            if (err != SUCCESS_OK) return handle_error(err, "calculate_e/zvadd");
        } /* END-OF : if (do_add_CtildeY == TRUE)*/
    }
    return SUCCESS_OK;
}

error_t primal_feasibility(
        creal_t *F,
        creal_t *G,
        creal_t *cmin,
        creal_t *cmax,
        creal_t *umin,
        creal_t *umax,
        creal_t *xmin,
        creal_t *xmax,
        creal_t *x,
        creal_t *u,
        cus_t nx,
        cus_t nu,
        cus_t nc,
        cus_t np,
        real_t *g
        ) {
    /* Declarations */
    us_t i = 0, j = 0;
    real_t *FxGu = NULL; /* Store Fx + Gu*/
    error_t err;
    /* End of Declarations */

    if ((F == NULL && G == NULL && cmin == NULL && cmax == NULL && umin == NULL && umax == NULL) || np == 0) {
        g = NULL;
        return SUCCESS_OK;
    } else if (g == NULL) {
        return handle_error(NULL_OUTPUT, "g==NULL/primal_feasibility");
    }

    if (cmax != NULL || cmin != NULL) {
        FxGu = malloc(nc * DBL_SIZE);
        err = zmvmult(F, x, NULL, NULL, nc, nx, DO_NOT_ADD_OVER, FxGu);
        err = zmvmult(G, u, NULL, NULL, nc, nu, DO_ADD_OVER, FxGu);
        if (cmax != NULL) {
            for (j = 0; j < nc; j++) {
                g[j] = FxGu[j] - cmax[j];
            }
            i += nc;
        }
        if (cmin != NULL) {
            for (j = 0; j < nc; j++) {
                g[j + i] = cmin[j] - FxGu[j];
            }
            i += nc;
        }
    }
    if (xmax != NULL) {
        for (j = 0; j < nx; j++) {
            g[j + i] = x[j] - xmax[j];
        }
        i += nx;
    }
    if (xmin != NULL) {
        for (j = 0; j < nx; j++) {
            g[j + i] = xmin[j] - x[j];
        }
        i += nx;
    }
    if (umax != NULL) {
        for (j = 0; j < nu; j++) {
            g[j + i] = u[j] - umax[j];
        }
        i += nu;
    }
    if (umin != NULL) {
        for (j = 0; j < nu; j++) {
            g[j + i] = umin[j] - u[j];
        }
    }
    if (FxGu != NULL) free(FxGu);
    return SUCCESS_OK;
}

error_t dgrad(
        creal_t *A,
        creal_t *B,
        creal_t *f,
        creal_t *Q,
        creal_t *R,
        creal_t *S,
        creal_t *q,
        creal_t *r,
        creal_t *QN,
        creal_t *qN,
        creal_t *FN,
        creal_t *gN,
        creal_t *F,
        creal_t *G,
        creal_t *cmin,
        creal_t *cmax,
        creal_t *xmin,
        creal_t *xmax,
        creal_t *umin,
        creal_t *umax,
        creal_t *D,
        creal_t *L,
        creal_t *M,
        creal_t *K,
        creal_t *C,
        creal_t *s,
        creal_t *d,
        creal_t *Rchol,
        creal_t *w,
        creal_t *wN,
        creal_t *x0,
        creal_t alpha,
        cus_t N,
        cus_t nx,
        cus_t nu,
        cus_t nc,
        cus_t nf,
        real_t *Xstar,
        real_t *Ustar,
        real_t *slack,
        real_t *slackN,
        real_t *y_new,
        real_t *yN_new
        ) {

    /* Declarations */
    error_t err;
    us_t k, i, np = 0;
    real_t *e = NULL, *minus_one = NULL, *ytemp = NULL;
    cus_t isempty_cmin = (cmin == NULL), isempty_cmax = (cmax == NULL), isempty_xmin = (xmin == NULL),
            isempty_xmax = (xmax == NULL), isempty_umin = (umin == NULL), isempty_umax = (umax == NULL);
    /* End of Declarations */

    /* Check input variables: */
    if (A == NULL) return handle_error(NULL_INPUT, "A==NULL/dgrad");
    if (B == NULL) return handle_error(NULL_INPUT, "B==NULL/dgrad");
    if (K == NULL) return handle_error(NULL_INPUT, "K==NULL/dgrad");
    if (M == NULL) return handle_error(NULL_INPUT, "M==NULL/dgrad");
    if (FN == NULL) return handle_error(NULL_INPUT, "FN==NULL/dgrad");
    if (gN == NULL) return handle_error(NULL_INPUT, "gN==NULL/dgrad");
    if (Q == NULL) return handle_error(NULL_INPUT, "Q==NULL/dgrad");
    if (R == NULL) return handle_error(NULL_INPUT, "R==NULL/dgrad");
    if (w == NULL) return handle_error(NULL_INPUT, "y==NULL/dgrad");
    if (wN == NULL) return handle_error(NULL_INPUT, "yN==NULL/dgrad");
    if (x0 == NULL) return handle_error(NULL_INPUT, "x0==NULL/dgrad");
    if (Rchol == NULL) return handle_error(NULL_INPUT, "Rchol==NULL/dgrad");

    if (Xstar == NULL) return handle_error(NULL_OUTPUT, "Xstar==NULL/dgrad");
    if (Ustar == NULL) return handle_error(NULL_OUTPUT, "Ustar==NULL/dgrad");
    if (slack == NULL) return handle_error(NULL_OUTPUT, "slack==NULL/dgrad");
    if (slackN == NULL) return handle_error(NULL_OUTPUT, "slackN==NULL/dgrad");
    if (y_new == NULL) return handle_error(NULL_OUTPUT, "y_new==NULL/dgrad");
    if (yN_new == NULL) return handle_error(NULL_OUTPUT, "yN_new==NULL/dgrad");


    np = nc * (2 - isempty_cmax - isempty_cmin) + nx * (2 - isempty_xmax - isempty_xmin) +
            nu * (2 - isempty_umax - isempty_umin);
    if (np == 0) return handle_error(BAD_PARAM_VALUE, "np==0 (No constraints)");


    /* CALCULATE THE MATRIX e : */
    e = malloc(nx * (N + 1) * DBL_SIZE);
    err = calculate_e(FN, qN, K, L, C, s, w, wN, isempty_cmin, isempty_cmax,
            isempty_xmin, isempty_xmax, isempty_umin, isempty_umax, nx, nu, nc, nf, N, e);
    if (err != SUCCESS_OK) return handle_error(err, "dgrad/calculate_e");

    /* Solution Step: */
    memcpy(Xstar, x0, nx * DBL_SIZE);
    minus_one = malloc(DBL_SIZE);
    *minus_one = -1.0;

    ytemp = malloc(nu * DBL_SIZE);
    if (ytemp == NULL) return handle_error(ALLOCATION_ERROR, "dgrad/malloc ytemp");


    for (k = 0; k < N; k++) {
        err = zmvmult(K, Xstar + k*nx, NULL, NULL, nu, nx, DO_NOT_ADD_OVER, Ustar + k * nu); /* u = Kx */
        if (err != SUCCESS_OK) return handle_error(err, "dgrad/zmvmult$1");
        err = zmvmult(M, e + (k + 1) * nx, NULL, NULL, nu, nx, DO_ADD_OVER, Ustar + k * nu); /* u += M*e(k+1)*/
        if (err != SUCCESS_OK) return handle_error(err, "dgrad/zmvmult$2");
        i = 0;
        if (isempty_cmax == FALSE) {
            err = zmvmult(D, w + k * np + i, NULL, NULL, nu, nc, DO_ADD_OVER, Ustar + k * nu);
            if (err != SUCCESS_OK) return handle_error(err, "dgrad/zmvmult$3");
            i += nc;
        }
        if (isempty_cmin == FALSE) {
            err = zmvmult(D, w + k * np + i, minus_one, NULL, nu, nc, DO_ADD_OVER, Ustar + k * nu);
            if (err != SUCCESS_OK) return handle_error(err, "dgrad/zmvmult$4");
            i += nc;
        }
        i += nx * (2 - isempty_xmax - isempty_xmax);
        if (isempty_umax == FALSE) {
            err = zchlslv_cpy(Rchol, w + k * np + i, ytemp, nu);
            if (err != SUCCESS_OK) return handle_error(err, "dgrad/zchlslv_cpy$5");
            if (zvadd(nu, -1.0, ytemp, Ustar + k * nu) != SUCCESS_OK) return handle_error(err, "dgrad/zvadd$6");
            i += nu;
        }
        if (isempty_umin == FALSE) {
            zchlslv_cpy(Rchol, w + k * np + i, ytemp, nu);
            if (zvadd(nu, 1.0, ytemp, Ustar + k * nu) != SUCCESS_OK) return handle_error(err, "dgrad/zvadd$7");
        }

        primal_feasibility(F, G, cmin, cmax, umin, umax, xmin, xmax, Xstar + k*nx, Ustar + k*nu, nx, nu, nc, np, slack + k * np);
        if (slack + k * np == NULL) return handle_error(NULL_OUTPUT, "dgrad/primal_feasibility/slack==NULL$8");
        state_update(A, B, f, Xstar + k*nx, Ustar + k*nu, nx, nu, Xstar + (k + 1) * nx);

        for (i = 0; i < np; i++) {
            y_new[k * np + i] = xmax(w[k * np + i] + alpha * slack[k * np + i], 0.0);
        }
    }

    zmvmult(FN, Xstar + N*nx, NULL, NULL, nf, nx, DO_NOT_ADD_OVER, slackN);
    for (i = 0; i < nf; i++) {
        slackN[i] -= gN[i];
    }
    for (i = 0; i < nf; i++) {
        yN_new[i] = xmax(wN[i] + alpha * slackN[i], 0.0);
    }

    if (e != NULL) free(e);
    if (ytemp != NULL) free(ytemp);
    if (minus_one != NULL) free(minus_one);
    return SUCCESS_OK;
}



#ifdef GPAD_DEBUG

__INLINE__ void matlabPrint(const short mRows, const short nCols, const double *data) {
    mxArray * out[1], *in[1];
    mxArray *the_mxArray;
    if (data == 0) return;
    the_mxArray = mxCreateDoubleMatrix(mRows, nCols, mxREAL);
    memcpy(mxGetPr(the_mxArray), data, nCols * mRows * sizeof (double));
    in[0] = the_mxArray;
    mexCallMATLAB(0, out, 1, in, "disp"); /* call disp from Matlab */
    mxDestroyArray(the_mxArray);
}
#endif



/* See http://rosettacode.org/wiki/Cholesky_decomposition#C 
   Useful for the factor step*/

/*
real_t *cholesky(real_t *A, cus_t n) {
    real_t *L, s;
    us_t i, j, k;
    
    L = calloc(n * n, DBL_SIZE);
    for (i = 0; i < n; i++)
        for (j = 0; j < (i + 1); j++) {
            s = 0.0;
            for (k = 0; k < j; k++) { s += L[i * n + k] * L[j * n + k]; }
            if (i == j) {
                L[i * n + j] = sqrt(A[i * n + i] - s);
            } else {
                L[i * n + j] = (1.0 / L[j * n + j] * (A[i * n + j] - s));
            }
        }
    return L;
}
 */


error_t gpad_control(
        creal_t *A,
        creal_t *B,
        creal_t *f,
        creal_t *Q,
        creal_t *R,
        creal_t *S,
        creal_t *q,
        creal_t *r,
        creal_t *QN,
        creal_t *qN,
        creal_t *FN,
        creal_t *gN,
        creal_t *F,
        creal_t *G,
        creal_t *cmin,
        creal_t *cmax,
        creal_t *xmin,
        creal_t *xmax,
        creal_t *umin,
        creal_t *umax,
        creal_t *D,
        creal_t *L,
        creal_t *M,
        creal_t *K,
        creal_t *C,
        creal_t *s,
        creal_t *d,
        creal_t *Rchol,
        creal_t alpha,
        creal_t *x0,
        creal_t epsilon_g,
        creal_t epsilon_V,
        real_t *Xstar,
        real_t *Ustar,
        cus_t nx,
        cus_t nu,
        cus_t nc,
        cus_t nf,
        cus_t N,
        cus_t classQ,
        cus_t classR,
        cus_t max_iter,
        cus_t enforce_monotonicity,
        diagnostics_t *diagnostics
        ) {

    us_t np, doloop = TRUE, iter = 0, i;
    error_t err;
    real_t theta = 1.0, theta_p = 1.0, theta2 = 0.0, beta = 0.0, max_viol = VERY_BIG_DBL, max_viol_bar = VERY_BIG_DBL, gap;
    real_t *y = NULL, *yN = NULL, *y_p = NULL, *yN_p = NULL, *slack = NULL, *slackN = NULL, *slack_bar = NULL,
            *slackN_bar = NULL, *w = NULL, *wN = NULL, *minus_one = NULL, *primal_cost = NULL, *dual_cost, rel_gap = VERY_BIG_DBL;
    cus_t isempty_cmin = (cmin == NULL), isempty_cmax = (cmax == NULL), isempty_xmin = (xmin == NULL),
            isempty_xmax = (xmax == NULL), isempty_umin = (umin == NULL), isempty_umax = (umax == NULL);

    minus_one = malloc(DBL_SIZE);
    if (minus_one == NULL) return handle_error(ALLOCATION_ERROR, "gpad$95");
    primal_cost = malloc(DBL_SIZE);
    if (primal_cost == NULL) return handle_error(ALLOCATION_ERROR, "gpad$96");
    dual_cost = malloc(DBL_SIZE);
    if (dual_cost == NULL) return handle_error(ALLOCATION_ERROR, "gpad$97");
    *minus_one = -1.0;
    *primal_cost = 0.0;
    *dual_cost = 0.0;


    np = nc * (2 - isempty_cmax - isempty_cmin) + nx * (2 - isempty_xmax - isempty_xmin) +
            nu * (2 - isempty_umax - isempty_umin);

    y = calloc(np*N, DBL_SIZE);
    if (y == NULL) return handle_error(ALLOCATION_ERROR, "gpad$101");
    w = malloc(np * N * DBL_SIZE);
    if (w == NULL) return handle_error(ALLOCATION_ERROR, "gpad$102");
    y_p = malloc(np * N * DBL_SIZE);
    if (y_p == NULL) return handle_error(ALLOCATION_ERROR, "gpad$103");
    yN = calloc(nf, DBL_SIZE);
    if (yN == NULL) return handle_error(ALLOCATION_ERROR, "gpad$104");
    wN = malloc(nf * DBL_SIZE);
    if (wN == NULL) return handle_error(ALLOCATION_ERROR, "gpad$105");
    yN_p = malloc(nf * DBL_SIZE);
    if (yN_p == NULL) return handle_error(ALLOCATION_ERROR, "gpad$106");
    slack = calloc(np * N, DBL_SIZE);
    if (slack == NULL) return handle_error(ALLOCATION_ERROR, "gpad$107");
    slackN = calloc(nf, DBL_SIZE);
    if (slackN == NULL) return handle_error(ALLOCATION_ERROR, "gpad$108");
    slack_bar = calloc(np * N, DBL_SIZE);
    slackN_bar = calloc(nf, DBL_SIZE);

#if (defined(GPAD_DEBUG) && GPAD_DEBUG >= DEBUG_VERBOSE)
    printf("log(MV)         log(MVB)\n");
#endif
    if (diagnostics != NULL) {
        diagnostics -> msg = "Maximum Iterations reached";
        diagnostics -> point_1 = 0;
        diagnostics -> point_2 = 0;
        diagnostics -> point_3 = 0;
    }

    while (doloop == TRUE && iter <= max_iter) {
        /* Update w: */
        beta = theta * (1 / theta_p - 1);
        for (i = 0; i < np * N; i++) {
            w[i] = y[i] + beta * (y[i] - y_p[i]);
        }
        for (i = 0; i < nf; i++) {
            wN[i] = yN[i] + beta * (yN[i] - yN_p[i]);
        }
        memcpy(y_p, y, np * N * DBL_SIZE);
        memcpy(yN_p, yN, nf * DBL_SIZE);


        err = dgrad(A, B, f, Q, R, S, q, r, QN, qN, FN, gN, F, G, cmin, cmax,
                xmin, xmax, umin, umax, D, L, M, K, C, s, d, Rchol, y, yN, x0, alpha, N,
                nx, nu, nc, nf, Xstar, Ustar, slack, slackN, y, yN);
        if (err != SUCCESS_OK) return handle_error(err, "gpad/dgrad$109");
        max_viol = xmax(max(slack, np * N), max(slackN, nf));

        /* Check the termination criterion */
        /* slack_bar  = (1-theta)*slack_bar+theta*slack; */
        for (i = 0; i < np * N; i++) {
            slack_bar[i] = (1 - theta) * slack_bar[i] + theta * slack[i];
        }
        for (i = 0; i < nf; i++) {
            slackN_bar[i] = (1 - theta) * slackN_bar[i] + theta * slackN[i];
        }
        max_viol_bar = xmax(max(slack_bar, np * N), max(slackN_bar, nf));
        if (diagnostics != NULL
                && diagnostics -> log_max_viol != NULL
                && diagnostics -> log_max_viol_bar != NULL) {
            if (max_viol <= 0) {
                *((diagnostics -> log_max_viol) + iter) = -123.456;
            } else {
                *((diagnostics -> log_max_viol) + iter) = log(max_viol);
            }
            if (max_viol_bar <= 0) {
                *((diagnostics -> log_max_viol_bar) + iter) = -123.456;
            } else {
                *((diagnostics -> log_max_viol_bar) + iter) = max_viol_bar;
            }
        }
#if (defined(GPAD_DEBUG) && GPAD_DEBUG >= DEBUG_VERBOSE)
        printf("%7.3f%15.3f\n", log(max_viol), log(max_viol_bar));
#endif
        if (max_viol_bar <= epsilon_g) { /* POINT #1 */
            if (diagnostics != NULL) diagnostics -> msg = "max_viol_bar <= eg";
            doloop = FALSE;
        } else if (max_viol <= epsilon_g) { /* POINT #2 */
            if (min(wN, nf) >= 0 && min(w, np * N) >= 0) {
                if (diagnostics != NULL) diagnostics -> point_1++;
                gap = 0.0;
                for (i = 0; i < N; i++) {
                    zcdot(slack + i * np, w + i * np, np, minus_one, DO_ADD_OVER, &gap);
                }
                zcdot(slackN, wN, nf, minus_one, DO_ADD_OVER, &gap);
                if (gap <= epsilon_V) { /* POINT #3 */
                    if (diagnostics != NULL) diagnostics -> msg = "gap <= eV";
                    doloop = FALSE;
                } else { /* POINT #4 */
                    if (diagnostics != NULL) diagnostics -> point_2++;
                    err = pcost(Q, R, S, q, r, QN, qN, N, Xstar, Ustar, nx, nu, classQ, classR, primal_cost);
                    if (SUCCESS_OK != err) return handle_error(err, "gpad/dgrad$188");
                    if (gap <= (epsilon_V / (1 + epsilon_V))* (*primal_cost)) {
                        if (diagnostics != NULL) diagnostics -> msg = "gap <= eV/(1+eV)*primal_cost";
                        doloop = FALSE;
                    }
                }
            } else { /* POINT #5: if min(wN)>=0  && min(min(w))>=0 */
                if (diagnostics != NULL) diagnostics -> point_3++;
                err = pcost(Q, R, S, q, r, QN, qN, N, Xstar, Ustar, nx, nu, classQ, classR, primal_cost);
                if (SUCCESS_OK != err) return handle_error(err, "gpad/pcost$201");
                err = dcost(A, B, f, Q, R, S, q, r, QN, qN, FN, gN, F, G, cmin, cmax, xmin, xmax, umin,
                        umax, D, L, M, K, C, s, d, Rchol, y, yN, x0, N, nx, nu, nc, nf, classQ, classR,
                        Xstar, Ustar, dual_cost);
                if (SUCCESS_OK != err) return handle_error(err, "gpad/dcost$202");
                rel_gap = (*primal_cost - *dual_cost) / xmax(*dual_cost, 1.00);
                if (rel_gap <= epsilon_V) {
                    doloop = FALSE;
                    if (diagnostics != NULL) diagnostics -> msg = "real_gap <= eV";
                }
            }
        }
        /* Update theta */
        theta_p = theta;
        theta2 = theta*theta;
        theta = (square_root(theta2 * theta2 + 4 * theta2, 0.001, NULL) - theta2) / 2;
        iter++;
    }
    if (primal_cost != NULL) free(primal_cost);
    if (minus_one != NULL) free(minus_one);
#if (defined(GPAD_DEBUG) && GPAD_DEBUG >= DEBUG_VERBOSE)
    printf("iterations = %d\n", iter);
#endif
    if (diagnostics != NULL) {
        diagnostics -> iters = iter - 1;
    }
    return NOT_IMPLEMENTED_YET;
}
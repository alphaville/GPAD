/* 
 * File: mx_gpad.c
 * Created on: 5 July, 2013
 * Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
 * Institute: IMT Lucca, Lucca Italy
 *    __  __      __  
 *   / _ |__) /\ |  \ 
 *   \__)|   /--\|__/ 
 *
 */
#include "mex.h"
#include "../gpad.h"

/* [Ustar, Xstar, Diagnostics] = gpad(A, B, f, Q, R, S, q, r, QN, qN, ...
 *                                       FN, gN, F, G, cmin, cmax, xmin, xmax, umin, ...
 *                                       umax, D, L, M, K, C, s, d, Rchol, ...
 *                                       x0, N, alpha, eps_g, eps_V, max_iter); */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

    /* Declarations */

    us_t N, nx = 0, nu = 0, nc = 0, nf = 0, np = 0, i;
    const mxArray *A_mx = prhs[0],
            *B_mx = prhs[1],
            *f_mx = prhs[2],
            *Q_mx = prhs[3],
            *R_mx = prhs[4],
            *S_mx = prhs[5],
            *q_mx = prhs[6],
            *r_mx = prhs[7],
            *QN_mx = prhs[8],
            *qN_mx = prhs[9],
            *FN_mx = prhs[10],
            *gN_mx = prhs[11],
            *F_mx = prhs[12],
            *G_mx = prhs[13],
            *cmin_mx = prhs[14],
            *cmax_mx = prhs[15],
            *xmin_mx = prhs[16],
            *xmax_mx = prhs[17],
            *umin_mx = prhs[18],
            *umax_mx = prhs[19],
            *D_mx = prhs[20],
            *L_mx = prhs[21],
            *M_mx = prhs[22],
            *K_mx = prhs[23],
            *C_mx = prhs[24],
            *s_mx = prhs[25],
            *d_mx = prhs[26],
            *Rchol_mx = prhs[27],
            *x0_mx = prhs[28],
            *N_mx = prhs[29],
            *alpha_mx = prhs[30],
            *epsilon_g_mx = prhs[31],
            *epsilon_V_mx = prhs[32],
            *max_iter_mx = prhs[33];
    real_t *A = NULL, *B = NULL, *f = NULL, *Q = NULL, *R = NULL, *S = NULL, *q = NULL, *r = NULL, *QN = NULL, *qN = NULL, *FN = NULL, *gN = NULL,
            *F = NULL, *G = NULL, *cmin = NULL, *cmax = NULL, *xmin = NULL, *xmax = NULL, *umin = NULL, *umax = NULL, *D = NULL, *L = NULL, *M = NULL,
            *K = NULL, *C = NULL, *s = NULL, *d = NULL, *Rchol = NULL, *x0 = NULL, *Ndbl = NULL, *iter_ptr = NULL, *point1_ptr = NULL, *point2_ptr = NULL,
            *point3_ptr = NULL, *epsilon_g = NULL, *epsilon_V = NULL, *Xstar = NULL, *Ustar = NULL, *alpha = NULL, *out = NULL, *max_iter = NULL,
            *log_max_viol_ptr = NULL, *log_max_viol_bar_ptr = NULL;
    const char *fieldnames[2];
    mxArray *iter_MATLAB = NULL, *msg_MATLAB = NULL, *point1_MATLAB = NULL, *point2_MATLAB = NULL, *point3_MATLAB = NULL,
            *log_max_viol_MATLAB = NULL, *log_max_viol_bar_MATLAB = NULL;
    diagnostics_t *diagn = NULL;
    us_t num_diagn_fields = 5;

    /* End of Declarations */

    if (nrhs != 34) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "gpad requires exactly 32 input arguments");

    A = mxGetPr(A_mx);
    B = mxGetPr(B_mx);
    f = mxGetPr(f_mx);
    Q = mxGetPr(Q_mx);
    R = mxGetPr(R_mx);
    S = mxGetPr(S_mx);
    q = mxGetPr(q_mx);
    r = mxGetPr(r_mx);
    QN = mxGetPr(QN_mx);
    qN = mxGetPr(qN_mx);
    FN = mxGetPr(FN_mx);
    gN = mxGetPr(gN_mx);
    F = mxGetPr(F_mx);
    G = mxGetPr(G_mx);
    cmin = mxGetPr(cmin_mx);
    cmax = mxGetPr(cmax_mx);
    xmin = mxGetPr(xmin_mx);
    xmax = mxGetPr(xmax_mx);
    umin = mxGetPr(umin_mx);
    umax = mxGetPr(umax_mx);
    D = mxGetPr(D_mx);
    L = mxGetPr(L_mx);
    M = mxGetPr(M_mx);
    K = mxGetPr(K_mx);
    C = mxGetPr(C_mx);
    s = mxGetPr(s_mx);
    d = mxGetPr(d_mx);
    Rchol = mxGetPr(Rchol_mx);
    x0 = mxGetPr(x0_mx);
    Ndbl = mxGetPr(N_mx);
    alpha = mxGetPr(alpha_mx);
    epsilon_g = mxGetPr(epsilon_g_mx);
    epsilon_V = mxGetPr(epsilon_V_mx);
    max_iter = mxGetPr(max_iter_mx);

    /* Check for null input variables */
    if (epsilon_g == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "epsilon_g cannot be empty");
    if (epsilon_V == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "epsilon_V cannot be empty");
    if (A == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "A cannot be empty");
    if (B == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "B cannot be empty");
    if (Q == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "Q cannot be empty");
    if (R == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "R cannot be empty");
    if (QN == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "QN cannot be empty");
    if (FN == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "FN cannot be empty");
    if (gN == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "gN cannot be empty");
    if (Ndbl == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "The prediction horizon N cannot be empty");
    if (alpha == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "alpha cannot be empty");
    if (x0 == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:badInput", "x0 cannot be empty");

    /* Get Sizes */
    N = (cus_t) round(* Ndbl);
    nx = mxGetM(A_mx);
    nu = mxGetN(B_mx);
    if (F != NULL) nc = mxGetM(F_mx);
    nf = mxGetM(FN_mx);
    np = nc * (2 - (cmin == NULL)-(cmax == NULL)) + nx * (2 - (xmin == NULL)-(xmax == NULL)) +
            nu * (2 - (umin == NULL)-(umax == NULL));

    /* Allocation */
    Xstar = malloc(nx * (N + 1) * DBL_SIZE);
    if (Xstar == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:allocationError", "Xstar cannot be allocated");
    memcpy(Xstar, x0, nx * DBL_SIZE);
    Ustar = malloc(nu * N * DBL_SIZE);
    if (Ustar == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:allocationError", "Ustar cannot be allocated");

    /* Do le magic */
    if (nlhs == 3) {
        diagn = malloc(sizeof (diagnostics_t));
        if (diagn == NULL) mexErrMsgIdAndTxt("MATLAB:gpad:allocationError", "diagn cannot be allocated");
        diagn -> msg = malloc(50 * sizeof (char));
        diagn -> log_max_viol = calloc((us_t) (*max_iter) + 1, DBL_SIZE);
        diagn -> log_max_viol_bar = calloc((us_t) (*max_iter) + 1, DBL_SIZE);
    }

    gpad_control(A, B, f, Q, R, S, q, r, QN, qN, FN, gN, F, G, cmin, cmax, xmin, xmax,
            umin, umax, D, L, M, K, C, s, d, Rchol, *alpha, x0, *epsilon_g, *epsilon_V, Xstar,
            Ustar, nx, nu, nc, nf, N, MAT_DENSE, MAT_DENSE, (us_t) (*max_iter), FALSE, diagn);

    if (nlhs == 0) return;
    if (nlhs >= 1) {
        plhs[0] = mxCreateDoubleMatrix(nu, N, mxREAL);
        out = mxGetPr(plhs[0]);
        memcpy(out, Ustar, nu * N * DBL_SIZE);
    }
    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleMatrix(nx, N + 1, mxREAL);
        out = mxGetPr(plhs[1]);
        memcpy(out, Xstar, nx * (N + 1) * DBL_SIZE);
    }

    if (nlhs == 3 && diagn != NULL) {
        if (diagn -> iters > 0 && diagn -> log_max_viol != NULL && diagn -> log_max_viol_bar != NULL) {
            num_diagn_fields += 2;
        }
        for (i = 0; i < num_diagn_fields; i++) fieldnames[i] = (char*) mxMalloc(20);


        iter_MATLAB = mxCreateDoubleMatrix(1, 1, mxREAL);
        point1_MATLAB = mxCreateDoubleMatrix(1, 1, mxREAL);
        point2_MATLAB = mxCreateDoubleMatrix(1, 1, mxREAL);
        point3_MATLAB = mxCreateDoubleMatrix(1, 1, mxREAL);
        if (diagn -> iters > 0 && diagn -> log_max_viol != NULL && diagn -> log_max_viol_bar != NULL) {
            log_max_viol_MATLAB = mxCreateDoubleMatrix(diagn -> iters, 1, mxREAL);
            log_max_viol_bar_MATLAB = mxCreateDoubleMatrix(diagn -> iters, 1, mxREAL);
        }

        msg_MATLAB = mxCreateString(diagn->msg);

        iter_ptr = mxGetPr(iter_MATLAB);
        *iter_ptr = (real_t) diagn -> iters;

        point1_ptr = mxGetPr(point1_MATLAB);
        *point1_ptr = (real_t) diagn -> point_1;

        point2_ptr = mxGetPr(point2_MATLAB);
        *point2_ptr = (real_t) diagn -> point_2;

        point3_ptr = mxGetPr(point3_MATLAB);
        *point3_ptr = (real_t) diagn -> point_3;

        if (diagn -> iters > 0 && diagn -> log_max_viol != NULL && diagn -> log_max_viol_bar != NULL) {
            log_max_viol_ptr = mxGetPr(log_max_viol_MATLAB);
            log_max_viol_bar_ptr = mxGetPr(log_max_viol_bar_MATLAB);
            for (i = 0; i < diagn -> iters; i++) {
                log_max_viol_ptr[i] = diagn -> log_max_viol[i];
                log_max_viol_bar_ptr[i] = diagn -> log_max_viol_bar[i];
            }
        }

        memcpy((void *) fieldnames[0], "iter", sizeof ("iter"));
        memcpy((void *) fieldnames[1], "message", sizeof ("message"));
        memcpy((void *) fieldnames[2], "point1", sizeof ("point1"));
        memcpy((void *) fieldnames[3], "point2", sizeof ("point2"));
        memcpy((void *) fieldnames[4], "point3", sizeof ("point3"));
        if (diagn -> iters > 0 && diagn -> log_max_viol != NULL && diagn -> log_max_viol_bar != NULL) {
            memcpy((void *) fieldnames[5], "log_max_viol", sizeof ("log_max_viol"));
            memcpy((void *) fieldnames[6], "log_max_viol_bar", sizeof ("log_max_viol_bar"));
        }

        plhs[2] = mxCreateStructMatrix(1, 1, num_diagn_fields, fieldnames);

        for (i = 0; i < num_diagn_fields; i++) mxFree((void *) fieldnames[i]);

        mxSetFieldByNumber(plhs[2], 0, 0, iter_MATLAB);
        mxSetFieldByNumber(plhs[2], 0, 1, msg_MATLAB);
        mxSetFieldByNumber(plhs[2], 0, 2, point1_MATLAB);
        mxSetFieldByNumber(plhs[2], 0, 3, point2_MATLAB);
        mxSetFieldByNumber(plhs[2], 0, 4, point3_MATLAB);
        if (diagn -> iters > 0 && diagn -> log_max_viol != NULL && diagn -> log_max_viol_bar != NULL) {
            mxSetFieldByNumber(plhs[2], 0, 5, log_max_viol_MATLAB);
            mxSetFieldByNumber(plhs[2], 0, 6, log_max_viol_bar_MATLAB);
        }


    }


}


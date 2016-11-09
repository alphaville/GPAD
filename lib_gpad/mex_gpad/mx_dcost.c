/* 
 * File: mx_dcost.c
 * Created on: 15 June, 2013
 * Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
 * Institute: IMT Lucca, Lucca Italy
 *    __  __      __  
 *   / _ |__) /\ |  \ 
 *   \__)|   /--\|__/ 
 *
 */
#include "mex.h"
#include "../gpad.h"
#include <math.h>
#include <stdio.h>


/* [d, Ustar, Xstar] = dcost(A, B, f, Q, R, S, q, r, QN, qN, FN, gN, F, G, cmin, cmax, xmin, xmax, umin, umax, D, L, M, K, C, s, d, Rchol, y, yN, x0, N); */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

    /* Declarations */
    error_t err;
    char error_message_to_matlab[45];
    us_t N, nx=0, nu=0, nc=0, nf=0;
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
            *y_mx = prhs[28],
            *yN_mx = prhs[29],
            *x0_mx = prhs[30],
            *N_mx = prhs[31];
    real_t *A=NULL, *B=NULL, *f=NULL, *Q=NULL, *R=NULL, *S=NULL, *q=NULL, *r=NULL, *QN=NULL, *qN=NULL, *FN=NULL, *gN=NULL,
            *F=NULL, *G=NULL, *cmin=NULL, *cmax=NULL, *xmin=NULL, *xmax=NULL, *umin=NULL, *umax=NULL, *D=NULL, *L=NULL, *M=NULL,
            *K=NULL, *C=NULL, *s=NULL, *d=NULL, *Rchol=NULL, *y=NULL, *yN=NULL, *x0=NULL, *Ndbl=NULL, *dual_cost=NULL, 
            *Xstar=NULL, *Ustar=NULL, *out=NULL;
    
    /* End of Declarations */

    if (nrhs != 32) mexErrMsgIdAndTxt("MATLAB:dcost:badInput", "dcost requires exactly 32 input arguments"); 

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
    y = mxGetPr(y_mx);
    yN = mxGetPr(yN_mx);
    x0 = mxGetPr(x0_mx);
    Ndbl = mxGetPr(N_mx);
    
    /* Check for null input variables */
    if (A==NULL) mexErrMsgIdAndTxt("MATLAB:dcost:badInput", "A cannot be empty"); 
    if (B==NULL) mexErrMsgIdAndTxt("MATLAB:dcost:badInput", "B cannot be empty"); 
    if (Q==NULL) mexErrMsgIdAndTxt("MATLAB:dcost:badInput", "Q cannot be empty"); 
    if (R==NULL) mexErrMsgIdAndTxt("MATLAB:dcost:badInput", "R cannot be empty"); 
    if (QN==NULL) mexErrMsgIdAndTxt("MATLAB:dcost:badInput", "QN cannot be empty"); 
    if (FN==NULL) mexErrMsgIdAndTxt("MATLAB:dcost:badInput", "FN cannot be empty"); 
    if (gN==NULL) mexErrMsgIdAndTxt("MATLAB:dcost:badInput", "gN cannot be empty"); 
    if (Ndbl==NULL) mexErrMsgIdAndTxt("MATLAB:dcost:badInput", "The prediction horizon N cannot be empty"); 
    
    /* Get Sizes */
    N = (cus_t) round(* Ndbl);
    nx = mxGetM(A_mx);
    nu = mxGetN(B_mx);
    if (F!=NULL) nc = mxGetM(F_mx);
    nf = mxGetM(FN_mx);
    
    Xstar = malloc(nx*(N+1)*DBL_SIZE);
    Ustar = malloc(nu*N*DBL_SIZE);
    
    dual_cost = malloc(DBL_SIZE);
    *dual_cost = 0.0;
    err =  dcost(A, B, f, Q, R, S, q, r, QN, qN, FN, gN, F, G, cmin, cmax, xmin, xmax, 
        umin, umax, D, L, M, K, C, s, d, Rchol, y, yN, x0, N, nx, nu, nc, nf, MAT_DENSE, MAT_DENSE,
        Xstar,Ustar, dual_cost);    
    if (err != SUCCESS_OK){
        sprintf(error_message_to_matlab, "Function dcost - Error %d.", err );
        mexErrMsgIdAndTxt("MATLAB:dcost:failure", error_message_to_matlab);
    }
    
    if (nlhs==0) return;
    if (nlhs>=4) mexErrMsgIdAndTxt("MATLAB:dcost:bad_output", "DCOST outputs up to 3 variables");
    if (nlhs>=1){
    	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    	out = mxGetPr(plhs[0]);
    	memcpy(out, dual_cost, DBL_SIZE);
    }
    if (nlhs>=2){
    	plhs[1] = mxCreateDoubleMatrix(nu, N, mxREAL);
    	out = mxGetPr(plhs[1]);
    	memcpy(out, Ustar, nu*N*DBL_SIZE);
    } else { if (Ustar!=NULL) free(Ustar); }
    if (nlhs==3){
    	plhs[2] = mxCreateDoubleMatrix(nx, N+1, mxREAL);
    	out = mxGetPr(plhs[2]);
    	memcpy(out, Xstar, nx*(N+1)*DBL_SIZE);
    } else { if (Xstar!=NULL) free(Xstar); }

}
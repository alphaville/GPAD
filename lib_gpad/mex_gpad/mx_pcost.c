/* 
 * File: mx_pcost.c
 * Created on: 15 June, 2013
 * Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
 * Institute: IMT Lucca, Lucca Italy
 *    __  __      __  
 *   / _ |__) /\ |  \ 
 *   \__)|   /--\|__/ 
 *
 */
#include "../gpad.h"
#include "mex.h"



/*
 * Syntax (from MATLAB):
 * J = pcost(Q, R, S, q, r, QN, qN, N, X, U, classQ, classR)
 */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

    /* Declarations: */
    error_t error_pcost;
    real_t *Q, *R, *S, *q, *r, *QN, *qN, *X, *U, classQ, classR, *cost=NULL;
    real_t *intTemp = NULL, *output=NULL;
    us_t N, nx, nu;
    real_t *NdataPtr;
    dynsys_t *dyn;
    mpc_t *mpc;
    char error_message_to_matlab[40];
    
    /* Raw input data from MATLAB: */
    const mxArray *Q_mx = prhs[0], *R_mx = prhs[1], *S_mx = prhs[2], *q_mx = prhs[3], *r_mx = prhs[4], *QN_mx = prhs[5],
            *qN_mx = prhs[6], *N_mx = prhs[7], *X_mx = prhs[8], *U_mx = prhs[9], *classQ_mx = prhs[10], *classR_mx = prhs[11];    

    /* Check the input data: */
    if (nlhs == 0) return;
    NdataPtr = (real_t *) mxGetPr(N_mx);
    
    N = (int) (*NdataPtr);
    if (nlhs > 1) mexErrMsgIdAndTxt("MATLAB:pcost:numOutputs", "PCOST exports only 1 variable."
            "\nType 'help pcost' for more information.");
    if (nrhs < 10 || nrhs > 12) mexErrMsgIdAndTxt("MATLAB:pcost_mx:numInputs",
            "PCOST admits 10 to 12 input parameters.");
    if (mxGetN(N_mx) == 0 || mxGetM(N_mx) == 0) mexErrMsgIdAndTxt("MATLAB:pcost_mx:badInput",
            "pcost(Q,R,S,q,r,QN,qN,N,X,U,classQ,classR): N cannot be empty (null)!");    
    if (mxGetN(N_mx) != 1 || mxGetM(N_mx) != 1) mexErrMsgIdAndTxt("MATLAB:pcost_mx:badInput",
            "pcost(Q,R,S,q,r,QN,qN,N,X,U,classQ,classR): N should be a scalar - not a matrix!");
    if (mxGetN(X_mx) < N + 1) mexErrMsgIdAndTxt("MATLAB:pcost:badInput",
            "pcost(Q,R,S,q,r,QN,qN,N,X,U,classQ,classR): The number of columns of X must be >= N+1");
    if (mxGetN(U_mx) < N) mexErrMsgIdAndTxt("MATLAB:pcost:badInput",
            "pcost(Q,R,S,q,r,QN,qN,N,X,U,classQ,classR): The number of columns of U must be >= N+1");
    /* TODO Do more tests  */

    /* Get les input data from MATLAB: */
    Q = mxGetPr(Q_mx);
    R = mxGetPr(R_mx);
    S = mxGetPr(S_mx);
    q = mxGetPr(q_mx); 
    r = mxGetPr(r_mx);
    QN = mxGetPr(QN_mx); 
    qN = mxGetPr(qN_mx);
    X = mxGetPr(X_mx);
    U = mxGetPr(U_mx);
    
    output = 0;
    cost = calloc(1, DBL_SIZE);
    nx = mxGetM(X_mx); nu = mxGetM(U_mx);

    if (nrhs >= 11) {
        intTemp = mxGetPr(classQ_mx);
        classQ = (us_t) *intTemp;
    } else {
        classQ = MAT_DENSE;
    }
    if (nrhs == 12) {
        intTemp = mxGetPr(classR_mx);
        classR = (us_t) * intTemp;
    } else {
        classR = MAT_DENSE;
    }
    intTemp = mxGetPr(N_mx);
    

    /* Do le math: */
    dyn = dynsys_factory_meagre(nx, nu, 0);
    mpc = mpc_factory_simple(Q, R, S, q, r, QN, qN, classQ, classR, (cus_t) N, 0, dyn);
    error_pcost = pcost_i(mpc, X, U, cost);

    if (error_pcost != SUCCESS_OK){
        sprintf(error_message_to_matlab, "Function pcost - Error %d.", error_pcost );
        mexErrMsgIdAndTxt("MATLAB:pcost:failure", error_message_to_matlab);
    }


    /* Output le result back to MATLAB: */
    if (nlhs == 1) {
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        output = mxGetPr(plhs[0]);
        output[0] = *cost;
    }
}
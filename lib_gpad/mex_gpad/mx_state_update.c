/* 
 * File: mx_state_update.c
 * Created on: 15 June, 2013
 * Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
 * Institute: IMT Lucca, Lucca Italy
 *    __  __      __  
 *   / _ |__) /\ |  \ 
 *   \__)|   /--\|__/ 
 */

#include "mex.h"
#include "../gpad.h"

/* x = state_update(A, B , f, x, u); */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    
    /* Declarations */
    const mxArray *A_mx = prhs[0], *B_mx = prhs[1], *f_mx = prhs[2];
    const mxArray *x_mx = prhs[3], *u_mx = prhs[4];
    dynsys_t *sys;
    error_t err;
    real_t *xnew = NULL, *out = NULL, *A = NULL;
    real_t *B = NULL, *f = NULL, *x = NULL, *u = NULL;
    us_t nx, nu, colsA;
    char error_message_to_matlab[40];
    
    if (nrhs!=5) mexErrMsgIdAndTxt("MATLAB:state_update:badInput", "state_update requires exactly 5 input arguments");
    
    nx = (us_t) mxGetM(A_mx); nu = (us_t) mxGetN(B_mx); colsA = (us_t) mxGetN(A_mx);
    if (nx==0 || colsA ==0) mexErrMsgIdAndTxt("MATLAB:state_update:badInput", "A cannot be empty");
    if (nx!=colsA) mexErrMsgIdAndTxt("MATLAB:state_update:badInput", "A must be square");
    if (nu==0 || mxGetM(B_mx)==0) mexErrMsgIdAndTxt("MATLAB:state_update:badInput", "B cannot be empty");
    if (mxGetM(x_mx)==0 || mxGetN(x_mx)==0) mexErrMsgIdAndTxt("MATLAB:state_update:badInput", "x cannot be empty");
    if (mxGetM(u_mx)==0 || mxGetN(u_mx)==0) mexErrMsgIdAndTxt("MATLAB:state_update:badInput", "u cannot be empty");
    if (nu!=mxGetM(u_mx)) mexErrMsgIdAndTxt("MATLAB:state_update:badInput", "B and u are of incompatible dimensions");
    if (nx!=mxGetM(x_mx)) mexErrMsgIdAndTxt("MATLAB:state_update:badInput", "A and x are of incompatible dimensions");
    
    A = mxGetPr(A_mx);
    B = mxGetPr(B_mx);
    f = mxGetPr(f_mx);
    x = mxGetPr(x_mx);
    u = mxGetPr(u_mx);
    
    xnew = malloc(nx * DBL_SIZE);
    if (xnew==NULL) mexErrMsgIdAndTxt("MATLAB:state_update:allocation", "Cannot allocate memory for xnew");
    sys = dynsys_factory(A, B, f, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, (cus_t) nx, (cus_t) nu, 0);
    err = state_update_i(sys, x, u, xnew);    
    if (err!=SUCCESS_OK){
        sprintf(error_message_to_matlab, "Function state_update_i - Error %d.", err );
        mexErrMsgIdAndTxt("MATLAB:state_update_i:error", error_message_to_matlab);
    }
    if (nlhs >= 1) {
        plhs[0] = mxCreateDoubleMatrix(nx, 1, mxREAL);
        out = mxGetPr(plhs[0]);
        memcpy(out, xnew, nx*DBL_SIZE);
    }
}
/* 
 * File: mx_calculate_e.c
 * Created on: 29 June, 2013
 * Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
 * Institute: IMT Lucca, Lucca Italy
 *    __  __      __  
 *   / _ |__) /\ |  \ 
 *   \__)|   /--\|__/ 
 */

#include "mex.h"
#include "../gpad.h"
#include <stdio.h>

/* e = calculate_e(FN, qN, L, C, K, s, isempty_cmin, isempty_cmax, isempty_umin, isempty_umax, isempty_xmin, isempty_xmax, N, w, wN); */
/* MEX Interface for testing only */
/* calculate_e_i is used internally in C in a more efficient way */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {    
    const mxArray *FN_mx = prhs[0], *qN_mx = prhs[1], *L_mx = prhs[2], *C_mx = prhs[3], 
            *K_mx = prhs[4], *s_mx = prhs[5], *isempty_cmin_mx = prhs[6], *isempty_cmax_mx = prhs[7],
            *isempty_umin_mx = prhs[8], *isempty_umax_mx = prhs[9], *isempty_xmin_mx = prhs[10], *isempty_xmax_mx = prhs[11],
            *N_mx = prhs[12], *w_mx = prhs[13], *wN_mx = prhs[14];
    real_t *FN=NULL, *qN=NULL, *L=NULL, *C=NULL, *K=NULL, *s=NULL, *w=NULL,
            *wN=NULL, *Ndbl=NULL, *e, *out = NULL;
    real_t *isempty_cmin=NULL,*isempty_cmax=NULL, *isempty_umin=NULL, *isempty_umax=NULL, 
            *isempty_xmin=NULL, *isempty_xmax=NULL;
    us_t N, nx=0, nu=0, nc=0, nf=0;
    error_t err;
    char error_message_to_matlab[45];
    
    if (nrhs!=15) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e requires exactly 15 input arguments");    
    if (nlhs > 1) mexErrMsgIdAndTxt("MATLAB:calculate_e:badOutput", "calculate_e calculates only one output variable");
    
    FN = mxGetPr(FN_mx);
    qN = mxGetPr(qN_mx);
    L = mxGetPr(L_mx);
    C = mxGetPr(C_mx);
    K = mxGetPr(K_mx);
    s = mxGetPr(s_mx);
    
    isempty_cmin = mxGetPr(isempty_cmin_mx);
    isempty_cmax = mxGetPr(isempty_cmax_mx);
    isempty_xmin = mxGetPr(isempty_xmin_mx);
    isempty_xmax = mxGetPr(isempty_xmax_mx);
    isempty_umin = mxGetPr(isempty_umin_mx);
    isempty_umax = mxGetPr(isempty_umax_mx);
    
    w = mxGetPr(w_mx);
    wN = mxGetPr(wN_mx);
    
    Ndbl = mxGetPr(N_mx);
    
    if (Ndbl==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : N cannot be empty");
    N = (us_t) *Ndbl;
    if (FN==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : FN cannot be empty");
    if (L==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : L cannot be empty");
    /* C may be NULL - No problem! */
    if (K==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : K cannot be empty");
    if (w==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : w cannot be empty");
    if (wN==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : wN cannot be empty");
    if (isempty_cmin==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : isempty_cmin cannot be empty");
    if (isempty_cmax==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : isempty_cmax cannot be empty");
    if (isempty_xmin==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : isempty_xmin cannot be empty");
    if (isempty_xmax==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : isempty_xmax cannot be empty");
    if (isempty_umin==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : isempty_umin cannot be empty");
    if (isempty_umax==NULL) mexErrMsgIdAndTxt("MATLAB:calculate_e:badInput", "calculate_e : isempty_umax cannot be empty");
    
    nx = mxGetM(L_mx);
    nu = mxGetM(K_mx);
    nf = mxGetM(FN_mx); 
    nc = mxGetN(C_mx); /* If C=[], nc=0 */
#if (defined(GPAD_DEBUG) && (GPAD_DEBUG==CHATTERBOX_VERBOSE))
    mexPrintf("Dimensions: \n");
    mexPrintf("nx=%d, nu=%d, nc=%d, nf=%d\n", nx, nu, nc, nf);
#endif
        
    /* Do le magic: */
    e = malloc(nx*(N+1)*DBL_SIZE);
    err = calculate_e(FN, qN, K, L, C, s, w, wN, (cus_t) *isempty_cmin, (cus_t) *isempty_cmax, 
            (cus_t) *isempty_xmin, (cus_t) *isempty_xmax, (cus_t) *isempty_umin, 
            (cus_t) *isempty_umax, nx, nu, nc, nf, N, e);        
    if (err != SUCCESS_OK){
        sprintf(error_message_to_matlab, "Function calculate_e - Error %d.", err );
        mexErrMsgIdAndTxt("MATLAB:calculate_e:failure", error_message_to_matlab);
    }
    
    if (nlhs==0) return; /* Nothing more to do */
    /* Output le result */
    if (nlhs == 1) {
        plhs[0] = mxCreateDoubleMatrix(nx, N+1, mxREAL);
        out = mxGetPr(plhs[0]);
        memcpy(out, e, nx*(N+1)*DBL_SIZE);
    }
    
    
}

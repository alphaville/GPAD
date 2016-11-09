#include "mex.h"
#include "../gpad.h"
#include <stdio.h>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {    
    us_t n;    
    error_t err;
    const mxArray *L_mx = prhs[0], *x_mx = prhs[1];
    real_t *L=NULL, *x=NULL, *y=NULL, *out=NULL;
    char error_message[30];           
    
    L = mxGetPr(L_mx);
    x = mxGetPr(x_mx);
    if (L==NULL) mexErrMsgIdAndTxt("MATLAB:testcholslv:badInput", "testcholslv - L cannot be null");
    if (x==NULL) mexErrMsgIdAndTxt("MATLAB:testcholslv:badInput", "testcholslv - x cannot be null");
    n = (us_t) mxGetM(L_mx);
    
    y=malloc(n * DBL_SIZE);
    err = zchlslv_cpy(L, x, y, n);
    
    if (err!=SUCCESS_OK) {
        sprintf(error_message, "testcholslv - error %d\n", err);
        mexErrMsgIdAndTxt("MATLAB:testcholslv:error", error_message);
        return;
    }
    if (nlhs==1){
        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
        out = mxGetPr(plhs[0]);
        memcpy(out, y, n*DBL_SIZE);
    }
    
}

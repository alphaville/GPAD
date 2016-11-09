/* 
 * File: mx_test_utilities.c
 * Created on: 31 October, 2013
 * Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
 * Institute: IMT Lucca, Lucca Italy
 *    __  __      __  
 *   / _ |__) /\ |  \ 
 *   \__)|   /--\|__/ 
 *
 */
#include "mex.h"
#include "../gpad.h"

#define NORMAL_ZLOWSLV_MODE 0
#define TRANSPOT_ZLOWSLV_MODE 1 

/* z = zlowslv(L, x, mode) */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    /*
     * It solves the lower-triangular system Lz=x
     */
    const mxArray *L_mx = prhs[0], *x_mx = prhs[1], *option_mx = NULL;
    double *L = NULL, *x = NULL, *x_copy = NULL, *out = NULL, *option = NULL;
    int n;
    us_t i, what_to_do = NORMAL_ZLOWSLV_MODE;
    error_t err;
    char error_message[30];
         
    if (nrhs == 3){
        option_mx = prhs[2];
        option = mxGetPr(option_mx);
        what_to_do = (us_t) (*option);
    }

    L = mxGetPr(L_mx);
    x = mxGetPr(x_mx);
    n = mxGetM(x_mx);

    x_copy = malloc(n * DBL_SIZE);
    for (i = 0; i < n; i++) {
        x_copy[i] = x[i];
    }

    if (what_to_do == NORMAL_ZLOWSLV_MODE) {
        err = zlowslv(L, x_copy, n);
    }else {
        err = zlowtrslv(L, x_copy, n);
    }

    if (err != SUCCESS_OK) {
        sprintf(error_message, "testcholslv - error %d\n", err);
        mexErrMsgIdAndTxt("MATLAB:testcholslv:error", error_message);
        return;
    }


    if (nlhs >= 1) {
        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
        out = mxGetPr(plhs[0]);
        memcpy(out, x_copy, n * DBL_SIZE);
    }

}

#include "mex.h"
#include "../gpad.h"

/* g = mx_primal_feasibility(F,G,cmin,cmax,umin,umax,xmin,xmax,x,u) */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    us_t nx, nu, nc, np;
    error_t err;
    const mxArray *F_mx = prhs[0],
            *G_mx = prhs[1],
            *cmin_mx = prhs[2],
            *cmax_mx = prhs[3],
            *umin_mx = prhs[4],
            *umax_mx = prhs[5],
            *xmin_mx = prhs[6],
            *xmax_mx = prhs[7],
            *x_mx = prhs[8],
            *u_mx = prhs[9];
    real_t *F = NULL, *G = NULL, *cmin = NULL, *cmax = NULL, *umin = NULL, *umax = NULL,
            *xmin = NULL, *xmax = NULL, *x = NULL, *u = NULL, *g = NULL, *out = NULL;
    char error_message[45];

    F = mxGetPr(F_mx);
    G = mxGetPr(G_mx);
    cmin = mxGetPr(cmin_mx);
    cmax = mxGetPr(cmax_mx);
    umin = mxGetPr(umin_mx);
    umax = mxGetPr(umax_mx);
    xmax = mxGetPr(xmax_mx);
    xmin = mxGetPr(xmin_mx);
    x = mxGetPr(x_mx);
    u = mxGetPr(u_mx);

    if (x == NULL) mexErrMsgIdAndTxt("MATLAB:primal_feasibility:badInput", "primal_feasibility - x cannot be null");
    if (u == NULL) mexErrMsgIdAndTxt("MATLAB:primal_feasibility:badInput", "primal_feasibility - u cannot be null");


    nx = mxGetM(x_mx);
    nu = mxGetM(u_mx);
    if (F != NULL) {
        nc = mxGetM(F_mx);
    } else {
        nc = 0;
    }
    np = nc * ((cmin != NULL)+(cmax != NULL)) + nx * ((xmin != NULL)+(xmax != NULL)) +
            nu * ((umin != NULL)+(umax != NULL));
    
    g = calloc(np, DBL_SIZE);


    err = primal_feasibility(F, G, cmin, cmax, umin, umax, xmin, xmax, x, u, nx, nu, nc, np, g);
    
    if (err != SUCCESS_OK){
        sprintf(error_message, "Function primal_feasibility - Error %d.", err );
        mexErrMsgIdAndTxt("MATLAB:primal_feasibility:failure", error_message);
    }
    
    if (nlhs==1){
        if (g!=NULL){
        plhs[0] = mxCreateDoubleMatrix(np, 1, mxREAL);
    	out = mxGetPr(plhs[0]);
    	memcpy(out, g, np*DBL_SIZE);
        }else{
            plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        }
    }
    if (nlhs>1){
        mexErrMsgIdAndTxt("MATLAB:primal_feasibility:badOutput", "primal_feasibility exports only one output");
    }

}


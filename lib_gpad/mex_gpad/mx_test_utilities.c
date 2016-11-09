/* 
 * File: mx_test_utilities.c
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


/* MX_TEST_UTILITIES: */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	
	/* Declarations: */
	const mxArray *Q_mx = prhs[0], *x1_mx = prhs[1], *x2_mx = prhs[2], 
	*A_mx = prhs[3], *alpha_mx = prhs[4], *beta_mx = prhs[5], *b_mx = prhs[6];
        double 	*Q=NULL, *x1=NULL, *x2=NULL, *A=NULL, *alpha=NULL, *beta=NULL, *b=NULL, 
                    *alpha_Qx1_pb=NULL, *alpha_Qx1_pb_test=NULL, *alpha_x1tQx2=NULL,
                    *alpha_Qtx2_px1=NULL, *alpha_Qtx2_px1_test=NULL, *zinnerProd=NULL, 
                    *out=NULL, norm_difference, *difference;
        int colsQ, rowsQ;
	error_t err0;
        us_t i;

	/* Read input: */	
	colsQ = mxGetN(Q_mx);
        rowsQ = mxGetM(Q_mx); 
	Q = mxGetPr(Q_mx); x1 = mxGetPr(x1_mx); x2 = mxGetPr(x2_mx); 
	A = mxGetPr(A_mx); alpha = mxGetPr(alpha_mx); beta = mxGetPr(beta_mx), 
	b=mxGetPr(b_mx);
	alpha_Qx1_pb = calloc(rowsQ, DBL_SIZE);
	alpha_Qx1_pb_test = calloc(rowsQ, DBL_SIZE);
	alpha_x1tQx2 = calloc(1, DBL_SIZE);
	alpha_Qtx2_px1 = calloc(colsQ, DBL_SIZE);
	alpha_Qtx2_px1_test = malloc(colsQ*DBL_SIZE);
	zinnerProd = calloc(1, DBL_SIZE);
	
	
	/*
	 * TEST ZQUAD
	 */
	err0 = zquad(x1, Q, x2, rowsQ, colsQ, alpha, 0, NULL);
	if (err0!=NULL_OUTPUT) mexErrMsgIdAndTxt("MATLAB:assertionError", 
		"should have failed (null output)");	
	
	err0 = zquad(0, Q, x2, rowsQ, colsQ, alpha, 0, alpha_x1tQx2);
	if (err0!=NULL_INPUT) mexErrMsgIdAndTxt("MATLAB:assertionError", 
		"should have failed (null input)");

	*alpha_x1tQx2 = *beta;
	err0 = zquad(x1, Q, x2,  rowsQ, colsQ,  alpha, 1, alpha_x1tQx2);
	if (err0!=0) mexErrMsgIdAndTxt("MATLAB:assertionError", "zquad threw an exception");

	/*
	 * TEST ZVMULT
	 *
	 */
	zmvmult(Q, x1, alpha, b, rowsQ, colsQ, 0, alpha_Qx1_pb); /* this works!	*/
	memcpy(alpha_Qx1_pb_test, b, rowsQ * DBL_SIZE);
	zmvmult(Q, x1, alpha, 0, rowsQ, colsQ, 1, alpha_Qx1_pb_test);
	
	difference = malloc(rowsQ * DBL_SIZE);
	for (i=0; i<rowsQ; i++) {
	 	difference[i] = alpha_Qx1_pb_test[i]-alpha_Qx1_pb[i];
	}
	zvnrmsq(difference, rowsQ, NULL, 0, &norm_difference);
	if (norm_difference>1e-8) mexErrMsgIdAndTxt("MATLAB:assertionError", 
		"ZVMULT: norm_difference>1e-8");

	/*
	 * TEST ZVMULT_TRANS
	 */
	zmvmult_trans(Q, x2, alpha, x1, rowsQ, colsQ, 0, alpha_Qtx2_px1);
	*alpha_Qtx2_px1_test = *x1;
	zmvmult_trans(Q, x2, alpha, 0, rowsQ, colsQ, 1, alpha_Qtx2_px1_test);
	difference = malloc(colsQ * DBL_SIZE);
	for (i=0; i<colsQ; i++) {
	 	difference[i] = alpha_Qtx2_px1_test[i]-alpha_Qtx2_px1_test[i];
	}
	zvnrmsq(difference, colsQ, NULL, 0, &norm_difference);
	if (norm_difference>1e-8) mexErrMsgIdAndTxt("MATLAB:assertionError", 
		"ZVMULT_TRANS: norm_difference>1e-8");

	
	/*
	 * TEST ZCDOT
	 */
	 /* Perform the multiplication z = <x2, alpha_Qx1_pb> */
	err0 = zcdot( x2, alpha_Qx1_pb, rowsQ, NULL, 0, zinnerProd);
	if (err0!=0) mexErrMsgIdAndTxt("MATLAB:assertionError", "zdot threw an error");


	/*
	 * OUTPUT THE RESULTS
	 */
	if (nlhs>=1){
    	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    	out = mxGetPr(plhs[0]);
    	memcpy(out, alpha_x1tQx2, DBL_SIZE);
	}
	if (nlhs>=2){
    	plhs[1] = mxCreateDoubleMatrix(rowsQ,1, mxREAL);
    	out = mxGetPr(plhs[1]);
    	memcpy(out, alpha_Qx1_pb, rowsQ * DBL_SIZE);
	}
	if (nlhs>=3){
		plhs[2] = mxCreateDoubleMatrix(colsQ,1, mxREAL);
		out = mxGetPr(plhs[2]);
		memcpy(out, alpha_Qtx2_px1, colsQ * DBL_SIZE);
	}
	if (nlhs>=4){
		plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL);
		out = mxGetPr(plhs[3]);
		memcpy(out, zinnerProd, DBL_SIZE);
	}

	/*
	 * FREE THE ALLOCATED MEMORY
	 */
	if(zinnerProd!=NULL) free(zinnerProd);
	if (alpha_Qx1_pb!=NULL) free(alpha_Qx1_pb);
	if (alpha_Qx1_pb_test!=NULL) free(alpha_Qx1_pb_test);
	if (alpha_x1tQx2!=NULL) free(alpha_x1tQx2);
	if (alpha_Qtx2_px1!=NULL) free(alpha_Qtx2_px1);
	if (difference!=NULL) free(difference);


}
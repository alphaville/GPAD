/* 
 * File: gpad.h
 * Created on: 15 June, 2013
 * Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
 * Institute: IMT Lucca, Lucca Italy
 *    __  __      __  
 *   / _ |__) /\ |  \ 
 *   \__)|   /--\|__/ 
 */

/* Use only PLAIN C - No libraries at all, so that the code is transferable
   to Arduino! */

/*
 * NOTE:
 * All matrices are column-packed, i.e., if double *A is an n-by-m matrix, then 
 * its (i,j)-element (where i=1,...,n and j=1,...,m) is accessed by:
 * 
 * double *element_ij = A + (j-1)*n + (i-1)
 */


#ifndef GPAD_H
#define GPAD_H

/*  Includes:  */
#include <stdlib.h> /* Necessary for malloc and calloc */
#include <string.h> /* Necessary for memcpy and memset */
#include <math.h>   

#ifdef GPAD_DEBUG
#include "mex.h"
#include <stdio.h>
#endif


/* TYPES DEFINITIONS */

/* Consider changing to unsigned int or even unsigned long
   for very big problems. */

/* */
#if !(defined (__GNUG__) && defined (size_t))
#define __size_type__ unsigned short 
#else
#define __size_type__ size_t
#endif

/** Unsigned short type - used for sizes and indices. */
typedef __size_type__ us_t;

/** 
 * The corresponding constant type of us_t; used
 * as an input argument to functions.
 */
typedef const us_t cus_t;
/** Real Number */
typedef double real_t;
/** Constant real_t */
typedef const real_t creal_t;
/** Error Type (Unsigned Short) */
typedef unsigned short error_t;

/*
 * MATRIX TYPES
 */
/** Dense/Full Matrix */
#define MAT_DENSE 0 
/** Matrix in the form alpha*I */
#define MAT_aI 1 
/** Diagonal Matrix */
#define MAT_DIAG 2 
/** Arbitrary Sparse Matrix */
#define MAT_SPARSE 3 

/*
 * SIZES
 */
/** Size of real_t */
#define DBL_SIZE sizeof(real_t)
/** Size of int */
#define INT_SIZE sizeof(int)
/** Size of Unsigned Int */
#define UINT_SIZE sizeof(unsigned int)
/** Size of us_t */
#define USHORT_SIZE sizeof(us_t)
/** A very big double */
#define VERY_BIG_DBL 1e15 

/*
 * ERROR CODES
 *  1XX : User Errors (Bad request)
 *  2XX : Allocation Errors
 *  5XX : Other errors
 */
/** Operation Succeeded */
#define SUCCESS_OK 0
/** Input variable was not allocated (is NULL) */
#define NULL_INPUT 100
/** Output variable was not allocated (is NULL) */
#define NULL_OUTPUT 101
/** Bad input (generic error) */
#define BAD_PARAM_VALUE 102
/** Resource cannot be allocated */
#define ALLOCATION_ERROR 200
/** This is not implemented yet */
#define NOT_IMPLEMENTED_YET 505
/** Negative quantity in a square root */
#define NEGATIVE_IN_SQRT 140
/** Division by zero. */
#define DEVISION_BY_ZERO 141


/** Add the result to the output variable */
#define DO_ADD_OVER 1
/** Copy the result to the output variable */
#define DO_NOT_ADD_OVER 0

/** Absolutely nothing is printed on the output. (Recommended/Faster) */
#define SILENT_MODE 0
/** 
 * Normal Verbose Mode: GPAD to print iteration information. 
 */
#define NORMAL_VERBOSE 1 
/** 
 * To be used only for debugging. (Rather verbose mode). 
 * Stack-trace info printed! 
 */
#define DEBUG_VERBOSE 2
/** Extreme verbose mode - Too much detail printed */
#define CHATTERBOX_VERBOSE 3

#ifdef __INLINE__
#undef __INLINE__
#endif
#ifdef  __cplusplus
extern "C" {
#define __INLINE__ inline
#else
#define __INLINE__
#endif

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
    
#ifndef xmax
#define xmax( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef xmin
#define xmin( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

/* Functions from BLAS - Activated only using the     */
/* preprocessor directive -DBLAS                      */
#if defined(BLAS)
#ifndef CBLAS_H

#ifndef CBLAS_ENUM_DEFINED_H

#define CBLAS_H
#define CBLAS_ENUM_DEFINED_H

    enum CBLAS_ORDER {
        CblasRowMajor = 101,
        CblasColMajor = 102
    };

    enum CBLAS_TRANSPOSE {
        CblasNoTrans = 111,
        CblasTrans = 112,
        CblasConjTrans = 113,
        AtlasConj = 114
    };

    enum CBLAS_UPLO {
        CblasUpper = 121,
        CblasLower = 122
    };

    enum CBLAS_DIAG {
        CblasNonUnit = 131,
        CblasUnit = 132
    };

    enum CBLAS_SIDE {
        CblasLeft = 141,
        CblasRight = 142
    };
#endif /* END: CBLAS_ENUM_DEFINED_H */


    __INLINE__ real_t max(
                       creal_t *array, 
                       cus_t n
                      );
    
    __INLINE__ real_t min(
                       creal_t *array, 
                       cus_t n
                      );
    
    /**
     * Multiplies a matrix by a vector (double precision).
     * This function multiplies A * X (after transposing A, if needed) and multiplies 
     * the resulting matrix by alpha. It then multiplies vector Y by beta. 
     * It stores the sum of these two products in vector Y.
     * 
     * @param Order
     *      Specifies row-major (C) or column-major (Fortran) data ordering.
     * @param TransA
     *      Specifies whether to transpose matrix A.
     * @param M
     *      Number of rows in matrix A.
     * @param N
     *      Number of columns in matrix A.
     * @param alpha
     *      Scaling factor for the product of matrix A and vector X.
     * @param A
     *      Matrix A.
     * @param lda
     *      The size of the first dimension of matrix A; if you are 
     *      passing a matrix A[m][n], the value should be m.
     * @param X
     *      Vector X.
     * @param incX
     *      Stride within X. For example, if incX is 7, every 7th element is used.
     * @param beta
     *      Scaling factor for vector Y.
     * @param Y
     *      Vector Y
     * @param incY
     *      Stride within Y. For example, if incY is 7, every 7th element is used.
     */
    extern void cblas_dgemv(
            const enum CBLAS_ORDER Order,
            const enum CBLAS_TRANSPOSE TransA,
            const int M,
            const int N,
            const double alpha,
            const double *A,
            const int lda,
            const double *X,
            const int incX,
            const double beta,
            double *Y,
            const int incY
            );


    /**
     * Computes the double-precision dot product of a pair of 
     * single-precision vectors.
     * @param N
     *      The number of elements in the vectors.
     * @param X
     *      Vector X.
     * @param incX
     *      Stride within X. For example, if incX is 7, every 7th element is used.
     * @param Y
     *      Vector Y.
     * @param incY
     *      Stride within Y. For example, if incY is 7, every 7th element is used.
     * @return 
     *      The dot product x'y = y'x.
     */
    extern double cblas_ddot(
            const int N,
            const double *X,
            const int incX,
            const double *Y,
            const int incY
            );

    /**
     * Computes a constant times a vector plus a vector (double-precision).
     * On return, the contents of vector Y are replaced with the result. 
     * The value computed is (alpha * X[i]) + Y[i].
     * 
     * @param N
     *      Number of elements in the vectors.
     * @param alpha
     *      Scaling factor for the values in X.
     * @param X
     *      Input vector X.
     * @param incX
     *      Stride within X. For example, if incX is 7, every 7th element is used.
     * @param Y
     *      Input vector Y.
     * @param incY
     *      Stride within Y. For example, if incY is 7, every 7th element is used.
     */
    extern void cblas_daxpy(
            const int N,
            const double alpha,
            const double *X,
            const int incX,
            double *Y,
            const int incY
            );


    extern void cblas_dtrsv (
            const enum CBLAS_ORDER Order,
            const enum CBLAS_UPLO Uplo,
            const enum CBLAS_TRANSPOSE TransA,
            const enum CBLAS_DIAG Diag,
            const int N,
            const double *A,
            const int lda,
            double *X, /* the result */
            const int incX
         );
#endif /* END: CBLAS_H */
#endif /* END: BLAS */

#ifdef GPAD_DEBUG
    /**
     * Displays a matrix or vector in MATLAB using the MATLAB command
     * 'disp'.
     * 
     * @param mRows
     *  Number of rows.
     * @param nCols
     *  Number of dolumns.
     * @param data
     *  The actual data as a pointer to double.
     */
    __INLINE__ void matlabPrint(
            const short mRows,
            const short nCols,
            const double *data
            );
#endif

    /**
     * A dynamical system with state equation: x+ = Ax + Bu + f.
     * 
     */
    typedef struct {
        real_t *A;
        real_t *B;
        real_t *f; 
        real_t *F;
        real_t *G;
        real_t *cmin; 
        real_t *cmax;
        real_t *xmax;
        real_t *xmin;
        real_t *umax;
        real_t *umin;
        us_t nx;
        us_t nu;
        us_t nc;
    } dynsys_t;

    /**
     * A structure for an MPC problem.
     */
    typedef struct {
        real_t  *Q;
        real_t  *R;
        real_t  *S;
        real_t  *q;
        real_t  *r; 
        real_t  *QN; 
        real_t  *qN;
        real_t  *K;
        real_t  *M; 
        real_t  *D;
        real_t  *L;
        real_t  *C;
        real_t  *s;
        real_t  *FN;
        real_t  *gN;
        us_t classQ;
        us_t classR;
        us_t N;
        us_t nf;
        dynsys_t *sys;
    } mpc_t;
    
    /**
     * Diagnostic information for the GPAD solver.
     */
    typedef struct {
        us_t iters;
        real_t  *log_max_viol;
        real_t  *log_max_viol_bar;
        real_t  *dcost;
        real_t  *pcost;
        real_t *duality_gap;
        char *msg;    
        us_t point_1;
        us_t point_2;
        us_t point_3;        
    } diagnostics_t;


    /* Function Declarations...                     */
    /* Auxiliary Functions start with a 'z'.        */
    
    /**
     * Calculates the square root in double precision using the Newton method.
     * @param x
     *          The number for which we need to calculate the square root.
     * @param epsilon
     *          Desired accuracy.
     * @param max_iter
     *          Maximum number of iterations.
     * @return 
     *          The square root of a positive number. This function has indeterminate
     *          behaviour for negative numbers.
     */
    __INLINE__ real_t square_root(
              creal_t x,
              creal_t epsilon,
              cus_t *max_iter
            );

    /**
     * DYNSYS_FACTORY is a factory method for creating dynsys_t variables.
     */
    dynsys_t *dynsys_factory(
              creal_t *A,
              creal_t *B,
              creal_t *f,
              creal_t *F,
              creal_t *G,
              creal_t *cmin,
              creal_t *cmax,
              creal_t *xmin,
              creal_t *xmax,
              creal_t *umin,
              creal_t *umax,
              cus_t nx,
              cus_t nu,
              cus_t nc
            );

    dynsys_t *dynsys_factory_meagre(
             cus_t nx,
             cus_t nu,
             cus_t nc
            );


    mpc_t *mpc_factory(
             creal_t *Q,
             creal_t *R,
             creal_t *S,
             creal_t *q,
             creal_t *r,
             creal_t *QN,
             creal_t *qN,
             creal_t *K,
             creal_t *M,
             creal_t *D,
             creal_t *L,
             creal_t *C,
             creal_t *s,
             creal_t *FN,
             creal_t *gN,
             cus_t classQ,
             cus_t classR,
             cus_t N,
             cus_t nf,
             const dynsys_t *sys
            );


    mpc_t *mpc_factory_simple(
             creal_t *Q,
             creal_t *R,
             creal_t *S,
             creal_t *q,
             creal_t *r,
             creal_t *QN,
             creal_t *qN,
             cus_t classQ,
             cus_t classR,
             cus_t N,
             cus_t nf,
             const dynsys_t *sys
            );


    error_t state_update(
             creal_t *A, 
             creal_t *B, 
             creal_t *f, 
             creal_t *x,
             creal_t *u, 
             cus_t nx, 
             cus_t nu, 
             real_t *x_new
            );
    
    /**
     * STATE_UPDATE calculates the new state of a dynamical system according to
     * the dynamics of the system, i.e., x+ = Ax + Bu + f, or in case f is NULL,
     * then it is x+ = Ax + Bu.
     * 
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t state_update_i(
             dynsys_t *sys,
             creal_t *x,
             creal_t *u,
             real_t *x_new
            );



    /**
     * Solves the triangular system Lx=b where L is a lower triangular 
     * square matrix. This function performs the replacement x = inv(L)*x.
     * @param L
     *          A lower triangular matrix (column packed).
     * @param x
     *          The initial value of x is the vector b. This method updates x
     *          with the solution of the system Lx=b.
     * @param n
     *          The row and column-dimension of L, as well as the dimension
     *          of the vectors x and b.
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t zlowslv(
             creal_t *L,
             real_t *x, 
             cus_t n
            );
    
    /**
     * Solves the triangular system L'*x=b, where L is a lower triangular 
     * matrix.
     * @param L
     *          A lower triangular matrix (column packed).
     * @param x
     *          The initial value of x is the vector b. This method updates x
     *          with the solution of the system L'x=b (L' is the transpose of L).
     * @param n
     *          The row and column-dimension of L, as well as the dimension
     *          of the vectors x and b.
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t zlowtrslv(
             creal_t *L, 
             real_t *x, 
             cus_t n
            );
    
    /**
     * Solves the system Ax=b where A=LL'. A must be a positive definite matrix 
     * and L is its lower-triangular Cholesky factor which must be provided as 
     * input to this function.
     * 
     * @param L
     * @param b
     * @param x
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t zchlslv(
             creal_t *L,
             real_t *x,
             cus_t n
            );
    
    
    /**
     * 
     * @param L
     * @param x
     * @param y
     * @param n
     * @return 
     */
    error_t zchlslv_cpy(
             creal_t *L,
             creal_t *x,
             real_t *y,
             cus_t n
            );
    

    /**
     * Returns the quadratic x'Sx where x is an n-vector and S an n-n symmetric
     * matrix. The result is a scalar of type real_t.
     *      
     * @param x
     *          The vector x.
     * @param S
     *          The square symmetric matrix S.
     * @param n
     *          Length of x.
     * @param L
     *          The result (by reference).
     * 
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t zsym_quad(
             creal_t *x,
             creal_t *S,
             cus_t n,
             cus_t doAddOver,
             real_t *L
            );

    /**
     * ZMVMULT_TRANS performs the algebraic operation `y += alpha*A'*x + b`, where
     * `A` is an n-by-m square matrix, `x` is an m-vector and `b` is an 
     * m-vector and alpha is a scalar. The resulting vector `y` is an m-vector.
     * 
     * @param A
     * @param x
     * @param alpha
     * @param b
     * @param nRows
     * @param nCols
     * @param doAddOver
     * @param y
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t zmvmult_trans(
             creal_t *A,
             creal_t *x,
             creal_t *alpha,
             creal_t *b,
             cus_t nRows,
             cus_t nCols,
             cus_t doAddOver,
             real_t *y
            );
    
    /**
     * ZMVMULT performs the algebraic operation y += alpha*A*x + b, where
     * A is an n-by-m square matrix, x is an m-vector and b is an 
     * n-vector and alpha is a scalar. The resulting vector y is an m-vector.
     * 
     * @param A
     * @param x
     * @param alpha
     * @param b
     * @param nRows
     * @param nCols
     * @param doAddOver
     * @param y
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise. 
     */
    error_t zmvmult(
             creal_t *A,
             creal_t *x,
             creal_t *alpha,
             creal_t *b,
             cus_t nRows,
             cus_t nCols,
             cus_t doAddOver,
             real_t *y
            );
    
    /**
     * 
     * ZQUAD performs the operation z=alpha*x'*Q'*y where x is an nx-vector, 
     * y is an ny-vector and Q is an ny-by-nx matrix.
     * 
     * @param x
     * @param Q
     * @param y
     * @param nx
     * @param ny
     * @param alpha
     * @param doAddOver
     * @param z
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t zquad(
             creal_t *x,
             creal_t *Q,
             creal_t *y,
             cus_t nx,
             cus_t ny,
             creal_t *alpha,
             cus_t doAddOver,
             real_t *z
            );

    /**
     * Calculates the squared Euclidean norm ||x||^2 of a given vector x
     * and updates the variable norm by adding this value, i.e. norm += alpha*||x||^2.
     * It is important that in all cases the pointer `norm` is allocated.
     * The parameter doAddOver, if set to 0, the returned value is actually
     * norm = alpha*||x||^2 - any previous value is discarded. Set `doAddOver=1`
     * for the value of `norm` to be updated by addition, i.e., norm += alpha*||x||^2.
     * 
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t zvnrmsq(
             creal_t *x,
             cus_t n,
             creal_t *alpha,
             cus_t doAddOver,
             real_t *norm
            );

    /**
     * ZCDOT perdorms the operation z=alpha*x'*y where alpha is a scalar, 
     * and x and y are nx-vectors.
     * 
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t zcdot(
             creal_t *x,
             creal_t *y,
             cus_t nx,
             creal_t *alpha,
             cus_t doAddOver,
             real_t *z
            );

    /**
     * Add two vectors. Perform the operation y = y + alpha*x. 
     * @param nx
     *          Size of x and y.
     * @param alpha
     *          Scalar factor.
     * @param X
     *          The first vector (constant).
     * @param Y
     *          The second vector (output).
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t zvadd(
             cus_t nx,
             creal_t alpha,
             creal_t *x,
             real_t *y
            );

    /**
     * Performs the operation : y += alpha*x'*Q*x where Q is a diagonal
     * matrix and x is a vector of dimension nx.
     * The parameter `doAddOver` determines whether the value of `y` should
     * be updated according to y += alpha*x'*Q*x, or the value of `y` should
     * be disregarded, thus setting y = alpha*x'*Q*x.
     * Setting alpha=0 (null), the value *alpha=1.0 will be used.
     * 
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t zvdmquad(
             creal_t *x,
             cus_t nx,
             creal_t *Q,
             creal_t *alpha,
             cus_t doAddOver,
             real_t *y
            );

      
    /**
     * PCOST calculates the primal cost for a given sequence of states X and
     * inputs U over a given prediction horizon N.
     *
     * PCOST returns the primal cost given by:
     * J = V_f(xN) + sum_k [xk' uk']'[Q S'; S R][xk; uk] + [q' r'][xk; uk]
     * where Vf is the terminal cost function given by:
     * Vf(x) = x' QN x + qN'x
     * for a given sequence of states X and a sequence of inputs U.
     *
     * @param Q
     *          The matrix Q mentioned above. If Q has the form
     *          Q = alphaQ*eye(nx) (where nx is the number of states) 
     *          and alphaQ>=0, then we may set classQ=1 and let Q=alphaQ. 
     *          If Q is a diagonal matrix, we may set classQ=2 and let Q be 
     *          a vector with the diagonal elements thereof
     * @param R
     *          The matrix R.
     * @param S
     *          The matrix S.
     * @param q
     *          The vector q.
     * @param r
     *          The vector r.
     * @param QN
     *          The matrix QN.
     * @param qN
     *          The vector qN.
     * @param N
     *          The prediction horizon.
     * @param X
     *          A sequence of states represented as a nx-by-K matrix where nx is 
     *          the dimension of the state vector and K is at least N+1 (The elements of
     *          X in all columns after the N+1-th - if any - will not be used). Normally,
     *          X should be nx-by-(N+1) and x0=X(:,1), x1=X(:,2), and so on.
     * @param U
     *          This is a sequence of input values and is represented as a nu-by-K matrix
     *          where K is at least equal to N.
     * @param nx
     *          The dimension of the state vector.
     * @param nu
     *          The dimension of the input vector.
     * @param classQ
     *          The parameters classQ, classR are integer values determine the 
     *          class of the Matrices Q and R respectively. See MAT_DENSE, MAT_aI,
     *          MAT_DIAG and MAT_SPARSE for more information.
     * @param classR
     *          The class of the matrix R.
     * @param J
     *          The calcualted optimal primal cost.
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
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
            );

    /**
     * Similar method to PCOST, but the matrices of the MPC problem
     * are provided wrapped in a structure.
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t pcost_i(
             const mpc_t *mpc, 
             creal_t *X, 
             creal_t *U, 
             real_t *J
            );


    /**
     * 
     * DCOST is an implementation of Algorithm 4 in [1]. 
     * 
     * @param A The matrices A, B and f are the data of the dynamical 
     *          system: x+ = Ax + Bu + f
     * @param B
     *          The matrix B.
     * @param f
     *          The vector f or NULL if f is empty.
     * @param Q 
     *          The matrices Q, R and S are the parameters of the quadratic 
     *          term of the stage cost.
     * @param R
     *          The positive definite matrix R.
     * @param S
     *          The matrix S or NULL is S is empty.
     * @param q 
     *          The vectors q and r are the linear terms of the stage cost.
     * @param r
     *          The vector r or NULL if r is empty.
     * @param QN 
     *          Matrix of the quadratic term of the terminal cost
     * @param qN 
     *          Vector of the linear term of the terminal cost
     * @param FN 
     *          FN and gN are the matrices of the terminal constraint.
     * @param gN
     *          The vector gN.
     * @param F 
     *          Generic state-input constraints of the form cmin &lt; Fx+Gu &lt; cmax
     * @param G
     *          The matrix G.
     * @param cmin
     *          The vector cmin or NULL if cmin is empty.
     * @param cmax
     *          The vector cmax or NULL if cmax is empty.
     * @param xmin state lower bounds
     *           The vector xmin or NULL if xmin is empty.
     * @param xmax state upper bounds
     *           The vector xmax or NULL if xmax is empty.
     * @param umin input lower bounds
     *           The vector umin or NULL if umin is empty.
     * @param umax input upper bounds
     *           The vector umax or NULL if umax is empty.
     * @param D 
     *          These are matrices calculated by the factorisation algorithm 
     *          offline (for LTI systems).
     * @param L
     *          The matrix L.
     * @param M
     *          The matrix M.
     * @param K
     *          The matrix K.
     * @param C
     *          The matrix C.
     * @param s
     *          The vector s.
     * @param d
     *          The vector d.
     * @param Rchol 
     *          The Cholesky factorisation of Rbar.
     * @param y 
     *          A matrix of dimension np-by-N with the dual vectors that correspond
     *          to the stage-constraints.
     * @param yN 
     *          The dual vector that corresponds to the terminal constraints
     * @param x0 
     *          The current state
     * @param N 
     *          The prediction horizon
     * @param nx 
     *          The dimension of the state vector
     * @param nu 
     *          The dimension of the input vector
     * @param nc 
     *          The row-dimension of F, or 0 if F is empty
     * @param nf 
     *          The number of rows of FN
     * @param dual_cost 
     *          The calculated optimal dual cost
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     * 
     * References:
     * [1] P. Patrinos and A. Bemporad, "An Accelerated Dual
     * Gradient-Projection Algorithm for Embedded Linear Model Predictive
     * Control", 2012, Under review.
     */
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
            );
    
   /**
    * 
    * @param A
    *           The system matrix A.
    * @param B
    *           The system matrix B.
    * @param f
    *           The system vector f or NULL if it does not exist.
    * @param Q
    *           The quadratic cost matrix for the state.
    * @param R
    *           The weighting matrix for the quadratic cost on the inputs.
    * @param S
    *           This is a mixed state-input weighting factor which gives rise
    *           to the term x'Su
    * @param q
    *           A linear cost term for the state, i.e., q'*x.
    * @param r
    *           A linear cost term for the input, i.e., r'u
    * @param QN
    *           Matrix of the quadratic part of the terminal cost. The terminal
    *           cost is given by a function of the form Vf(x)=1/2 x' QN x + qN' x
    * @param qN
    *           The linear term of the terminal cost function; see also
    *           the argumenr QN of this function.
    * @param FN
    *           The matrix FN and the vector qN define the terminal set of the
    *           MPC problem. This is FN xN &lt;= gN
    * @param gN
    *           See FN
    * @param F
    *           The matrices F and G and the vector cmin and cmax define the
    *           double-sided mixed state-input constraints of the form:
    *           cmin &lt;= Fx + Gu &lt;= cmax
    * @param G See F           
    * @param cmin See F
    * @param cmax See F
    * @param xmin Lower bound on the state
    * @param xmax Upper bound on the state
    * @param umin Lower bound on the input
    * @param umax Upper bound on the input
    * @param D 
    * @param L
    * @param M
    * @param K
    * @param C
    * @param s
    * @param d
    * @param Rchol
    * @param y
    * @param yN
    * @param x0
    * @param alpha
    * @param N
    * @param nx
    * @param nu
    * @param nc
    * @param nf
    * @param Xstar
    * @param Ustar
    * @param slack
    * @param slackN
    * @param y_new
    * @param yN_new
    * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
    */
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
             creal_t *y, 
             creal_t *yN, 
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
            );

    /**
     * Similar to dcost, but here all matrices are provided in a structured manner.
     * @param mpc 
     *          An MPC structure of type mpc_t
     * @param Rchol 
     *          The Cholesky factorisation of the matrix Rbar
     * @param y 
     *          The dual vector
     * @param yN 
     *          The dual vector corresponding to the terminal state constraints
     * @param x0 
     *          The state
     * @param dual_cost 
     *          The output of this function, i.e., the value of the optimal
     *          dual cost.
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t dcost_i(
             const mpc_t *mpc,
             creal_t *Rchol,
             creal_t *y,
             creal_t *yN,
             creal_t *x0,
             real_t *dual_cost
            );

    /**
     * Updates theta.
     * @param theta
     * @return This function returns SUCCESS_OK if the computation succeeds or 
     *         an error code otherwise.
     */
    error_t update_theta(
             real_t *theta
            );

    /**
     * 
     * @param FN 
     *          Matrix of the terminal constraints: cmin &lt; FN*xN &lt; cmax
     * @param qN 
     *          Linear terminal cost term
     * @param K 
     *          The matrix K.
     * @param L
     *          The matrix L.
     * @param C
     *          The matrix C.
     * @param s
     *           The vector s.
     * @param w 
     *          The dual vectors w paced in a matrix column-wise.
     * @param wN 
     *          Dual vector for the terminal constraints
     * @param isempty_cmin 
     *          Whether cmin is empty (TRUE/FALSE)
     * @param isempty_cmax 
     *          Whether cmax is empty (TRUE/FALSE)
     * @param isempty_xmin 
     *          Whether xmin is empty (TRUE/FALSE)
     * @param isempty_xmax 
     *          Whether xmax is empty (TRUE/FALSE)
     * @param isempty_umin 
     *          Whether umin is empty (TRUE/FALSE)
     * @param isempty_umax 
     *          Whether umax is empty (TRUE/FALSE)
     * @param nx 
     *          State Dimension
     * @param nu 
     *          Input Dimension
     * @param nc 
     *          Number of rows of the matrix F, as in Fx+Gu 
     * @param nf 
     *          Number of terminal constraints
     * @param N 
     *          Prediction horizon
     * @param e 
     *          The resulting matrix e, of dimension nx times N+1.
     * @return  This function returns SUCCESS_OK if the computation succeeds or 
     *          an error code otherwise.
     */
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
            );
    
    /**
     * Calculates the value of the function g(x,u) of the primal constraints:
     * g(x,u)&lt;0.
     * 
     * @param F
     *          The matrix F as in cmin &lt; Fx+Gu &lt; cmax
     * @param G
     *          The matrix G as above of dimensions nc-by-nu.
     * @param cmin
     *          The vector cmin of dimension nc.
     * @param cmax
     *          The vector cmax.
     * @param umin
     *          The lower bound on the input variable.
     * @param umax
     *          The upper bound on u.
     * @param xmin
     *          The lower bound on x.
     * @param xmax
     *          The upper bound on x.
     * @param x
     *          The state for which g(x,u) will be calculated.
     * @param u
     *          The corresponding input vector.
     * @param nx
     *          Dimension of the state vector.
     * @param nu
     *          Dimension of the input vector.
     * @param nc
     *          Row-count of F and G, or 0 if there are no mixed state-input
     *          constraints.
     * @param g
     *          The result.
     * @return  This function returns SUCCESS_OK if the computation succeeds or 
     *          an error code otherwise.
     */
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
            );

    
    /**
     * 
     * @param A
     * @param B
     * @param f
     * @param Q
     * @param R
     * @param S
     * @param q
     * @param r
     * @param QN
     * @param qN
     * @param FN
     * @param gN
     * @param F
     * @param G
     * @param cmin
     * @param cmax
     * @param xmin
     * @param xmax
     * @param umin
     * @param umax
     * @param D
     * @param L
     * @param M
     * @param K
     * @param C
     * @param s
     * @param d
     * @param Rchol
     * @param x0
     * @param Xstar
     * @param Ustar
     * @param nx
     * @param nu
     * @param nc
     * @param classQ
     * @param classR
     * @param enforce_monotonicity
     * @return 
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
          );


#ifdef __cplusplus
}
#endif


#endif /* GPAD_H */
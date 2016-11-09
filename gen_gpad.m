N=10;
cmax=g;
cmin=-g;
np = nc + 2*(nx+nu);
FN = FN + rand(size(FN))*0.01;

filename = './lib_gpad/gpad_data';
fid = fopen([filename '.h'], 'w', 'n', 'UTF-8');
if fid==-1
    error('Cannot open file. Mode: w');
end
fprintf(fid, ['/* \n * File: ' filename '.h\n']);
fprintf(fid, ' * Automatically Created on: 4 July, 2013\n');
fprintf(fid, ' * using GPADToolbox\n');
fprintf(fid, ' * Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>\n');
fprintf(fid, ' * Institute: IMT Lucca, Lucca Italy\n');
fprintf(fid, ' *    __  __      __  \n');
fprintf(fid, ' *   / _ |__) /\\ |  \\ \n');
fprintf(fid, ' *   \\__)|   /--\\|__/ \n');
fprintf(fid, ' *\n */\n\n\n');
fprintf(fid, '#ifndef GPAD_DATA_H\n#define GPAD_DATA_H\n\n');

MAT_CLASS_DENSE = 0;
MAT_CLASS_aI = 1;
MAT_CLASS_DIAGONAL = 2;
fprintf(fid, '/* MAT_CLASS_DENSE: A Matrix is dense */\n');
fprintf(fid, '#define MAT_CLASS_DENSE 0\n');
fprintf(fid, '/* MAT_CLASS_aI: A Matrix has the form alpha*I, where I is the identity matrix */\n');
fprintf(fid, '#define MAT_CLASS_aI 1\n');
fprintf(fid, '/* MAT_CLASS_DENSE: A Matrix is diagonal */\n');
fprintf(fid, '#define MAT_CLASS_DIAGONAL 2\n');
fprintf(fid, '/* MAT_CLASS_DENSE: A Matrix has some other sparsity pattern */\n');
fprintf(fid, '#define MAT_CLASS_SPARSE 3\n\n');

% Define a function `isdiag` which returns 1 if a matrix is diagonal 
% and 0 otherwise
isdiag=@(A) ~nnz(A-spdiags(diag(A),0,speye(size(A,1))));
% Explore the structure of Q
fprintf(fid, '/* CLASS_Q: The class/type of the matrix Q */\n');
if isdiag(Q),
    if (all(diag(Q)==Q(1))),        
        fprintf(fid, '#define CLASS_Q MAT_CLASS_aI\n');
        class_Q = MAT_CLASS_aI;
    else
        fprintf(fid, '#define CLASS_Q MAT_CLASS_DIAGONAL\n');
        class_Q = MAT_CLASS_DIAGONAL;
    end
else
    fprintf(fid, '#define CLASS_Q MAT_CLASS_DENSE\n');
    class_Q = MAT_CLASS_DENSE;
end

fprintf(fid, '/* CLASS_Q: The class/type of the matrix R */\n');
if isdiag(R),
    if (all(diag(R)==R(1))),
        fprintf(fid, '#define CLASS_R MAT_CLASS_aI\n');
        class_R = MAT_CLASS_aI;
    else
        fprintf(fid, '#define CLASS_R MAT_CLASS_DIAGONAL\n');
        class_R = MAT_CLASS_DIAGONAL;
    end
else
    fprintf(fid, '#define CLASS_R MAT_CLASS_DENSE\n');
    class_R = MAT_CLASS_DENSE;
end

fprintf(fid, '\n');
fprintf(fid, '#ifdef FALSE\n#undef FALSE\n#endif\n#define FALSE 0\n');
fprintf(fid, '#ifdef TRUE\n#undef TRUE\n#endif\n#define TRUE 1\n');
fprintf(fid, '#define IS_EMPTY_CMIN FALSE\n');
fprintf(fid, '#define IS_EMPTY_CMAX FALSE\n');
fprintf(fid, '#define IS_EMPTY_XMIN FALSE\n');
fprintf(fid, '#define IS_EMPTY_XMAX FALSE\n');
fprintf(fid, '#define IS_EMPTY_UMIN FALSE\n');
fprintf(fid, '#define IS_EMPTY_UMAX FALSE\n');
fprintf(fid, '#define IS_EMPTY_qN TRUE\n');
fprintf(fid, '\n');

fprintf(fid, '\n/* Typedefs */\n');
fprintf(fid, '#if !(defined (__GNUG__) && defined (size_t))\n');
fprintf(fid, '#define gpad_size_type unsigned short\n');
fprintf(fid, '#else\n#define gpad_size_type size_t\n#endif\n\n');
fprintf(fid, 'typedef gpad_size_type us_t;\ntypedef const us_t cus_t;\n');
fprintf(fid, 'typedef double real_t;\n');
fprintf(fid, 'typedef const real_t creal_t;\n');

fprintf(fid,'\n/* Functions */\n');
fprintf(fid,'void calculate_e(void);\n');
fprintf(fid,'void zvadd(cus_t n, creal_t alpha, creal_t *x, real_t *y);\n');
fprintf(fid,'void zmvmult(creal_t a[], creal_t x[], creal_t alpha, creal_t b[], cus_t nRows, cus_t nCols, real_t y[]);\n');
fprintf(fid,'void zmvmult_trans(creal_t a[], creal_t x[], creal_t alpha, creal_t b[], cus_t nRows, cus_t nCols, real_t y[]);\n');
fprintf(fid,'void zlowslv(creal_t low_triang[], real_t x[], cus_t n);\n');
fprintf(fid,'void zlowtrslv(creal_t low_triang[], real_t x[], cus_t n);\n');
fprintf(fid,'real_t zvdmquad(creal_t x[], cus_t nx, creal_t q[]);\n');
fprintf(fid,'real_t zvnrmsq(creal_t x[], cus_t n, creal_t alpha);\n');
fprintf(fid,'real_t zcdot(creal_t x[], creal_t y[], cus_t nx, creal_t alpha);\n');
fprintf(fid,'real_t zsym_quad(creal_t x[], creal_t s[], cus_t n);\n');
fprintf(fid,'real_t zquad(creal_t x[], creal_t q[], creal_t y[], cus_t nx, cus_t ny, creal_t alpha);\n');
fprintf(fid,'real_t square_root(creal_t x, creal_t epsilon, cus_t max_iter);\n');


nxNp1 = nx*(N+1);
nuN = nu*N;
nxnu = nx*nu;
nxnx = nx*nx;
npN = np*N;
fprintf(fid, '\n/* Sizes */\n');
fprintf(fid, 'cus_t NX_ = %d;\n', nx);
fprintf(fid, 'cus_t NU_ = %d;\n', nu);
fprintf(fid, 'cus_t NC_ = %d;\n', nc);
fprintf(fid, 'cus_t NF_ = %d;\n', nf);
fprintf(fid, 'cus_t NP_ = %d;\n', np);
fprintf(fid, 'cus_t NX_NX_ = %d;\n', nxnx);
fprintf(fid, 'cus_t NX_NU_ = %d;\n', nxnu);
fprintf(fid, 'cus_t N_ = %d;\n', N);
fprintf(fid, 'cus_t NP_N_ = %d;\n',  npN);
fprintf(fid, 'cus_t NP_Np1_ = %d;\n',  nxNp1);
fprintf(fid, 'cus_t NX_N_ = %d;\n',  nx*N);
fprintf(fid, 'cus_t NX_Np1_ = %d;\n',  nx*N+1);
fprintf(fid, 'cus_t NU_N_ = %d;\n',  nuN);

fprintf(fid, '\n/* Data - Scalars */\n');
fprintf(fid, 'creal_t ALPHA_ = %g;\n', alpha);
fprintf(fid, 'creal_t EPSILON_V_ = %g;\n', e_V);
fprintf(fid, 'creal_t EPSILON_g_ = %g;\n', e_g);
fprintf(fid, 'creal_t MAX_ITER_ = %d;\n', max_iter);

fprintf(fid, '\n/* Pre-Allocated Memory */\n');
fprintf(fid, 'real_t gTHETA = 1.0, gTHETA_P = 1.0, gTHETA2 = 0.0, gBETA = 0.0;\n');
gen_print_mat(fid, zeros(npN,1),'gY',[],'real_t');
fprintf(fid, 'real_t gY_PREV[%d];\n', npN);
fprintf(fid, 'real_t gE[%d];\n', nxNp1);

gen_print_mat(fid, zeros(nc,1),'gFXGU',[],'real_t'); % Only if cmin!=NULL || cmax!=NULL
gen_print_mat(fid, zeros(npN,1),'gW',[],'real_t');
gen_print_mat(fid, zeros(nf,1),'gWN',[],'real_t');
gen_print_mat(fid, zeros(nf,1),'gYN',[],'real_t');
gen_print_mat(fid, zeros(nf,1),'gYN_PREV',[],'real_t');
gen_print_mat(fid, zeros(npN,1),'gSLACK',[],'real_t');
gen_print_mat(fid, zeros(nf,1),'gSLACK_N',[],'real_t');
gen_print_mat(fid, zeros(npN,1),'gSLACK_BAR',[],'real_t');
gen_print_mat(fid, zeros(nf,1),'gSLACK_N_BAR',[],'real_t');
fprintf(fid, 'real_t gPR_FEAS[%d];\n', np);

fprintf(fid, 'real_t gXstar[%d];\n', nxNp1);
fprintf(fid, 'real_t gUstar[%d];\n', nuN);

fprintf(fid, '\n\n/* Data - Matrices */\n');
gen_print_mat(fid, x0,'X0_',[],'real_t');


if class_Q == MAT_CLASS_DENSE,
    gen_print_mat(fid, Q,'Q_');
elseif class_Q == MAT_CLASS_aI,
    fprintf(fid, 'creal_t Q_ = %g;\n', Q(1));
elseif class_Q == MAT_CLASS_DIAGONAL,
    gen_print_mat(fid, diag(Q),'Q_');
end

if class_R == MAT_CLASS_DENSE,
    gen_print_mat(fid, R,'R_');
elseif class_R == MAT_CLASS_aI,
    fprintf(fid, 'creal_t R_ = %g;\n', R(1));
elseif class_R == MAT_CLASS_DIAGONAL,
    gen_print_mat(fid, diag(R),'R_');
end

gen_print_mat(fid, Pf,'QN_');
gen_print_mat(fid, A,'A_');
gen_print_mat(fid, B,'B_');
gen_print_mat(fid, F,'F_');
gen_print_mat(fid, G,'G_');
gen_print_mat(fid, FN,'FN_');
gen_print_mat(fid, gN,'GN_');
gen_print_mat(fid, xmin,'XMIN_');
gen_print_mat(fid, xmax,'XMAX_');
gen_print_mat(fid, umin,'UMIN_');
gen_print_mat(fid, umax,'UMAX_');
gen_print_mat(fid, cmax,'CMAX_');
gen_print_mat(fid, cmin,'CMIN_');
gen_print_mat(fid, K,'K_');
gen_print_mat(fid, L,'L_');
gen_print_mat(fid, M,'M_');
gen_print_mat(fid, mpc.D,'D_');
gen_print_mat(fid, mpc.C,'C_');
gen_print_lowtr(fid, mpc.RbarChol,'R_BAR_CHOL_');



fprintf(fid, '\n\n\n#endif /* GPAD_DATA_H */\n');
fclose(fid);


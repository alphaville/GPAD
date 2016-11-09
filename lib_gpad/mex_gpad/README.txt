GPAD-Toolbox (MEX Interfaces)

An implementation of GPAD in ANSI-C with interfaces for MATLAB.
You can either use this C-library in your C code or invoke it from MATLAB
using the MEX interfaces.
 
  Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
  Institute: IMT Lucca, Lucca Italy
     __  __      __  
    / _ |__) /\ |  \ 
    \__)|   /--\|__/ 


Contents:
mx_pcost.c          - MATLAB Interface method (MEX) to calculate the primal cost
mx_dcost.c          - MATLAB Interface method to calculate the optimal dual cost
mx_state_update.c   - Calcuate the next state of an LTI system
mx_calculate_e.c    - Calcuate the matrix e as in Algorithm 4
mx_test_utilities.c - A test interface file
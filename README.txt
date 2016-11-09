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
./lib_gpad    - The ANSI-C implementation of GPAD and the MEX interfaces
./test        - Some unit test files
make.m        - A MATLAB Make-file that allows you to compile the MEX interfaces
./*.m         - MATLAB files with documentation (for built-in functions)
checksums     - SHA-Checksums of all files in the project
verifysums.sh - Shell script to verify the checksums


Additionally, the following classes are available to facilitate the
use of the toolbox:
DynamicalSystem   - A general dynamical system (abstract).
LTISystem         - Define your own LTI system of the form x+=Ax+Bu+f and
                    the accompanying constraints.
MPCGPAD           - An MPC controller powered by GPAD.


Requirements:
CBLAS (optional)  - CBLAS is an optional dependency. You may specify the
                    option `no_blas` during compilation if you need to compile
                    without BLAS. Otherwise, specify the path to your CBLAS
                    installation in make.m.
                    
                    
Please cite this work as follows:
P. Patrinos, A. Bemporad, "An accelerated dual gradient-projection algorithm
for embedded linear model predictive control," IEEE Transactions on Automatic 
Control, 59 (1) (2014), pp. 18–33.

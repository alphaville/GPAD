function make(varargin)
%MAKE is a make-function for easy compilation of the GPAD-Toolbox.
%
%Syntax:
% >> make [[target(s)]] [[option(s)]]
%
%Targets:
% * clean
% * all 
% * pcost
% * dcost
% * testutil
% * state_update
% * calculate_e
% * primal_feasibility
% * zlowslv
% 
%Options:
% + debug (debug mode)
% + verbose (print compilation info)
% + no_blas (do not use BLAS)
%
%Examples of use:
% >> make clean all verbose
% >> make dcost debug
%

all = false;
clean = false;
pcost=false;
dcost=false;
testutil=false;
state_update=false;
calculate_e=false;
primal_feasibility=false;
dgrad=false;
gpad=false;
zlowslv=false;

% Modify these constants if necessary:
lib_blas_path = '/usr/lib/libcblas.dylib'; % PATH of CBLAS
my_c_flags = '-Wall -O4 -pedantic'; % C Flags
lib_gpad_path='./lib_gpad/'; % Where the C files are

is_verbose=false;       %Whether to print details
do_use_blas = true;     %Whether to use blas
is_debug_mode = false;  %What to use debug mode (the flag -DGPAD_DEBUG will
                        %be used in compilation in this case)


for i=1:length(varargin)
    varargin_i = varargin{i};
    if strcmp(varargin_i,'all'), all = true; 
    elseif strcmp(varargin_i,'clean'), clean = true; 
    elseif strcmp(varargin_i,'pcost'), pcost = true; 
    elseif strcmp(varargin_i,'dcost'), dcost = true;
    elseif strcmp(varargin_i,'dgrad'), dgrad = true;
    elseif strcmp(varargin_i,'gpad'), gpad = true;
    elseif strcmp(varargin_i,'zlowslv'), zlowslv = true;        
    elseif strcmp(varargin_i,'testutil'), testutil = true; 
    elseif strcmp(varargin_i,'state_update'), state_update = true; 
    elseif strcmp(varargin_i,'calculate_e'), calculate_e = true; 
    elseif strcmp(varargin_i,'verbose'), is_verbose = true; 
    elseif strcmp(varargin_i,'no_blas'), do_use_blas = false; 
    elseif strcmp(varargin_i,'debug'), disp('DEBUG=ON'); is_debug_mode = true; 
    elseif strcmp(varargin_i,'primal_feasibility'), primal_feasibility = true; 
        
    else error(['Unrecognised target/option : ' varargin{i}]);
    end
end
if (~all&&~pcost&&~testutil&&~state_update&&~calculate_e&&~dcost...
        &&~primal_feasibility&&~dgrad&&~gpad&&~zlowslv), all = true; end
if clean
    ! rm -rf *~
    ! rm -rf *.mexmaci64
    ! rm -rf *.o
    disp('Cleaning');
end

if all || isempty(varargin)
    compile_pcost();
    compile_dcost();
    compile_testutil();
    compile_state_update();
    compile_calculate_e();
    compile_primal_feasibility();
    compile_dgrad();
    compile_gpad();
    compile_zlowslv();
end

if (~all && pcost), compile_pcost(); end
if (~all && zlowslv), compile_zlowslv(); end
if (~all && dcost), compile_dcost(); end
if (~all && testutil), compile_testutil(); end
if (~all && state_update), compile_state_update();end
if (~all && calculate_e), compile_calculate_e();end
if (~all && dgrad), compile_dgrad();end
if (~all && gpad), compile_gpad();end
if (~all && primal_feasibility), compile_primal_feasibility();end

    function compile_zlowslv()
        disp('Compiling  : ZLOWSLV');
        mex_command = ['mex CFLAGS=''$CFLAGS ' my_c_flags];
        if (is_debug_mode), mex_command = [mex_command ' -DGPAD_DEBUG']; end;
        if (do_use_blas), mex_command = [mex_command ' -DBLAS']; end;
        mex_command = [mex_command ''' -largeArrayDims ' ...
            '-g -O -output zlowslv ' lib_gpad_path 'gpad.c ' ...
            lib_gpad_path 'mex_gpad/mx_zlowslv.c '];
        if (do_use_blas), mex_command = [mex_command lib_blas_path]; end
        if (is_verbose), disp([' -- ' mex_command]); disp(' '); end
        eval(mex_command);
    end

    function compile_pcost()
        disp('Compiling  : PCOST');
        mex_command = ['mex CFLAGS=''$CFLAGS ' my_c_flags];
        if (is_debug_mode), mex_command = [mex_command ' -DGPAD_DEBUG']; end;
        if (do_use_blas), mex_command = [mex_command ' -DBLAS']; end;
        mex_command = [mex_command ''' -largeArrayDims ' ...
            '-g -O -output pcost ' lib_gpad_path 'gpad.c ' ...
            lib_gpad_path 'mex_gpad/mx_pcost.c '];
        if (do_use_blas), mex_command = [mex_command lib_blas_path]; end
        if (is_verbose), disp([' -- ' mex_command]); disp(' '); end
        eval(mex_command);
    end

    function compile_dcost()
        disp('Compiling  : DCOST');
        mex_command = ['mex CFLAGS=''$CFLAGS ' my_c_flags];
        if (is_debug_mode), mex_command = [mex_command ' -DGPAD_DEBUG']; end;
        if (do_use_blas), mex_command = [mex_command ' -DBLAS']; end;
        mex_command = [mex_command ''' -largeArrayDims ' ...
            '-g -O -output dcost ' lib_gpad_path 'gpad.c ' ...
            lib_gpad_path 'mex_gpad/mx_dcost.c '];
        if (do_use_blas), mex_command = [mex_command lib_blas_path]; end
        if (is_verbose), disp([' -- ' mex_command]); disp(' '); end
        eval(mex_command);
    end

    function compile_dgrad()
        disp('Compiling  : DGRAD');
        mex_command = ['mex CFLAGS=''$CFLAGS ' my_c_flags];
        if (is_debug_mode), mex_command = [mex_command ' -DGPAD_DEBUG']; end;
        if (do_use_blas), mex_command = [mex_command ' -DBLAS']; end;
        mex_command = [mex_command ''' -largeArrayDims ' ...
            '-g -O -output dgrad ' lib_gpad_path 'gpad.c ' ...
            lib_gpad_path 'mex_gpad/mx_dgrad.c '];
        if (do_use_blas), mex_command = [mex_command lib_blas_path]; end
        if (is_verbose), disp([' -- ' mex_command]); disp(' '); end
        eval(mex_command);
    end

    function compile_gpad()
        disp('Compiling  : GPAD!');
        mex_command = ['mex CFLAGS=''$CFLAGS ' my_c_flags];
        if (is_debug_mode), mex_command = [mex_command ' -DGPAD_DEBUG=2']; end;
        if (do_use_blas), mex_command = [mex_command ' -DBLAS']; end;
        mex_command = [mex_command ''' -largeArrayDims ' ...
            '-g -O -output gpad ' lib_gpad_path 'gpad.c ' ...
            lib_gpad_path 'mex_gpad/mx_gpad.c '];
        if (do_use_blas), mex_command = [mex_command lib_blas_path]; end
        if (is_verbose), disp([' -- ' mex_command]); disp(' '); end
        eval(mex_command);
    end

    function compile_testutil()
        disp('Compiling  : TESTUTIL');
        mex_command = ['mex CFLAGS=''$CFLAGS ' my_c_flags];      
        if (is_debug_mode), mex_command = [mex_command ' -DGPAD_DEBUG']; end;
        if (do_use_blas), mex_command = [mex_command ' -DBLAS']; end;
        mex_command = [mex_command ''' -largeArrayDims ' ...
            '-g -O -output testutil ' lib_gpad_path 'gpad.c ' ...
            lib_gpad_path 'mex_gpad/mx_test_utilities.c '];
        if (do_use_blas), mex_command = [mex_command lib_blas_path]; end
        if (is_verbose), disp([' -- ' mex_command]); disp(' '); end
        eval(mex_command);
    end

    function compile_state_update()
        disp('Compiling  : STATE_UPDATE');
        mex_command = ['mex CFLAGS=''$CFLAGS ' my_c_flags];    
        if (is_debug_mode), mex_command = [mex_command ' -DGPAD_DEBUG=2']; end;
        if (do_use_blas), mex_command = [mex_command ' -DBLAS']; end;
        mex_command = [mex_command ''' -largeArrayDims ' ...
            '-g -O -output state_update ' lib_gpad_path 'gpad.c ' ...
            lib_gpad_path 'mex_gpad/mx_state_update.c '];
        if (do_use_blas), mex_command = [mex_command lib_blas_path]; end
        if (is_verbose), disp([' -- ' mex_command]); disp(' '); end
        eval(mex_command);
    end

    function compile_calculate_e()
        disp('Compiling  : COMPILE_E');
        mex_command = ['mex CFLAGS=''$CFLAGS ' my_c_flags];    
        if (is_debug_mode), mex_command = [mex_command ' -DGPAD_DEBUG=2']; end;
        if (do_use_blas), mex_command = [mex_command ' -DBLAS']; end;
        mex_command = [mex_command ''' -largeArrayDims ' ...
            '-g -O -output calculate_e ' lib_gpad_path 'gpad.c ' ...
            lib_gpad_path 'mex_gpad/mx_calculate_e.c '];
        if (do_use_blas), mex_command = [mex_command lib_blas_path]; end
        if (is_verbose), disp([' -- ' mex_command]); disp(' '); end
        eval(mex_command);
    end

    function compile_primal_feasibility()
        disp('Compiling  : PRIMAL_FEASIBILITY');
        mex_command = ['mex CFLAGS=''$CFLAGS ' my_c_flags];    
        if (is_debug_mode), mex_command = [mex_command ' -DGPAD_DEBUG=2']; end;
        if (do_use_blas), mex_command = [mex_command ' -DBLAS']; end;
        mex_command = [mex_command ''' -largeArrayDims ' ...
            '-g -O -output primal_feasibility ' lib_gpad_path 'gpad.c ' ...
            lib_gpad_path 'mex_gpad/mx_primal_feasibility.c '];
        if (do_use_blas), mex_command = [mex_command lib_blas_path]; end
        if (is_verbose), disp([' -- ' mex_command]); disp(' '); end
        eval(mex_command);
    end

end
% Author: Pantelis Sopasakis
% Created on: 15 June, 2013
% Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
% Institute: IMT Lucca, Lucca Italy
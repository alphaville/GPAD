clc
disp('COMPILATION');
make all no_blas

disp(' ');

disp('TEST');
test_factor;
test_dgrad;
test_dcost;
test_pcost;
test_utilities;
test_state_update;
test_00;
test_calculate_e;
test_cholslv;
test_dual;
test_prfeas;
test_gpad;

disp(' ');
disp('ALL TESTS PASSED!');
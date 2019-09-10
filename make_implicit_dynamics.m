if isunix
    mex GCC=/usr/bin/g++ -output implicit_dynamics matlab_interface.c implicitDynamics.c lie_theory.c -lgsl -lgslcblas -lm -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
    mex GCC=/usr/bin/g++ -output implicit_dynamics_step matlab_interface_step.c implicitDynamics.c lie_theory.c -lgsl -lgslcblas -lm -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
elseif ispc
    setenv('MW_MINGW64_LOC','C:\TDM-GCC-64');
    mex -setup cpp
    mex -output implicit_dynamics -I'C:\MinGW\msys\1.0\local\include\' matlab_interface.c implicitDynamics.c lie_theory.c 'C:\MinGW\msys\1.0\local\lib\libgsl.a' 'C:\MinGW\msys\1.0\local\lib\libgslcblas.a' -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
    mex -output implicit_dynamics_step -I'C:\MinGW\msys\1.0\local\include\' matlab_interface_step.c implicitDynamics.c lie_theory.c 'C:\MinGW\msys\1.0\local\lib\libgsl.a' 'C:\MinGW\msys\1.0\local\lib\libgslcblas.a' -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
else
    error(['Haven''t tested compilation on this platform. Try something like:', newline, 'mex -output implicit_dynamics matlab_interface.c implicitDynamics.c -lgsl -lgslcblas -lm']);
end
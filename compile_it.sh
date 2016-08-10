#!/bin/bash -e
echo "bash: compile"
gfortran -I/usr/include \
const_mod.f90 \
sys_state_mod.f90 \
plans_mod.f90 \
trafo_mod.f90 \
nabla_mod.f90 \
init_mod.f90 \
exit_mod.f90 \
IO_mod.f90 \
pdgl_mod.f90 \
timestepping_mod.f90 \
test_mod.f90 \
main.f90 \
-lfftw3 \
-lm \
-O3 \
-g
echo "bash: compilation...	 done"

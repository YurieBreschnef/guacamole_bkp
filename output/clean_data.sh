#!/bin/bash

echo "bash: clean_data.."

cd data/u/
yes | rm *.dat
cd ..

cd u_f/
yes | rm *.dat
cd ..

cd abs_u/
yes | rm *.dat
cd ..

cd buo/
yes | rm *.dat
cd ..

cd temp/
yes | rm *.dat
cd ..

cd temp_f/
yes | rm *.dat
cd ..

cd chem/
yes | rm *.dat
cd ..

cd chem_f/
yes | rm *.dat
cd ..

cd vort/
yes | rm *.dat
cd ..

cd div/
yes | rm *.dat
cd ..

cd u_stat/
yes | rm *.dat
cd ..
cd E_stat/
yes | rm *.dat
cd ..

cd C_stat/
yes | rm *.dat
cd ..

cd sys_stat/
yes | rm *.dat
cd ..

cd T_stat/
yes | rm *.dat
cd ..

echo "bash: clean_data..	done"


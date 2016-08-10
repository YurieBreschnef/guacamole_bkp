#!/bin/bash

cd output/
bash clean_data.sh
bash clean_visual.sh
cd ..


bash compile_it.sh

echo "bash: starting run"
export OMP_NUM_THREADS=4
./a.out
echo "bash: run... done"

echo "bash: plotting"
cd output/
gnuplot plotall_2D.gp
echo "bash: plotting...	 done"

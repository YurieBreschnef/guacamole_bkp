#!/bin/bash

echo "bash: clean_visual.."

cd visual/u/
yes | rm *.png
cd ..      

cd abs_u/    
yes | rm *.png
cd ..      
          
cd u_f/    
yes | rm *.png
cd ..      
           
cd temp/   
yes | rm *.png
cd ..     
          
cd temp_f/
yes | rm *.png
cd ..      
           
cd chem/   
yes | rm *.png
cd ..      
           
cd chem_f/ 
yes | rm *.png
cd ..

cd stat/
yes | rm *.png
cd ..

cd buo/
yes | rm *.png
cd ..

cd div/
yes | rm *.png
cd ..

cd combo/
yes | rm *.png
cd ..

cd vort/
yes | rm *.png
cd ..

cd temp_chem_combo/
yes | rm *.png
cd ..

echo "bash: clean_visual.. done"


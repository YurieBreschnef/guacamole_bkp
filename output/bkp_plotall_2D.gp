load './gnuplot-palettes-master/jet.pal'    
aspect_ratio = 2
no_of_img = 299 

set terminal pngcairo size 800,800 enhanced font 'Verdana,10'

    set output './visual/statistics/Ekin_vs_Epot.png'
		set multiplot layout 2,2
    set title 'simulation time vs. mean kinetic energy'
    plot './data/statistics/measures.dat' using 2:8 title "E_{pot}" #with lines 
    set title 'simulation time vs. potential ernergy'
    plot './data/statistics/measures.dat' using 2:3 title "E_{kin}" #with lines 
    set title 'simulation time vs. maximum velocity (all directions)'
    plot './data/statistics/measures.dat' using 2:9 title "v_{max}" #with lines 
    set title 'simulation time vs. dt'
    plot './data/statistics/measures.dat' using 2:10 title "dt" #with lines 
    unset multiplot


    set output './visual/statistics/Temp.png'
    set title 'simulation time vs. mean Temperature'
    plot './data/statistics/measures.dat' using 2:4 title "mean temp" with lines  , \
         './data/statistics/measures.dat' using 2:5 title "max temp" with lines  

    set output './visual/statistics/Chem.png'
    set title 'simulation time vs. mean chemical concentration'
    plot './data/statistics/measures.dat' using 2:6 title "mean chem" with lines  , \
         './data/statistics/measures.dat' using 2:7 title "max chem" with lines  

do for [i=0:no_of_img] {
    set size ratio aspect_ratio 

    set terminal pngcairo size 800,800 enhanced font 'Verdana,10'

  	set output './visual/temp/'.i.'.png'
  	set title 'temperature field'
      	plot './data/temp/'.i.'.temp.dat' using 1:2:3  with image

  	set output './visual/u/'.i.'.png'
  	set title 'velocity field'
      	plot './data/u/'.i.'.u.dat' using 1:2:($3*50):($4*50)  with vectors

  	set output './visual/abs_u/'.i.'.png'
  	set title 'absolute magnitude  ofvelocity field'
      	plot './data/abs_u/'.i.'.abs_u.dat' using 1:2:3 with image

  	set output './visual/u_f/'.i.'.png'
   unset size
   set title 'absolute of fourier spectrum'
      	plot './data/u_f/'.i.'.u_f.dat' using 1:2:3  with image

    set output './visual/chem/'.i.'.png'
    set size ratio aspect_ratio 
  	set title 'chemical field'
      	plot './data/chem/'.i.'.chem.dat' using 1:2:3  with image

 #	set output './visual/buo/'.i.'.png'
 #	set title 'buoyancy field (B_therm*T_z - B_comp * S_z)'
 #   	plot './data/buo/'.i.'.buo.dat' using 1:2:3  with image

  	set output './visual/div/'.i.'.png'
  	set title 'real part of div (u)'
      	plot './data/div/'.i.'.div.dat' using 1:2:3 with image


    # Multiplot for better visibility:-----------------------------------
    
    set terminal pngcairo size 1280,1024 enhanced font 'Verdana,10'
    set output './visual/combo/'.i.'.png'

		set multiplot layout 2,3

    #set cbrange [-0.0:5.0]
		set title 'absolute of u '
     	plot './data/abs_u/'.i.'.abs_u.dat' using 1:2:3  with image
    #unset cbrange

    #set cbrange [-0.8:0.8]
		set title 'z component of u '
      	plot './data/u/'.i.'.u.dat' using 1:2:4  with image
    #unset cbrange

		set title 'vorticity field'
     	plot './data/vort/'.i.'.vort.dat' using 1:2:3  with image

	#	set title 'velocity field'
  #   	plot './data/u/'.i.'.u.dat' using 1:2:3:4  with vectors

		#set title 'real part of div (u)'
    #  	plot './data/div/'.i.'.div.dat' using 1:2:3 with image

   #set cbrange [-3.0:3.0]
		set title 'buoyancy field (B_therm*T_z - B_comp * S_z)'
     	plot './data/buo/'.i.'.buo.dat' using 1:2:3  with image
    #unset cbrange

   #set cbrange [-300.0:300.0]
		set title 'temperature field'
      	plot './data/temp/'.i.'.temp.dat' using 1:2:3  with image
    #load './gnuplot-palettes-master/Greys.pal'    
   #set cbrange [-30.0:30.0]
		set title 'chemical field'
      	plot './data/chem/'.i.'.chem.dat' using 1:2:3  with image
    #load './gnuplot-palettes-master/jet.pal'    
   # unset cbrange


		unset multiplot

    #MULTIPLOT temp/chem--------------------------------------------------
    set terminal pngcairo size 1024,1024 enhanced font 'Verdana,10'
    set output './visual/temp_chem_combo/'.i.'.png'

  	set multiplot layout 2,2

  	set title 'temperature field'
    plot './data/temp/'.i.'.temp.dat' using 1:2:3 with image
  	set title 'fourier spectrum of temp field'
   	plot './data/temp_f/'.i.'.temp_f.dat' using 1:2:3  with image
  	set title 'chemical field'
   	plot './data/chem/'.i.'.chem.dat' using 1:2:3 with image
  	set title 'fourier spectrum of chem field'
   	plot './data/chem_f/'.i.'.chem_f.dat' using 1:2:3 with image

  	unset multiplot
    unset size
}

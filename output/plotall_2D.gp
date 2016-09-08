#load './gnuplot-palettes-master/spectral.pal'    
#load './gnuplot-palettes-master/blues.pal'    
load './gnuplot-palettes-master/jet.pal'    

aspect_ratio = 1
Lx = 32 
Ly = 32  
#Lx = 8 
#Ly = 8  

no_of_img = 299 

    ####################STATISTICS###############################
    set terminal pngcairo size 800,600 enhanced font 'Verdana,10'

    #set output './visual/stat/Energies.png'
    #set title 'simulation time [arb] vs. E_kin [arb]'
    #plot './data/E_stat/E_stat.dat' using 2:3 title "E_{kin}"   ,\
    #     './data/E_stat/E_stat.dat' using 2:4 title "E_{pot}"   ,\
    #     './data/E_stat/E_stat.dat' using 2:5 title "E_{tot}"  
    #
    #set output './visual/stat/Velocities.png'
    #set title 'simulation [arb] time vs. velocities [arb]'
    #plot './data/u_stat/u_stat.dat' using 2:3 title "v_{max}"   ,\
    #     './data/u_stat/u_stat.dat' using 2:4 title "v_{rms}"  
    #
    #
    #set output './visual/stat/Temp.png'
    #set title 'simulation time vs.Temperature measures [arb]'
    #plot './data/T_stat/T_stat.dat' using 2:3 title "max temp"  , \
    #     './data/T_stat/T_stat.dat' using 2:4 title "mean temp"   
    #
    #set output './visual/stat/Chem.png'
    #set title 'simulation time vs. Chemical field measures [arb]'
    #plot './data/C_stat/C_stat.dat' using 2:3 title "max Chem"  , \
    #     './data/C_stat/C_stat.dat' using 2:4 title "mean Chem"   
    #set output './visual/stat/Chem.png'
    #set title 'simulation time vs. Chemical field measures [arb]'
    #plot './data/C_stat/C_stat.dat' using 2:3 title "max Chem"  , \
    #     './data/C_stat/C_stat.dat' using 2:4 title "mean Chem"   
    #########Multiplot###############
    set output './visual/stat/stat_combo.png'
    set multiplot layout 4,2
   #set xrange [0.0:4.0]

        set title 'simulation time vs.Temperature measures [arb]'
        plot './data/T_stat/T_stat.dat' using 2:3 title "max temp"  , \
             './data/T_stat/T_stat.dat' using 2:4 title "mean temp" , \
             './data/T_stat/T_stat.dat' using 2:5 title "min temp"   
        set title 'simulation time vs. Chemical field measures [arb]'
        plot './data/C_stat/C_stat.dat' using 2:3 title "max Chem"  , \
             './data/C_stat/C_stat.dat' using 2:4 title "mean Chem" , \
             './data/C_stat/C_stat.dat' using 2:5 title "min Chem"   


        set title 'simulation [arb] time vs. velocities (all directions) [arb]'
        plot './data/u_stat/u_stat.dat' using 2:3 title "v_{max}"   ,\
             './data/u_stat/u_stat.dat' using 2:4 title "v_{rms}"  
        set title 'simulation [arb] time vs. velocities (specific directions) [arb]'
        plot './data/u_stat/u_stat.dat' using 2:5 title "vx_{max}"   ,\
             './data/u_stat/u_stat.dat' using 2:6 title "vx_{min}"   ,\
             './data/u_stat/u_stat.dat' using 2:7 title "vy_{max}"   ,\
             './data/u_stat/u_stat.dat' using 2:8 title "vy_{min}"  

        set title 'simulation time [arb] vs. E_kin [arb]'
        plot './data/E_stat/E_stat.dat' using 2:3 title "E_{kin}"   ,\
             './data/E_stat/E_stat.dat' using 2:4 title "E_{pot}"   ,\
             './data/E_stat/E_stat.dat' using 2:5 title "E_{tot}"  
        set title 'simulation time vs. maximum of divergence field (should be near zero)'
        plot './data/sys_stat/sys_stat.dat' using 2:3 title "max div"             ,\
             './data/sys_stat/sys_stat.dat' using 2:6 title "max div (brucker)"

        set title 'simulation time vs. shearstrength '
        plot './data/sys_stat/sys_stat.dat' using 2:4 title "shear strength [arb]"
        set title 'simulation time vs. stepwidth dt [arb] '
        plot './data/sys_stat/sys_stat.dat' using 2:5 title "dt [arb]"

        #set title 'simulation time vs. average vorticity [arb] '
        #plot './data/sys_stat/sys_stat.dat' using 2:7 title "average vort [arb]"
    unset multiplot
    ##########################################################################


do for [i=0:no_of_img] {
    set xrange [0:Lx]
    set yrange [0:Ly]
    set size ratio aspect_ratio 

   set terminal pngcairo size 1200,800 enhanced font 'Verdana,10'
#
#  	set output './visual/temp/'.i.'.png'
#  	set title 'temperature field'
#      	plot './data/temp/'.i.'.temp.dat' using 1:2:3  with image notitle
#
#  	set output './visual/u/'.i.'.png'
#  	set title 'velocity field'
#      	plot './data/u/'.i.'.u.dat' using 1:2:($3*50):($4*50)  with vectors
#
#  	set output './visual/abs_u/'.i.'.png'
#  	set title 'absolute magnitude  ofvelocity field'
#      	plot './data/abs_u/'.i.'.abs_u.dat' using 1:2:3 with image notitle
#
  unset xrange
  unset yrange
 	set output './visual/u_f/'.i.'.png'
  set title 'absolute of fourier spectrum'
 	plot './data/u_f/'.i.'.u_f.dat' using 1:2:3  with image notitle
  set xrange [0:Lx]
  set yrange [0:Ly]
#
#    set output './visual/chem/'.i.'.png'
#  	set title 'chemical field'
#      	plot './data/chem/'.i.'.chem.dat' using 1:2:3  with image notitle

	set output './visual/buo/'.i.'.png'
	set title 'buoyancy field (B_therm*T_z - B_comp * S_z)'
 	plot './data/buo/'.i.'.buo.dat' using 1:2:3  with image notitle

    #MULTIPLOT FOR DIVERGENCE

    set terminal pngcairo size 600,400 enhanced font 'Verdana,10'
  	set output './visual/div/'.i.'.png'
		set multiplot layout 1,2
  	set title 'real part of div (u) (real space)'
      	plot './data/div/'.i.'.div.dat' using 1:2:3 with image notitle
  	set title 'real part brucker divergence kbar times u'
      	plot './data/div/'.i.'.div.dat' using 1:2:4 with image notitle
    unset multiplot


    # Multiplot for better visibility:-----------------------------------
    
    set size ratio aspect_ratio 
    set terminal pngcairo size 1200,800 enhanced font 'Verdana,10'
    set output './visual/combo/'.i.'.png'

		set multiplot layout 2,4
    # Axes
    set style line 11 lc rgb '#808080' lt 1
    set border 3 back ls 11
    set tics nomirror out scale 0.75

    set xrange [0:Lx]
    set yrange [0:Ly]

   #set cbrange [-300.0:300.0]
		set title 'temperature field'
      	plot './data/temp/'.i.'.temp.dat' using 1:2:3  with image notitle
    #load './gnuplot-palettes-master/Greys.pal'    
   #set cbrange [-30.0:30.0]
		set title 'chemical field'
      	plot './data/chem/'.i.'.chem.dat' using 1:2:3  with image notitle
    #load './gnuplot-palettes-master/jet.pal'    
   # unset cbrange

    #set cbrange [-0.0:5.0]
		set title 'absolute of u '
     	plot './data/abs_u/'.i.'.abs_u.dat' using 1:2:3  with image notitle
    #unset cbrange

    #set cbrange [-0.8:0.8]
		set title 'z component of u '
      	plot './data/u/'.i.'.u.dat' using 1:2:4  with image notitle
    #unset cbrange

    #set cbrange [-0.05:0.05]
		set title 'buoyancy field (B_therm*T_z - B_comp * S_z)'
     	plot './data/buo/'.i.'.buo.dat' using 1:2:3  with image notitle
    #unset cbrange

		set title 'velocity field'
     	plot './data/u/'.i.'.u.dat' using 1:2:3:4  with vectors

		set title 'vorticity field'
     	plot './data/vort/'.i.'.vort.dat' using 1:2:3  with image notitle


#		set title 'real part of div (u)'
#      	plot './data/div/'.i.'.div.dat' using 1:2:3 with image notitle

    unset xrange 
    unset yrange 

  	set title 'fourier spectrum of chem field'
   	plot './data/chem_f/'.i.'.chem_f.dat' using 1:2:3 with image notitle
    #unset parametric
    #set mapping cartesian
    #set view 60,30,1,1
    #set auto
    #set isosamples 60
    #set hidden3d
   	#splot './data/chem_f/'.i.'.chem_f.dat' using 1:2:3 with image notitle

		unset multiplot

    #MULTIPLOT temp/chem fourier ceck--------------------------------------------------
    set terminal pngcairo size 1200,1200 enhanced font 'Verdana,10'
    set output './visual/spec_inspect_combo/'.i.'.png'
  	set multiplot layout 3,3

    unset tics
    unset colorbox
    #unset border 
  	set title 'temperature field'
    plot './data/temp/'.i.'.temp.dat' using 1:2:3 with image notitle
  	#set title 'Real part of fourier spectrum of temp field'
   	#plot './data/temp_f/'.i.'.temp_f.dat' using 1:2:3  with image notitle
  	set title 'Imag part of fourier spectrum of temp field'
   	plot './data/temp_f/'.i.'.temp_f.dat' using 1:2:4  with image notitle
  	set title 'real cutof part of fourier spectrum of temp field'
   	plot './data/temp_f_remap/'.i.'.temp_f_remap.dat' using 1:2:4  with image notitle

  	set title 'chemical field'
   	plot './data/chem/'.i.'.chem.dat' using 1:2:3 with image notitle
  	set title 'Real part of fourier spectrum of chem field'
   	plot './data/chem_f/'.i.'.chem_f.dat' using 1:2:3 with image notitle
  	#set title 'Imag part of fourier spectrum of chem field'
   	#plot './data/chem_f/'.i.'.chem_f.dat' using 1:2:4 with image notitle
  	set title 'real cutof part of fourier spectrum of chem  field'
   	plot './data/chem_f_remap/'.i.'.chem_f_remap.dat' using 1:2:4  with image notitle

  	set title 'absolute magnitude  of velocity field'
   	plot './data/abs_u/'.i.'.abs_u.dat' using 1:2:3 with image notitle
  	set title 'Real part of fourier transform of u_x '
   	plot './data/u_f/'.i.'.u_f.dat' using 1:2:3 with image notitle
  	#set title 'Imag part of fourier transform of u_x '
   	#plot './data/u_f/'.i.'.u_f.dat' using 1:2:4 with image notitle
  	set title 'Real part of fourier transform of u_y '
   	plot './data/u_f/'.i.'.u_f.dat' using 1:2:5 with image notitle
  	#set title 'Imag part of fourier transform of u_y '
   	#plot './data/u_f/'.i.'.u_f.dat' using 1:2:6 with image notitle


  	unset multiplot
    unset size
}

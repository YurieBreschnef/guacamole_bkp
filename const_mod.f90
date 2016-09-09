module const
	! constants-module: contains simulation wide constants of parameter type
  !if variables are not flagged parameter, they may be temporarily changed (e.g. stepwidth)
	use ISO_C_BINDING 
  implicit none
	include "fftw3.f03"
	integer,parameter				      :: rp 	= 8					                !real-precision
	integer,parameter				      :: real_outp_precision 	= 4					!output precision
	integer,parameter				      :: ip 	= 4					!integer-precision

	integer,parameter				      :: fftw_plan_thoroughness = FFTW_MEASURE
	! possible also FFTW_MEASURE

	integer(kind=ip),parameter		:: xdim	        = 64 
	integer(kind=ip),parameter		:: ydim	        = 64  

	integer(kind = ip),parameter	:: seed 		    = 111111	! seed for random init
	integer(kind = ip),parameter	:: maxfiles 	  = 100 ! maximum no of output files per type
	integer(kind = ip),parameter	:: measure_every= 5  		  ! measure diagnostics every X steps
	integer(kind = ip),parameter	:: debuglevel 	= 1	  		
  ! level 0: no output, level 1: short, level 2: extensive
	real(kind = rp)   ,parameter 	:: pi 		    	= 3.1415926535897932384626433833_rp

	real(kind = rp) ,parameter 		:: Lx	          = 2.0_rp *pi !50.0_rp
	real(kind = rp) ,parameter 		:: Ly	          = 2.0_rp *pi !50.0_rp

	complex(kind = rp),parameter	:: imag		     	= (0.0_rp,1.0_rp)


	integer(kind = ip)	    			:: steps 		
	integer(kind = ip)	    			:: i,j,k,l,main_stp      !used for all kinds of loops

	real(kind = rp),parameter     :: tmax                      = 50.0_rp
	real(kind = rp)					      :: dt 	                     = 1.0e-4_rp

	real(kind = rp)					      :: dt_max                    = 1.0e-3_rp
	real(kind = rp)					      :: dt_min                    = 1.0e-6_rp
	real(kind = rp)					      :: max_step_error            = 1.0e-6_rp
	real(kind = rp)					      :: stepwidth_adjustment_rate = 0.05_rp
	real(kind = rp)					      :: write_intervall           = tmax/real(maxfiles,rp)
	real(kind = rp)					      :: last_written              = 0.0_rp

	real(kind = rp)					      :: dt_2           !(1/2) *dt
	real(kind = rp)					      :: dt_3           !(1/3) *dt
	real(kind = rp)					      :: dt_4           !(1/3) *dt
	real(kind = rp)					      :: dt_8           !(1/3) *dt
	real(kind = rp)					      :: dt_34          !(3/4) * dt
	real(kind = rp)					      :: dt_29          !(2/9) * dt
	real(kind = rp)					      :: dt_49          !(4/9) * dt
	real(kind = rp)					      :: dt_724         !(7/24) * dt

	real(kind = rp)               :: shear    = 0.10_rp
  integer(kind = ip)            :: shearing = 1
	real(kind = rp)               :: sheartime= 0.0_rp
	real(kind = rp)               :: T_rm 
  
  integer(kind = ip)            :: benchmarking = 0

  integer(kind = ip)             :: remapping = 0
  integer(kind = ip)             :: remapping_rate = 0 

	real(kind = rp),parameter     :: D_visc   = 0.7_rp 
	real(kind = rp),parameter			:: D_therm  = 0.010_rp
	real(kind = rp),parameter			:: D_comp   = 0.0001_rp
	real(kind = rp),parameter			:: B_therm  = 4.0_rp
	real(kind = rp),parameter			:: B_comp   = 2.0_rp
	real(kind = rp),parameter			:: S_therm  = 4.0_rp  
	real(kind = rp),parameter			:: S_comp   = 4.0_rp 

	!real(kind = rp),parameter     :: D_visc   = 0.0000000000014_rp 
	!real(kind = rp),parameter			:: D_therm  = 0.00000000000020_rp
	!real(kind = rp),parameter			:: D_comp   = 0.000000000000002_rp
	!real(kind = rp),parameter			:: B_therm  = 0.0000000001_rp
	!real(kind = rp),parameter			:: B_comp   = 0.0000000005_rp
	!real(kind = rp),parameter			:: S_therm  =0.0000000000005_rp  
	!real(kind = rp),parameter			:: S_comp   =0.0000000000001_rp 
end module

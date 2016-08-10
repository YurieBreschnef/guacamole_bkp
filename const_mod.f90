module const
	!constants-module: contains simulation wide constants of parameter type
  !if variables are not flagged parameter, they may be temporarily changed (e.g. stepwidth)
	use ISO_C_BINDING 
  implicit none
	include "fftw3.f03"
	integer,parameter				      :: rp 	= 8					                !real-precision
	integer,parameter				      :: real_outp_precision 	= 4					!output precision
	integer,parameter				      :: ip 	= 4					!integer-precision

	integer,parameter				      :: fftw_plan_thoroughness = FFTW_MEASURE
	! possible also FFTW_MEASURE

	integer(kind=ip),parameter		:: xdim	        = 256 
	integer(kind=ip),parameter		:: ydim	        = 512  

	integer(kind = ip),parameter	:: seed 		    = 123123	! seed for random init
	integer(kind = ip),parameter	:: maxfiles 	  = 300  		! maximum no of output files per type
	integer(kind = ip),parameter	:: measure_every= 5  		! measure diagnostics every X steps
	integer(kind = ip),parameter	:: debuglevel 	= 1	  		
  ! level 0: no output, level 1: short, level 2: extensive
	real(kind = rp)   ,parameter 	:: pi 		    	= 3.1415926535897932384626433833_rp

	real(kind = rp) ,parameter 		:: Lx	          = 100.0 !50.0_rp
	real(kind = rp) ,parameter 		:: Ly	          = 100.0 !50.0_rp

	complex(kind = rp),parameter	:: imag		     	= (0.0_rp,1.0_rp)


	integer(kind = ip)	    			:: steps 		
	integer(kind = ip)	    			:: i,j,k,l,main_stp      !used for all kinds of loops

	real(kind = rp)					      :: dt 	                     = 5.0e-3_rp
	real(kind = rp)					      :: dt_max                    = 5.0e-3_rp
	real(kind = rp)					      :: dt_min                    = 1.0e-6_rp
	real(kind = rp)					      :: max_step_error            = 1.0e-6_rp
	real(kind = rp)					      :: stepwidth_adjustment_rate = 0.05_rp

	real(kind = rp)					      :: tmax     = 600.0_rp
	real(kind = rp)					      :: dt_2           ! (1/2) *dt
	real(kind = rp)					      :: dt_3           ! (1/3) *dt
	real(kind = rp)					      :: dt_4           ! (1/3) *dt
	real(kind = rp)					      :: dt_8           ! (1/3) *dt
	real(kind = rp)					      :: dt_34          ! (3/4) * dt
	real(kind = rp)					      :: dt_29          ! (2/9) * dt
	real(kind = rp)					      :: dt_49          ! (4/9) * dt
	real(kind = rp)					      :: dt_724         ! (7/24) * dt

	real(kind = rp)               :: shear    = 0.00_rp
  integer(kind = ip)            :: shearing = 0

	real(kind = rp),parameter     :: D_visc   = 0.014_rp 

	real(kind = rp),parameter			:: D_therm  = 0.002_rp
	real(kind = rp),parameter			:: D_comp   = 0.00002_rp


	real(kind = rp),parameter			:: B_therm  = 5.0_rp
	real(kind = rp),parameter			:: B_comp   = 5.0_rp

	real(kind = rp),parameter			:: S_therm  = 5.00_rp
	real(kind = rp),parameter			:: S_comp   = 5.0_rp

end module

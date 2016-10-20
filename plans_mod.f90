module plans
	!plans- module: holds fftw-plans for program wide use
  ! integer types are plans for FFTW library
	use const
  implicit none

	integer( kind= 8)										      :: x_xf				! convert x pencil to xf
	integer( kind= 8)										      :: xf_x 			! ..and back	
	integer( kind= 8)										      :: y_yf				! analog
	integer( kind= 8)										      :: yf_y

	integer( kind= 8)										      :: y_yf2				! analog
	integer( kind= 8)										      :: yf_y2

	integer( kind= 8)										      :: full2D
	integer( kind= 8)										      :: ifull2D

  complex(kind = rp),dimension(0:xdim-1)		:: x_pen,x_pen_f	! 1D pencil in x-dir
  complex(kind = rp),dimension(0:ydim-1)		:: y_pen,y_pen_f	! 1D pencil in y-dir
  complex(kind = rp),dimension(0:(2*ydim)-1):: y_pen2,y_pen_f2	! 1D pencil in y-dir (zero padded)

end module

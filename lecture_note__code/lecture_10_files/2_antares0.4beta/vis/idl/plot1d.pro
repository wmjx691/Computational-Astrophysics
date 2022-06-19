;+
; Plot animate 
;
; Kuo-Chuan Pan Oct,19,2008
;
; argumetns:
;   istart   start number of frame
;   istep    step
;   iend     end number of frame
;   var
;-

pro plot1d, istart, istep, iend, var

 window, xsize=600, ysize=600
 name = findfile('den*.tab',count=N)

 For i=istart, iend, istep  do $
   begin
   file=name(i)
   readcol, file, R, rho, u, $
	    FORMAT='D,D,D'

   case var>0<4 of
  
   1:  begin   
       plot, R, rho, $
     	   ;xran=[0,256],yran=[0,1.1], $
	   xtit='R', ytit='Density', $
	   tit =' Sod test' ,psym=4
       end   
   2:  begin
       plot, R, u, $
     	   xran=[-5,5],yran=[-1.,0.1], $
	   xtit='R', ytit='velocity', $
	   tit =' CPA HW' ,psym=4
       end   
   3:  begin
       plot, R, P, $
     	   xran=[-5,5],yran=[-.1,1.1], $
	   xtit='R', ytit='Pressure', $
	   tit =' CPA HW' ,psym=4
       end   
   4:  begin 
       plot, R, e, $
     	   xran=[-5,5],yran=[-.1,2.1], $
	   xtit='R', ytit='internal energy', $
	   tit =' CPA HW' ,psym=4
       end   
   endcase	   
  endfor

end
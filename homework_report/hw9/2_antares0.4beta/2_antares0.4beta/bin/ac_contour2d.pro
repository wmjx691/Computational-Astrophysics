;+
;
; Plot 2D contour
;
;
;-

pro ac_contour2d,istart,istep,iend,nend,var

  name = findfile('dat*.bin',count=N)
  window,1,xsize=600,ysize=600

  For i = istart, iend, istep do $
	  begin
	  file=name(i)

	  data=dblarr(nend+4.,nend+4.,1,8)  ; for mhd
;	  data=dblarr(nend+4.,nend+4.,1,5)  ; for hydro
 
	  openr, lun, file, /GET_LUN, /F77_UNFORMATTED
	  readu, lun, data
	  FREE_LUN, lun

	  contour,data(*,*,0,0), nlevels=50   ; density


  endfor
end



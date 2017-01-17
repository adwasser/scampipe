;+
;given x and y coordinates and a chip and a date, this program will apply the distortion correction and return the new coordinates
;-


function conv, r, a

R1 = 0d
for i=0d,(n_elements(a)-1d) do R1 = R1 + a[i]*R^i
return, R1

end




function reverse_mapping, R, a=a, r1_actual=r1_actual

return, total(abs(conv(r,a) - R1_actual))

end






pro supcam_distcor_init, date, chip, xcen, ycen, xcen1, ycen1, t
if (date lt 2008.58) then begin
	CASE chip OF
	'w67c1'	:	BEGIN
					xcen = 5310.051604
					ycen = -23.288532
					xcen1 = 5249.05 + 0.5
					ycen1 = -26.7885 - 0.25
					t=0.001395
					yrange = [-5,5]
				END
	'w6c1'	:	BEGIN
					xcen = 3179.092082
					ycen = -30.946680
					xcen1 = 3167.09
					ycen1 = -29.4467 + 0.25
					t=-0.000157
				END
	'si005s':	BEGIN
					xcen = 1052.085930
					ycen = 4084.168020
					xcen1 = 1053.09
					ycen1 = 4056.17
					t=-0.000000
				END
	'si001s':	BEGIN
					xcen = -1073.302726
					ycen = 4086.766838
					xcen1 = -1063.30
					ycen1 = 4086.77 - 30
					t=-0.000000
				END
	'si002s'	:	BEGIN
					xcen = -1068.140361
					ycen = -37.513745
					xcen1 = -1066.14 - 0.6
					ycen1 = -29.5137 + 0.25
					t=0.002169
				END
	'si006s'	:	BEGIN
					xcen = 1052.567093
					ycen = -29.154287
					xcen1 = 1053.57
					ycen1 = -28.1543
					t=0.000382
				END
	'w93c2'	:	BEGIN
					xcen = 5305.048187
					ycen = 4084.741529
					xcen1 = 5244.05
					ycen1 = 4042.74
					t=-0.000998
				END
	'w9c2':	BEGIN
					xcen = 3176.996541
					ycen = 4087.737305
					xcen1 = 3164.00 + 0.6
					ycen1 = 4057.74 - 0.25
					t=-0.000294
				END
	'w4c5':	BEGIN
					xcen = -3201.998748
					ycen = 4077.553770
					xcen1 = -3158.00
					ycen1 = 4038.05
					t=0.001649
				END
	'w7c3'	:	BEGIN
					xcen = -3195.190437
					ycen = -36.348975
					xcen1 = -3161.
					ycen1 = -31.349
					t=0.000778
				END
	ENDCASE
ENDIF
end




;Given the coordinates from a distortion corrected image, this program returns those coordinates in the non-corrected image
;This isn't perfect, there is about a 1-2 pixel rms at worst
pro supcam_distcorr, x, y, chip, date, x1, y1

chips = ['w67c1', 'w6c1', 'si006s', 'si002s', 'w7c3', 'w93c2', 'w9c2', 'si005s', 'si001s', 'w4c5']
a = [0d,1d,7.16417d-08,3.03146d-10,5.69338d-14,-6.61572d-18,0d,0d]						;polynomial coefficients of the distortion

supcam_distcor_init, date, chip, xcen, ycen, xcen1, ycen1, t								;get the xcen,ycen,xcen1,ycen1,t parameters

x0 = x - xcen1																			;convert to the mosaic coordinate system
y0 = y - ycen1
r1_actual = sqrt(x0^2d + y0^2d)															;calclate the radius of each object

r1 = fltarr(n_elements(r1_actual))
quiet = 1
for j=0,(n_elements(r1_actual)-1l) do begin												;find the original radius of each object in the non-distorted image
	start = r1_actual[j]
	parms = tnmin('reverse_mapping', start, functargs={a:a, r1_actual:r1_actual[j]}, bestmin=bestmin, /autoderivative, quiet=quiet, maxiter=maxiter)
	R1[j] = parms	;conv(parms,a)
	if (bestmin gt 1d-5) then print, 'bestmin is bad: ' + strtrim(bestmin,2)
endfor

xyrot, x0*r1_actual/r1, y0*r1_actual/r1, t, x1, y1, trans=[xcen,ycen]					;derotate the coordinates and return to the chip coordinate system

;openw, lun, '1.reg', /get_lun
;printf, lun, [transpose(x1),transpose(y1)]
;free_lun, lun
;print, chip




end


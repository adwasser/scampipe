
pro jaa_astrom, x, y, ra, dec, file, error, nsig=nsig, maxiter=maxiter, noupdate=noupdate, result=result, quick=quick, weights=weights, quiet=quiet, start=start
	COMPILE_OPT idl2, HIDDEN 
	
	
	hdr = headfits(file)											;retrieve header from file
	hdr0 = hdr
	
	extast, hdr, astr
	
	
	if (n_elements(nsig) eq 0) then nsig = 3d
	if (n_elements(maxiter) eq 0) then maxiter = 0
	
	
	
	crval1  = [sxpar(hdr, 'CRVAL1'), sxpar(hdr, 'CRVAL2')]			;reference RA/DEC
	crpix1  = [sxpar(hdr, 'CRPIX1'), sxpar(hdr, 'CRPIX2')] + 1.		;corresponding reference x/y (have to subtract off 1 because fits indices start at 1)
	
	
	
	crval1 = [mean(minmax(float(ra))), mean(minmax(float(dec)))]
	crpix1 = [mean(minmax(float(x))), mean(minmax(float(y)))]
	astr.crval = crval1
	astr.crpix = crpix1
	
	
	wcssph2xy, [crval1[0], ra], [crval1[1], dec], xi, eta, 2, latpole = 90.0, crval=crval1, crxy=crpix1/(3600d/0.202d)		;tangent deprojection of ra, dec
	;if (n_elements(error) eq 0) then stop
	;error = error*3600d/0.202d				;convert from degrees to pixels
	
	;ind = where(error eq 0, count)
	;if (count gt 0) then error[ind] = 5.
	
	mxi = xi[0]														;first xi value corresponds to the reference ra value
	xi = xi[1:*]													;remove the reference value from the list
	
	meta = eta[0]													;first eta value corresponds to the reference dec value
	eta = eta[1:*]													;remove the reference value from the list
	
	
	
	
	use = indgen(n_elements(xi))									;indices of good objects
	
	nzones = [3,5]
	
	minxi = min(xi)
	mineta = min(eta)
	
	xi_zone = intarr(nzones[0])
	eta_zone = intarr(nzones[1])
	xi_inc = difference(minmax(xi))/float(nzones[0])
	eta_inc = difference(minmax(eta))/float(nzones[1])
	
	dstart = [(linfit(x[use] - crpix1[0], -(xi[use] - mxi)*3600d))[0],(linfit(y[use] - crpix1[1], (eta[use] - meta)*3600d, yfit=yfit))[0],0.202d,0.202d,0d]
	if (n_elements(start) gt 0) then begin
	;	if ((min(finite(start)) eq 0) or (abs(start[0]) gt 10) or (abs(start[1]) gt 10) or (abs(start[2] - 0.202) gt 0.05) or (abs(start[3] - 0.202) gt 0.05)) then start = dstart
		if ((min(finite(start)) eq 0) or (abs(start[2] - 0.202) gt 0.05) or (abs(start[3] - 0.202) gt 0.05)) then start = dstart
	endif else start = dstart
	;if ((abs(start[0]) gt 10) or (abs(start[1]) gt 10) or (abs(start[2] - 0.202) gt 0.05) or (abs(start[3] - 0.202) gt 0.05)) then start = [0d,0d,0.202d,0.202d,0d]
	
	
	
	for i=-1,(maxiter - 1) do begin									;for each iteration (atleast 1 will happen even if maxiter = 0)
	
		;now to relate xi,eta to x, y
		res = xyoff_rot((x[use] - crpix1[0]), y[use] - crpix1[1], -(xi[use] - mxi)*3600d, (eta[use] - meta)*3600d, error, start=start, fixed=[0,0,0,0,0],yfit=yfit, weights=weights, quiet=quiet)		;rotation and magnificaion coefficients
	
	
		xidev = -(xi[use] - mxi)*3600d - yfit[*,0]					;deviation of each xi value from the fit calculated above
		etadev = (eta[use] - meta)*3600d - yfit[*,1]				;deviation of each eta value from the fit calculated above
		dev = sqrt(xidev^2d + etadev^2d)							;total deviation
	
	
	
		;astrolib
		;!priv = 2
		;dbcreate, 'db_astrometry', 1, 1
		;dbopen, 'db_astrometry',1
	
	
	
		xirms = sqrt(total(xidev^2d) /float(n_elements(use)))			;rms of xi
		etarms = sqrt(total(etadev^2d) /float(n_elements(use)))			;rms of eta
	
	;	print, xirms, etarms
	
		dbname = 'db_astrometry.sav'
		dummy = file_search(dbname, count=count)
		dbentry = create_struct((strsplit(file,'.fits',/extract,/regex))[0], {xirms:xirms, etarms:etarms, xi:xi, xidev:xidev, eta:eta, etadev:etadev})
		if (count ne 0) then begin
			restore, dbname
			tagnames = tag_names(db)
			ind = where(tagnames+'.FITS' eq strupcase(file), nmatch)
			if (nmatch ne 0) then db = struct_trimtags(db, except=strupcase((strsplit(file,'.fits',/extract,/regex))[0]))
			if (type(db,/string) ne 'structure') then db = dbentry else db = struct_addtags(db, dbentry)
		endif else db = dbentry
		save, db, dbname, filename=dbname
	
	
		cleanplot, /silent
		erase
		!p.multi=[0,2,2]
		multiplot
		plot, xi[use], xidev, psym=1, ytitle='RMS Xi [arcsec]', YRANGE=[-3,3], ystyle=1
		oplot, (indgen(10)-5)*100, intarr(10), linestyle=2, color=fsc_color('red')
		multiplot
		plot, eta[use], xidev, psym=1, YRANGE=[-3,3], ystyle=1
		oplot, (indgen(10)-5)*100, intarr(10), linestyle=2, color=fsc_color('red')
		multiplot
		plot, xi[use], etadev, psym=1, xtitle='Xi [arcsec]', ytitle='RMS Eta [arcsec]', YRANGE=[-3,3], ystyle=1
		oplot, (indgen(10)-5)*100, intarr(10), linestyle=2, color=fsc_color('red')
		multiplot
		plot, eta[use], etadev, psym=1, xtitle='Eta [arcsec]', YRANGE=[-3,3], ystyle=1
		oplot, (indgen(10)-5)*100, intarr(10), linestyle=2, color=fsc_color('red')
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		if (i eq (maxiter - 1)) then break								;quit on last iteration
	
	
	
		nuse = n_elements(use)											;number of remaining objects from last iteration
	;	if (i eq 0) then ind = where( dev le max([3.,nsig])*max([0.05,sqrt(robust_sigma(xidev)^2d + robust_sigma(etadev)^2d)]) , count)	$;sigma clipped list, at least 3 sigma for first iteration
		if (i eq 0) then ind = where( dev le 0.6 , count)	$		;first big cut of outliers
					else ind = where( dev le nsig*max([0.05,sqrt(robust_sigma(xidev)^2d + robust_sigma(etadev)^2d)]) , count)	;sigma clipped list
		for z=0,(nzones[0] - 1) do begin
			dummy = where( (xi[use] ge (minxi + z*xi_inc)) and (xi[use] lt (minxi + (z+1)*xi_inc)), zcount)
			xi_zone[z] = zcount
		endfor
		for z=0,(nzones[1] - 1) do begin
			dummy = where( (eta[use] ge (mineta + z*eta_inc)) and (eta[use] lt (mineta + (z+1)*eta_inc)), zcount)
			eta_zone[z] = zcount
		endfor	
	
	
	;	use = use[where( (xidev le nsig*robust_sigma(xidev[use])) and (etadev le nsig*robust_sigma(etadev[use]))  , count)]	;sigma clipped list
			
		if (count eq nuse) then break									;quit if no objects were clipped
		if (count lt 10) then break										;quit if the object list is too small
		if (i ge 0) then begin
			if (min(xi_zone) le round(20/float(nzones[0]))) then break		;quit if the distribution of points is poorly distributed in xi							
			if (min(eta_zone) le round(20/float(nzones[1]))) then break		;quit if the distribution of points is poorly distributed in eta
			if (count lt 20) then break										;quit if the object list gets small, allow one iteration first
		endif
		start = res
	
	
		cleanplot, /silent
		erase
		!p.multi=[0,2,2]
		multiplot
		plot, xi[use], xidev, psym=1, ytitle='RMS Xi [arcsec]', YRANGE=[-2,2], ystyle=1
		multiplot
		plot, eta[use], xidev, psym=1, YRANGE=[-2,2], ystyle=1
		multiplot
		plot, xi[use], etadev, psym=1, xtitle='Xi [arcsec]', ytitle='RMS Eta [arcsec]', YRANGE=[-2,2], ystyle=1
		multiplot
		plot, eta[use], etadev, psym=1, xtitle='Eta [arcsec]', YRANGE=[-2,2], ystyle=1
	
		use = use[ind]
	
	
	endfor
	
	
	
	
	;cd = fltarr(2,2)
	;cd[0,0] = (cd1_1 = cos(res[4])/3600d*0.202d)									;cd matrix
	;cd[1,0] = (cd1_2 = -sin(res[4])/3600d*0.202d)
	;cd[0,1] = (cd2_1 = sin(res[4])/3600d*0.202d)
	;cd[1,1] = (cd2_2 = cos(res[4])/3600d*0.202d)
	
	
	;astrom1, ra, dec, indgen(n_elements(ra)), x, y, astr, err, 3600d/0.202d*!dtor, /display
	
	;print, xirms, etarms, res[2], res[3], nuse
	
	
	if (n_elements(noupdate) eq 0) then begin
	
		;astr.cdelt = [res[2],res[3]]/3600d
		;astr.crpix = astr.crpix + [res[0]/res[2], -res[1]/res[3]]
		;astr = struct_trimtags(astr, except='CD')
		;astr = struct_trimtags(astr, except='PV2')
		;astr.cd = cd
		;putast, hdr, astr, cd_type=0
		
		sxaddpar, hdr, 'CRVAL1', astr.crval[0], /saveComment
		sxaddpar, hdr, 'CRVAL2', astr.crval[1], /saveComment
		sxaddpar, hdr, 'CDELT1', -res[2]/3600d, /saveComment
		sxaddpar, hdr, 'CDELT2', res[3]/3600d, /saveComment
		sxaddpar, hdr, 'CRPIX1', astr.crpix[0] - res[0]/res[2] + 1, /saveComment
		sxaddpar, hdr, 'CRPIX2', astr.crpix[1] - res[1]/res[3] + 1, /saveComment
		sxaddpar, hdr, 'CROTA2', -res[4]*!radeg, /saveComment, BEFORE='HISTORY'
	
		sxdelpar, hdr, 'CROTA'
		sxdelpar, hdr, 'CD1_1'
		sxdelpar, hdr, 'CD1_2'
		sxdelpar, hdr, 'CD2_1'
		sxdelpar, hdr, 'CD2_2'
	
	;	sxaddpar, hdr, 'CD1_1' + alt, cd[0,0], degpix, 'HISTORY',/SaveC
	;	sxaddpar, hdr, 'CD2_1' + alt, cd[1,0], degpix, 'HISTORY',/SaveC
	;	sxaddpar, hdr, 'CD1_2' + alt, cd[0,1], degpix, 'HISTORY',/SaveC
	;	sxaddpar, hdr, 'CD2_2' + alt, cd[1,1], degpix, 'HISTORY',/SaveC
	
	;print, astr.crpix
		
	;	sxaddpar, hdr, 'CRPIX1', sxpar(hdr, 'CRPIX1') + res[0]/res[4], /saveComment
	;	sxaddpar, hdr, 'CRPIX2', sxpar(hdr, 'CRPIX2') - res[1]/res[5], /saveComment
	
	;	sxaddpar, hdr, 'CD1_1', cd1_1, /saveComment							;update cdmatrix in header
	;	sxaddpar, hdr, 'CD1_2', cd1_2, /saveComment
	;	sxaddpar, hdr, 'CD2_1', cd2_1, /saveComment
	;	sxaddpar, hdr, 'CD2_2', cd2_2, /saveComment
			
		
		;stop
		;djs_modfits, file, 0, hdr, /silent									;save new header to original file
		modfits, file, 0, hdr		;should i ditch djs_modfits?
	endif
	
	start = res
	result = {xscale:res[2], yscale:res[3], xirms:xirms, etarms:etarms, xrot:res[4]*!radeg, file:file, use:use, xidev:xidev, etadev:etadev}
	
	
end










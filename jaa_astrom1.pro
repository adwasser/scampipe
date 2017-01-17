
FUNCTION jaa_astrom_func, dum, arg, crpix1, crpix2, crota2, cdelt1, cdelt2, docdelt=docdelt

	CDELT1 = arg[0]
	CDELT2 = arg[1]
	CRVAL1 = double(arg[2])
	CRVAL2 = double(arg[3])
	crpix1 = arg[4]
	crpix2 = arg[5]

	;measure xi and eta given the inputs
	wcssph2xy, dum.ra, dum.dec, xi, eta, 2, latpole = 90.0, crval=[crval1,crval2]
	
	;measure the inferred crpix1, crpix2, crota2
	fixed = [0,0,1,1,0]	
	weights = dum.x*0 + 1d
	dstart = [-1*crpix1, -1*crpix2, cdelt1, cdelt2, 0d]
	quiet = 1
	

	res0 = xyrot_off(dum.x, dum.y, xi, eta, error, start=dstart, fixed=fixed, yfit=yfit, weights=weights, quiet=quiet, /newt)		;rotation and magnificaion coefficients
	resid = robust_sigma([xi,eta] - yfit)*3600d
	;print, resid
	res = xyrot_off(dum.x, dum.y, xi, eta, error, start=res0, fixed=[1,1,0,0,0], yfit=yfit, weights=weights, quiet=quiet, /newt)		;rotation and magnificaion coefficients
	resid = robust_sigma([xi,eta] - yfit)*3600d
	;print, resid
		
	crpix1 = -res[0]
	crpix2 = -res[1]
	cdelt1 = res[2]
	cdelt2 = res[3]
	crota2 = res[4]

	newXi = yfit[*,0]
	newEta = yfit[*,1]

	;convert xi and eta to ra and dec
	WCSXY2SPH, newXi, newEta, newRA, newDec, CRVAL = [crval1,crval2], latpole=90.
	new = struct_addtags(dum, arr_struct({xra:newRA, xdec:newDec}))
	return, new

END

;-------------------------------------------------------------------------------------------

PRO jaa_astrom1, x, y, ra, dec, file, error, nsig=nsig, maxiter=maxiter, noupdate=noupdate, result=result, $
		quick=quick, weights=weights, quiet=quiet, start=start
	COMPILE_OPT idl2, HIDDEN 
	
	hdr = headfits(file)											;retrieve header from file
	hdr0 = hdr
	extast, hdr, astr
	
	if (n_elements(nsig) eq 0) then nsig = 3d
	if (n_elements(maxiter) eq 0) then maxiter = 0
	
	crval1 = sxpar(hdr, 'CRVAL1') 
	crval2 = sxpar(hdr, 'CRVAL2')			;reference RA/DEC
	crpix1 = sxpar(hdr, 'CRPIX1')
	crpix2 = sxpar(hdr, 'CRPIX2')		;corresponding reference x/y (have to subtract off 1 because fits indices start at 1)

	CDELT1 = -5.6093104d-05							;initial guess at RA scale (deg/pixel)
	CDELT2 = 5.6095829d-05
	crota2 = 0d

	data = arr_struct({x:x, y:y, ra:ra, dec:dec})
	
	
	use = indgen(n_elements(data))									;indices of good objects
	nzones = [3,5]	
	minra = min(ra)
	mindec = min(dec)
	ra_zone = intarr(nzones[0])
	dec_zone = intarr(nzones[1])
	ra_inc = difference(minmax(data.ra))/float(nzones[0])
	dec_inc = difference(minmax(data.dec))/float(nzones[1])
	
	
	
	for i=-1,(maxiter - 1) do begin										;for each iteration (atleast 1 will happen even if maxiter = 0)
	
		arg = [cdelt1, cdelt2, crval1, crval2, crpix1, crpix2, crota2]

		res = jaa_astrom_func(data[use], arg, crpix1, crpix2, crota2, cdelt1, cdelt2)					;fit crpix1, crpix2
				
		radev = (res.ra - res.xRA)*3600d
		decdev = (res.dec - res.xDEC)*3600d
		dev = sqrt(raDev^2d + decDev^2d)
		rarms = sqrt(total(radev^2d) /float(n_elements(use)))			;rms of xi
		decrms = sqrt(total(decdev^2d) /float(n_elements(use)))			;rms of eta
	
	
		cleanplot, /silent
		erase
		!p.multi=[0,2,2]
		multiplot
		plot, -1*(data[use].ra - crval1)*60d, radev, psym=1, ytitle='RMS RA [arcmin]', YRANGE=[-3,3], ystyle=1, title=textoidl('\alpha_{RMS} = ')+roundx(RARMS,2)
		oplot, (indgen(10)-5)*1d5, intarr(10), linestyle=2, color=fsc_color('red')
		multiplot
		plot, (data[use].dec - crval2)*60d, decdev, psym=1, YRANGE=[-3,3], ystyle=1, title=textoidl('\delta_{RMS} = ')+roundx(DECRMS,2)
		oplot, (indgen(10)-5)*1d5, intarr(10), linestyle=2, color=fsc_color('red')
		multiplot
		plot, -1*(data[use].ra - crval1)*60d, radev, psym=1, xtitle='RA [arcmin]', ytitle='RMS DEC [arcmin]', YRANGE=[-3,3], ystyle=1
		oplot, (indgen(10)-5)*1d5, intarr(10), linestyle=2, color=fsc_color('red')
		multiplot
		plot, (data[use].dec - crval2)*60d, decdev, psym=1, xtitle='DEC [arcmin]', YRANGE=[-3,3], ystyle=1
		oplot, (indgen(10)-5)*1d5, intarr(10), linestyle=2, color=fsc_color('red')
		
		if (i eq (maxiter - 1)) then break								;quit on last iteration
	
	
		nuse = n_elements(use)											;number of remaining objects from last iteration
		if (i eq 0) then ind = where( dev le 0.6 , count)	$		;first big cut of outliers
					else ind = where( dev le nsig*max([0.05,sqrt(robust_sigma(radev)^2d + robust_sigma(decdev)^2d)]) , count)	;sigma clipped list
		for z=0,(nzones[0] - 1) do begin
			dummy = where( (data[use].ra ge (minra + z*ra_inc)) and (data[use].ra lt (minra + (z+1)*ra_inc)), zcount)
			ra_zone[z] = zcount
		endfor
		for z=0,(nzones[1] - 1) do begin
			dummy = where( (data[use].dec ge (mindec + z*dec_inc)) and (data[use].dec lt (mindec + (z+1)*dec_inc)), zcount)
			dec_zone[z] = zcount
		endfor	
			
		if (count eq nuse) then break										;quit if no objects were clipped
		if (count lt 10) then break											;quit if the object list is too small
		if (i ge 0) then begin
			if (min(ra_zone) le round(20/float(nzones[0]))) then break		;quit if the distribution of points is poorly distributed in xi							
			if (min(dec_zone) le round(20/float(nzones[1]))) then break		;quit if the distribution of points is poorly distributed in eta
			if (count lt 20) then break										;quit if the object list gets small, allow one iteration first
		endif
	
		use = use[ind]
		
	endfor
	
	
	dbname = 'db_astrometry.sav'
	dummy = file_search(dbname, count=count)
	dbentry = create_struct((strsplit(file,'.fits',/extract,/regex))[0], {rarms:rarms, decrms:decrms, ra:data[use].ra, radev:radev, dec:data[use].dec, decdev:decdev})
	if (count ne 0) then begin
		restore, dbname
		tagnames = tag_names(db)
		ind = where(tagnames+'.FITS' eq strupcase(file), nmatch)
		if (nmatch ne 0) then db = struct_trimtags(db, except=strupcase((strsplit(file,'.fits',/extract,/regex))[0]))
		if (type(db,/string) ne 'structure') then db = dbentry else db = struct_addtags(db, dbentry)
	endif else db = dbentry
	save, db, dbname, filename=dbname
		
	
	
	

	
	if (n_elements(noupdate) eq 0) then begin
		
		sxaddpar, hdr, 'CRVAL1', crval1, /saveComment
		sxaddpar, hdr, 'CRVAL2', crval2, /saveComment
		sxaddpar, hdr, 'CDELT1', cdelt1, /saveComment
		sxaddpar, hdr, 'CDELT2', cdelt2, /saveComment
		sxaddpar, hdr, 'CRPIX1', crpix1, /saveComment
		sxaddpar, hdr, 'CRPIX2', crpix2, /saveComment
		sxaddpar, hdr, 'CROTA2', crota2*!radeg, /saveComment, BEFORE='HISTORY'
	
		sxdelpar, hdr, 'CROTA'
		sxdelpar, hdr, 'CD1_1'
		sxdelpar, hdr, 'CD1_2'
		sxdelpar, hdr, 'CD2_1'
		sxdelpar, hdr, 'CD2_2'
		
		modfits, file, 0, hdr		;should i ditch djs_modfits?
	endif
	
	result = {xscale:cdelt1*3600d, yscale:cdelt2*3600d, rarms:rarms, decrms:decrms, xrot:crota2*!radeg, file:file, use:use, radev:radev, decdev:decdev}
	
	if ((RARMS gt 5) or (DECRMS gt 5)) then begin							;throw an error if the residual variance isn't small enough
		print, 'The residual variance on the astrometric solution is too big.'
		print, 'Something may have gone wrong...'
		choice = choose('Continue?', 'yes')
		if ((strupcase(choice) ne 'YES') and (strupcase(choice) ne 'Y')) then stop
	endif

END














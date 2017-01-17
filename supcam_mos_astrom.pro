FUNCTION supcam_mos_astrom_getcrpix, arg1
	COMPILE_OPT idl2, HIDDEN 
	common sharemos1, dum, crval1, crval2, cdelt1, cdelt2
	arg = arg1
	crpix1 = arg[0]
	crpix2 = arg[1]
	wcssph2xy, dum.ra, dum.dec, xi, eta, 2, latpole = 90.0, crval=[crval1,crval2]
	weights = dum.x*0 + 1d
	dstart = [0d, 0d, cdelt1, cdelt2, 0d]
	;quiet = 1
	res = xyoff_rot(dum.x - crpix1, dum.y - crpix2, xi, eta, error, start=dstart, fixed=[0,0,1,1,0], yfit=yfit, weights=weights, quiet=quiet, /newt)		;rotation and magnificaion coefficients
	return, res[0:1]
END
	
;-------------------------------------------------------------------------------------------
	
FUNCTION supcam_mos_astrom_getcrpix1, arg1
	COMPILE_OPT idl2, HIDDEN 
	common sharemos1
	arg = arg1
	crpix1 = arg[0]
	crpix2 = arg[1]
	wcssph2xy, dum.ra - 0*crval1, dum.dec - 0*crval2, xi, eta, 2, latpole = 90.0, crval=[crval1,crval2]
	weights = dum.x*0 + 1d
	dstart = [-1*crpix1, -1*crpix2, cdelt1, cdelt2, 0d]
	quiet = 1
	res = xyrot_off(dum.x, dum.y, xi, eta, error, start=dstart, fixed=[0,0,1,1,0], yfit=yfit, weights=weights, quiet=quiet, /newt)		;rotation and magnificaion coefficients
	return, res
END

;-------------------------------------------------------------------------------------------

FUNCTION supcam_mos_astrom_install, arg0
;--install these into the catalog
	COMPILE_OPT idl2, HIDDEN 
	common sharemos, data0, gind, crpix1, crpix2, crota2, counter, flag, results, crval2_0
	common sharemos1
	
	arg = arg0
	CDELT1 = arg[0]
	CDELT2 = arg[1]
	CRVAL1 = double(arg[2])
	CRVAL2 = double(arg[3])
	
	for i=0,(n_tags(data0)-1) do begin
		dum = data0.(i)
		arg1 = [crpix1[i], crpix2[i]]
		res = supcam_mos_astrom_getcrpix1(arg1)
		;stop
		;nres = newton(arg1, 'supcam_mos_astrom_getcrpix')					;find the crpix1,2 that minimize the xi,eta residuals
		crpix1[i] = -res[0]
		crpix2[i] = -res[1]
		crota2[i] = res[4]	
		xyrot, dum.x - crpix1[i], dum.y - crpix2[i], crota2[i], xi, eta, mult=[cdelt1,cdelt2]
		WCSXY2SPH, xi, eta, ra, dec, CRVAL = [crval1,crval2], latpole=90.
		dum = struct_addtags(dum, arr_struct({xra:ra, xdec:dec, chip:ra*0+i}))
		cat = struct_append(cat, dum)
	endfor
	if (n_elements(flag) gt 0) then begin
		cat = cat[gind]
	endif
	++counter
	return, cat
END

;-------------------------------------------------------------------------------------------

FUNCTION supcam_mos_astrom_newtify, arg
	COMPILE_OPT idl2, HIDDEN 
	common sharemos
	
	cat = supcam_mos_astrom_install(arg)						;install the guesses into the catalog
	
	result = sqrt( (robust_sigma(cat.ra - cat.xra, /zero)*cos(CRVAL2_0*!dtor))^2d + robust_sigma(cat.dec - cat.xdec, /zero)^2d )
	
	if (n_elements(results) eq 0) then results = result else results = [results,result]
	
	if ((counter mod 1) eq 0) then begin
		if 0 then begin
			cleanplot, /silent
			erase
			!p.multi=[0,2,2]
			yrange = [-2,2]
			multiplot
			plot, cat.ra, (cat.ra - cat.xra)*3600d, xrange=reverse(minmax(cat.ra)), yrange=yrange, xstyle=1, ystyle=1, psym=1, ytitle=textoidl('\Delta RA')
			oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
			multiplot
			plot, cat.dec, (cat.ra - cat.xra)*3600d, xrange=minmax(cat.dec), yrange=yrange, xstyle=1, ystyle=1, psym=1
			oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
			multiplot
			plot, cat.ra, (cat.dec - cat.xdec)*3600d, xrange=reverse(minmax(cat.ra)), yrange=yrange, xstyle=1, ystyle=1, psym=1, ytitle=textoidl('\Delta DEC'), xtitle='RA [degrees]'
			oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
			multiplot
			plot, cat.dec, (cat.dec - cat.xdec)*3600d, xrange=minmax(cat.dec), yrange=yrange, xstyle=1, ystyle=1, psym=1, xtitle='DEC [degrees]'
			oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
		endif else supcam_mos_astrom_plot_resid, cat
	
		window, 3
		cleanplot, /silent
		plot, cat.ra*cos(arg[3]*!dtor), cat.dec, psym=1, /iso, xrange=reverse(minmax(cat.ra)*cos(arg[3]*!dtor)), yrange=minmax(cat.dec), xstyle=1, ystyle=1
		oplot, cat.xra*cos(arg[3]*!dtor), cat.xdec, psym=1, color=fsc_color('red')
	
		window, 2
		cleanplot, /silent
		plot, [results]*3600d, psym=-1, /ylog, xtitle='# Sub-Iterations', ytitle='RMS residuals [arcsec]'
		wset, 0
	endif
	print, [arg, result]
	return, result
END

;-------------------------------------------------------------------------------------------

PRO supcam_mos_astrom_plot_resid, cat, ps=ps, outfile=outfile

	if (n_elements(ps) gt 0) then begin
		psopen, outfile, /encapsulated, /color, xs=8, ys=6, /inches
		!p.charsize = 1.2
		!p.charthick = 3
	endif

	plotsym, 0, 0.4, /fill
	cleanplot, /silent
	erase
	!p.multi=[0,2,2]
	yrange = [-1,1]*0.6
	multiplot
	plot, cat.ra, (cat.ra - cat.xra)*3600d, xrange=reverse(minmax(cat.ra)), yrange=yrange, xstyle=1, ystyle=1, psym=8, ytitle=textoidl('\Delta RA'), title=textoidl('\alpha_{RMS} = ')+roundx(robust_sigma(cat.ra - cat.xra)*3600d,3)+'"'
	oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
	sharpcorners, thick=3
	multiplot
	plot, cat.dec, (cat.ra - cat.xra)*3600d, xrange=minmax(cat.dec), yrange=yrange, xstyle=1, ystyle=1, psym=8, title=textoidl('\delta_{RMS} = ')+roundx(robust_sigma(cat.dec - cat.xdec)*3600d,3)+'"'
	oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
	sharpcorners, thick=3
	multiplot
	plot, cat.ra, (cat.dec - cat.xdec)*3600d, xrange=reverse(minmax(cat.ra)), yrange=yrange, xstyle=1, ystyle=1, psym=8, ytitle=textoidl('\Delta DEC'), xtitle='RA [degrees]'
	oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
	sharpcorners, thick=3
	multiplot
	plot, cat.dec, (cat.dec - cat.xdec)*3600d, xrange=minmax(cat.dec), yrange=yrange, xstyle=1, ystyle=1, psym=8, xtitle='DEC [degrees]'
	oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
	sharpcorners, thick=3

	if (n_elements(ps) gt 0) then psclose
	cleanplot, /silent
END

;-------------------------------------------------------------------------------------------

PRO supcam_mos_astrom_writeastrom, filenames, nfiles, astr, post=post
	COMPILE_OPT idl2, HIDDEN 
	if (n_elements(post) gt 0) then begin
		if (post eq '') then post = 'aS'
	endif else post = 'aS'

	for i=0,(nfiles-1) do begin
	
		img = mrdfits(filenames[i], 0, hdr, /silent)
	
		sxaddpar, hdr, 'CRVAL1', astr.crval1, /saveComment
		sxaddpar, hdr, 'CRVAL2', astr.crval2, /saveComment
		sxaddpar, hdr, 'CDELT1', astr.cdelt1, /saveComment
		sxaddpar, hdr, 'CDELT2', astr.cdelt2, /saveComment
		sxaddpar, hdr, 'CRPIX1', astr.crpix1[i], /saveComment
		sxaddpar, hdr, 'CRPIX2', astr.crpix2[i], /saveComment
		sxaddpar, hdr, 'CROTA2', astr.crota2[i]*!radeg, /saveComment, BEFORE='HISTORY'
	
		sxdelpar, hdr, 'CROTA'
		sxdelpar, hdr, 'CD1_1'
		sxdelpar, hdr, 'CD1_2'
		sxdelpar, hdr, 'CD2_1'
		sxdelpar, hdr, 'CD2_2'
		sxdelpar, hdr, 'PC001001'
		sxdelpar, hdr, 'PC001002'
		sxdelpar, hdr, 'PC002001'
		sxdelpar, hdr, 'PC002002'		
		
		ipos = strpos(filenames[i],'.gz')
		if (ipos gt -1) then outFile = post + strtrim(strmid(filenames[i],0,ipos),2)
		mwrfits, img, outFile, hdr
		;modfits, 'as'+filenames[i], 0, hdr
		;test = data.(i)
		;xyad, hdr, test.x, test.y, a1, d1
		;!p.multi=[0,1,2]
		;plothist, (test.ra - a1)*3600d, bin=0.1
		;oplot, intarr(10), indgen(10)*100, linestyle=2
		;plothist, (test.dec - d1)*3600d, bin=0.1
		;oplot, intarr(10), indgen(10)*100, linestyle=2
	
	endfor
	astr = struct_addtags(astr, {file:post+filenames})
END

;-------------------------------------------------------------------------------------------

FUNCTION supcam_mos_astrom, data, date, filenames=filenames, maxiter=maxiter, cat=cat, start=start, distcorr=distcorr, post=post
	COMPILE_OPT idl2, HIDDEN 
	if (n_elements(maxiter) eq 0) then maxiter = 2
	common sharemos
	
	nsig = 2.0d
	
	data0 = data

	if (n_elements(distcorr) gt 0) then begin
		for i=0,(n_tags(data)-1) do begin
			hdr = headfits(filenames[i])
			chip = strtrim(sxpar(hdr, 'DETECTOR'),2)
			supcam_distcorr, data.(i).x, data.(i).y, chip, date, x1, y1
			data0.(i).x = x1
			data0.(i).y = y1
		endfor
		;stop
	endif
	
	;initial guesses at the offsets and rotations
	
	if (n_elements(start) eq 0) then begin
		CDELT1 = -5.6093104d-05							;initial guess at RA scale (deg/pixel)
		CDELT2 = 5.6095829d-05
		;CRVAL1 = 151.29117d								;intial tangent projection RA
		;CRVAL2 = -7.7264511d								;initial tangent projection DEC
	
		for i=0,(n_tags(data)-1) do dhold = struct_append(dhold, data.(i))
		CRVAL1 = median(dhold.RA)
		CRVAL2 = median(dhold.DEC)
		undefine, dhold
	
		;CRPIX1 = 5342 - delx - min(delx)						;x-distance to tangent projection x
		;CRPIX2 = 3176 - dely - min(dely)						;y-distance to tangent projection y
		;CROTA2 = delx*0d - 0.0043
	
		if (date lt 2008.58) then begin
			if (strmid(filenames[0],0,1) eq 'n') then begin
				CRPIX1 = [5372.53,3285.94,1162.94,-956.959,-3044.84,5354.53,3269.39,1156.98,-966.578,-3053.81]	
				CRPIX2 = [-49.7044,-77.3883,-83.3438,-76.6404,-91.2516,4005.36,4019.26,4018.71,4004.06,3979.33]
				CROTA2 = -1*[-6.0914919e-03,-2.0648508,-4.9069205,4.2234785,1.3125425,1.5521885,4.2652558,4.6701419,-1.6867016e-04,0.71141694]	
			endif else begin
				CRPIX1 = [5246.5142,3165.0548,1054.4712,-1066.5537,-3161.3642,5244.9220,3163.8577,1051.1501,-1065.2115,-3159.0162]
				CRPIX2 = [35.490705,31.809698,33.976075,32.933892,30.762235,4105.4282,4120.2479,4118.9465,4118.8147,4102.6354]
				CROTA2 = -1*[4.5489500,4.8252375,4.3339790,4.5679359,4.3634360,4.0064291,4.0187368,3.9573857,4.4091170,4.6163363]*1d-3
			endelse
		endif else begin
			if (strmid(filenames[0],0,1) eq 'n') then begin
			endif else begin
			endelse
			CRPIX1 = [5425.8691,3351.0452,1241.0316,-871.95972,-2964.5720,5423.8774,3351.3420,1242.9849,-877.50110,-2966.1489]
			CRPIX2 = [111.32944,111.39816,111.15529,113.26266,113.70063,4315.5620,4330.3096,4338.3032,4330.0840,4318.5273]
			CROTA2 = [-0.0055500342,-0.0053783910,-0.0054863331,-0.0053572311,-0.0051588219,-0.0050752364,-0.0053784340,-0.0049788447,-0.0052406706,-0.0055207304]
		endelse
	endif else begin
		CDELT1 = start.CDELT1
		CDELT2 = start.CDELT2
		CRVAL1 = start.CRVAL1
		CRVAL2 = start.CRVAL2
		CRPIX1 = start.CRPIX1
		CRPIX2 = start.CRPIX2
		CROTA2 = start.CROTA2
	endelse
	CRVAL2_0 = CRVAL2

	
	;arg = [CDELT1, CDELT2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CROTA2]
	arg = [CDELT1, CDELT2, CRVAL1, CRVAL2]*1d
	
	k = 0
	del = arg
	while 1 do begin
	
		if (k eq 0) then ftol = 0.01 else ftol = 0.01
				
		counter = 0
		
		;dummy = supcam_mos_astrom_newtify(del)
	
		del = amoeba(ftol, function_name='supcam_mos_astrom_newtify', function_value=val, p0=del, scale=[1d-8,1d-8,1d-3,1d-3])
		if (n_elements(del) eq 1) then stop
	
		cat = supcam_mos_astrom_install(del)										;install the guesses into the catalog
		if (k eq 0) then count0 = n_elements(cat)
	
		print, robust_sigma(cat.ra - cat.xra)*3600d
		print, robust_sigma(cat.dec - cat.xdec)*3600d
		dev = sqrt( (cat.ra - cat.xra)^2d + (cat.dec - cat.xdec)^2d )*3600d			;rms value of the subset catalog (if k > 0)
	
		if (n_elements(flag) gt 0) then	undefine, flag								;undefine flag if its set
		cat = supcam_mos_astrom_install(del)										;get the whole catalog
		dev_all = sqrt( (cat.ra - cat.xra)^2d + (cat.dec - cat.xdec)^2d )*3600d		;rms value of the whole catalog
	
	
		window, 1
		cleanplot, /silent
		plot, cat.ra, dev_all, psym=1, title='Iter # '+strtrim(k+1,2), xrange=reverse(minmax(cat.ra))
		
		gind = -1
		for t=0,(n_tags(data)-1) do begin
			ind = where(cat.chip eq t,nMatch)
			if (nMatch eq 0) then begin
				print, 'there is something wrong with the photometry on this chip.'
				stop
			endif
			tgind = where(dev_all[ind] lt nsig*robust_sigma(dev[ind], /zero), nGood);isolate the good objects, could do this on a chip to chip basis
			if (nGood eq 0) then stop
			gind = [gInd, ind[tgind]]												;keep the indices of the good objects
		endfor
		gInd = gind[1:*]
		oplot, cat[gind].ra, dev_all[gind], psym=1, color=fsc_color('red')
		wset, 0
		count = n_elements(gind)	
		if (count eq count0) then break												;if count has not changed, then quit iterating
		count0 = count
		flag = 1
		++k
		if (k ge maxiter) then break	
	endwhile
	
	astr = {cdelt1:del[0]*1d, cdelt2:del[1]*1d, crval1:del[2]*1d, crval2:del[3]*1d, crpix1:crpix1*1d, crpix2:crpix2*1d, crota2:crota2*1d}
	
	if ((nfiles = n_elements(filenames)) gt 0) then supcam_mos_astrom_writeastrom, filenames, nfiles, astr, post=post
	
	supcam_mos_astrom_plot_resid, cat, /ps, outfile='astrometry_'+strmid(filenames[0],0,(strsplit(filenames[0],'_'))[-1]-1)
	
	dummy = file_search('dir_astrometry', /test_directory, count=nMatch)
	if ~nMatch then spawn, 'mkdir dir_astrometry'
	spawn, 'mv astrometry*eps dir_astrometry/'
	
	dellist = ['data0', 'gind', 'crpix1', 'crpix2', 'crota2', 'counter', 'flag', 'results']
	dellist = [dellist, ['dum', 'crval1', 'crval2', 'cdelt1', 'cdelt2', 'res']]
	for i=0,(n_elements(dellist)-1) do dummy = execute('undefine, '+dellist[i])
	
	return, astr

END



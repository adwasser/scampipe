
;-------------------------------------------------------------------------------------------

FUNCTION queryusno4, ra, dec, radius, epoch=epoch
	if (n_elements(epoch) eq 0) then epoch = 2000.0
	maxErrPos	= 500	;0.5 arcseconds				;320d			;mas
	maxPM	= 40			;mas/yr
	maxMupr		= 8				;60%
	maxErrFit	= 5d			;4 = 400 mas = 0.4 arcseconds
	
	;----------get/trim initial catalog of objects---------------------------------------------
	print
	print, 'Querying USNO-B 1.0 ...'
	; For some reason the canada tag was making it impossible to connect	
	astro = queryvizier('USNO-B1', [ra, dec], radius, /AllColumns)				;query USNO B 1.0
	;astro = queryvizier('USNO-B1', [ra, dec], radius, /AllColumns, /canada)				;query USNO B 1.0
	if (n_elements(astro) eq 1) then astro = queryvizier('USNO-B1', [ra, dec], radius, /AllColumns)				;query USNO B 1.0
	
	print, strtrim(n_elements(astro),1) + ' objects returned.'
	
	astro.raj2000 = astro.raj2000 + astro.pmra/1000d/3600d*(epoch - 2000.)										;precess the RA
	astro.dej2000 = astro.dej2000 + astro.pmde/1000d/3600d*(epoch - 2000.)										;	and DEC
	astro.e_raj2000 = round(sqrt((astro.e_raj2000)^2d + (astro.e_pmra*(epoch - astro.epoch))^2d))			;add the positional and proper motion errors in quadrature, convert to 
	astro.e_dej2000 = round(sqrt((astro.e_dej2000)^2d + (astro.e_pmde*(epoch - astro.epoch))^2d))
	
	
	astro = struct_addtags(astro, replicate({errpos:0d, pm:0d, errfit:0d}, n_elements(astro)))
	;---make preliminary cuts
	astro.errPos = sqrt( (astro.E_RAJ2000)^2d +  (astro.E_DEJ2000)^2d)	;Error on Position of star (mas)
	astro.errFit = sqrt( (astro.FIT_RA)^2d + (astro.FIT_DE)^2d)			;Error in the fitting of the RA and DEC, not sure exactly what this is
	astro.pm	= sqrt(astro.pmra^2d + astro.pmde^2d)					;total proper motion
	
	astrom = astro[where( ((finite(astro.B1mag) ne 0) or (finite(astro.R1mag) ne 0)) and ((finite(astro.B2mag) ne 0) or (finite(astro.R2mag) ne 0) or (finite(astro.IMAG) ne 0)) )]		;require 2 epochs of observation
	;astrom = astrom[where( (abs(astrom.pm) lt maxPM) )]																									;no more than 15 mas/yr
	
	astrom = astrom[where( (astrom.mupr ge 6))]				;greater than 60% probability of proper mtion being correct
	astrom = astrom[where(astrom.errPos le maxErrPos)]		;position is known to better than 0.5", this includes proper motion uncertainty
	astrom = astrom[where(astrom.errfit le maxErrFit)]		;this generally does the same thing as the mupr restriction
	astrom = astrom[where(astrom.pm le maxpm)]				;get rid of large proper motion stars to limit potential problems
	
	;	openw, lun, 'astrom.reg', /get_lun
	;	for i=0,(n_elements(astrom)-1) do printf, lun, astrom[i].raj2000, astrom[i].dej2000, format='(2(F-20.9))'
	;	free_lun, lun
	
	return, astrom

END

;-------------------------------------------------------------------------------------------

FUNCTION query2mass, ra, dec, radius
	
	;----------get/trim initial catalog of objects---------------------------------------------
	print
	print, 'Querying 2MASS ...'
	
	;ra, dec, radius
	;stop
	;astro = queryvizier('2MASS-PSC', [ra, dec], radius, /AllColumns)				;query USNO B 1.0
	
	astro = queryvizier('2MASS-PSC', [ra, dec], radius, /canada)
	if (type(astro,/string) eq 'longint') then return, -1l
	
	print, strtrim(n_elements(astro),1) + ' objects returned.'
	
	astro = astro[where( fix(astro.rflg) le 222) ]									;remove any saturated stars
	astro = astro[where( fix(astro.bflg) le 111) ]									;remove blended objects
	astro = astro[where( astro.cflg eq '000') ]										;remove bad objects
	astro = astro[where( astro.xflg eq 0) ]											;is a point-source
	astro = astro[where( astro.aflg eq 0) ]											;not associated with a solar system object
	
	return, astro

END

;-------------------------------------------------------------------------------------------

PRO makeSumPlot, data, info, astropar
	use = astropar.use
	info.hdr = headfits(info.image)
	adxy, info.hdr, data[use].ra, data[use].dec, newx, newy


	erase
	cleanplot, /silent
	psopen, 'diag_astrometry_0'+info.catalog, /encapsulated, /color, /heiles, xs=8, ys=6, /inches
	!p.multi=[0,1,2]
	!p.charsize = 1.3
	!p.charthick = 3
	multiplot
	plothist, (xval = (data[use].x - newx)*astropar.xscale), bin=0.05, xrange=[-1,1]*1, ytitle='Y'
	legend, 'x rms = ' + roundx(robust_sigma(xval),2), /top, /right
	sharpcorners, thick=3
	multiplot
	plothist, (yval = (data[use].y - newy)*astropar.yscale), bin=0.05, xrange=[-1,1]*1, ytitle='X', xtitle='RMS [arcseconds]'
	legend, 'y rms = ' + roundx(robust_sigma(yval),2), /top, /right
	sharpcorners, thick=3
	psclose

END

;-------------------------------------------------------------------------------------------

FUNCTION doRoughAstrom, imageName, make, info, rough=rough, astrom=astrom, astropar=astropar

	spawn, 'ls ' + imageName, result
	if ~(result) then stop	;print, 'Image found.' else stop
	img = mrdfits(imageName, 0, hdr, /silent)					;read in image plus header
	
	info = struct_addtags(info, { $
		xsize:sxpar(hdr, 'NAXIS1') - 1 , $			;width of image
		ysize:sxpar(hdr, 'NAXIS2') - 1 })			;height of image
	

	if (info.xsize lt 3000) then begin
		info.rows = 6
		info.cols = 1
	endif
	
	if (make eq 'yes') then begin
	
		info = struct_addtags(info, { $
			RA: sxpar(hdr, 'CRVAL1'), $
			DEC : sxpar(hdr, 'CRVAL2')})
	
		if (n_elements(batch) eq 0) then begin
	
			astrom = queryusno4(info.ra, info.dec, info.radius)



		
			print
			print, 'You will cycle through ' + strtrim(round(info.rows*info.cols),1) + ' frames to do a quick fix to the WCS.'
			print, 'First click on the green square, then click on the star that it should belong to.'
	
			xInc = round(info.xSize/info.cols)										;number of steps in x for the crude astrometry correction
			yInc = round(info.ySize/info.rows)										;number of steps in y for the crude astrometry correction
	
	
			list = replicate({x:0d, y:0d, ra:0d, dec:0d, era:0d, edec:0d}, round(info.rows*info.cols))	;
	
	
			window, 0, xsize=1000, ysize=1000
			cleanplot, /silent
			erase
			!p.multi=[0,1,1]
			loadct, 0, /silent
			minY = 0
			maxY = yInc
			xoff = 0
			yoff = 0
			k = 0
			for j=0,(info.rows - 1) do begin
	
				minX = 0
				maxX = xInc
	
				for i=0,(info.cols - 1) do begin
		
					subImg = img[minX:maxX, minY:maxY]						;subImage to display
	
					;display, subImg, top=254, /aspect, min=median(img)-100, max=median(img)+1000, /silent		;display subImage
					disp_med = median(subimg[(idisp = where((subimg ne -32768) and (finite(subimg) eq 1)))])
					disp_sig = robust_sigma(subimg[idisp])
					jdisplay, subimg, min=disp_med - 10.*disp_sig, max=disp_med + 50.*disp_sig, top=254, /aspect, /silent
	
					adxy, hdr, astrom.raj2000, astrom.dej2000, xcat, ycat
					xcat = xcat - minx + xoff
					ycat = ycat - miny + yoff


					cleanplot, /silent
					erase
					jdisplay, subimg, min=disp_med - 10.*disp_sig, max=disp_med + 50.*disp_sig, top=254, /aspect, /silent

					oplot, xCat, yCat, psym=6, color=fsc_color('green')						;overplot catalog objects




	
					cursor, xCCat, yCCat, /up													;select catalog object choice
	
					dummy = min(((xCat - xCCat)^2d + (yCat - yCCat)^2d), ind)		;closest catalog object to cursor location
	
					cursor, xPix, yPix, /up								;select location of corresponding star

					gcntrd, img, xPix, yPix, xCen, yCen, 5			;centroid star
	
					if (xcen eq -1) then xcen = xPix				;if centroid fails, adopt initial guess
					if (ycen eq -1) then ycen = yPix				
	
					xoff = xoff + (xcen - xcat[ind])
					yoff = yoff + (ycen - ycat[ind])
	
					xPix = xCen + minX								;convert back to mosaic coordinates
					yPix = yCen + minY
			
					minX = minX + xInc								;increment minimum and maximum x-pixels for next subImage
					maxX = maxX + xInc
		
					maxX = min([maxX, info.xSize])						;restrict to actual width of mosaic
			
					list[k].x 	= xPix								;store x location of selected object
					list[k].y 	= yPix								;store y location of selected object
					list[k].ra 	= astrom[ind].raj2000				;store ra of selected object
					list[k].dec = astrom[ind].dej2000				;store dec of selected object
					list[k].era = astrom[ind].e_raj2000				;store RA error
					list[k].edec = astrom[ind].e_dej2000			;store DEC error
					
					k++
				endfor
	
				minY = minY + yInc											;increment minimum and maximum y-pixels for next subImage
				maxY = maxY + yInc
		
				maxY = min([maxY, info.ySize])									;restrict to actual height of mosaic
	
			endfor
			undefine, subImg
			undefine, img
	
	
			iniFile = 'iniWCS.reg'									;name of initial WCS file containing x,y,ra,dec
	
			file_delete, iniFile, /allow_nonexistent
	
	
			openw, lun, iniFile, /get_lun							;save objects to file
				printf, lun, [transpose(list.x), transpose(list.y), transpose(list.ra), transpose(list.dec)]
			free_lun, lun
		
			str = strsplit(imageName, '.fits', /extract, /regex)
			newImage = info.pre + str[0] + '.fits'
			if (n_elements(str) gt 1) then newImage = newImage + str[-1]
		
			if (n_elements(overwrite) eq 0) then begin
				if (imageName ne NewImage) then begin
					file_delete, newImage, /allow_nonexistent
					spawn, 'cp ' + imageName + ' ' + newImage
				endif
			endif

			;jaa_astrom, list.x, list.y, list.ra, list.dec, newImage, error, nsig=2.5, maxiter=10	;[[list.era],[list.edec]]/1000d/3600d		
			jaa_astrom1, list.x, list.y, list.ra, list.dec, newImage, error, nsig=2.5, maxiter=10, result=astropar	;[[list.era],[list.edec]]/1000d/3600d
				
			hdrNew = headfits(newImage)
	
			;astrom = astro[use]
	
		endif else begin
			;newImage = pre + (strsplit(image, '.fits', /extract, /regex))[0] + '.fits'		
			str = strsplit(imageName, '.fits', /extract, /regex)
			newImage = pre + str[0] + '.fits'
			if (n_elements(str) gt 1) then newImage = newImage + str[-1]
			if (n_elements(overwrite) eq 0) then begin
				if (imageName ne NewImage) then begin
					file_delete, newImage, /allow_nonexistent
					spawn, 'cp ' + imageName + ' ' + newImage
				endif
			endif
			hdrNew = hdr
		endelse
	
		if (n_elements(rough) ne 0) then stop	;return			;this needs to return us to the calling program, so needs to be fixed

		catflag = 1
		while catflag do begin
			case strupcase(info.catalog) of
				'SDSS' 	:	BEGIN
								print, 'Querying SDSS...'
								sdss = queryvizier('II/294', [info.ra, info.dec], info.radius, /AllColumns)
								if (type(sdss,/string) eq 'longint') then begin
									print, 'No SDSS information, using USNO catalog for astrometry.'
									info.catalog = 'USNO'
									break
								endif else catflag = 0
								print, strtrim(n_elements(sdss),1) + ' objects returned.'
		
								use = where( (sdss.zs eq 1))			;classified as a star
								sdss = sdss[use]
								astrom = sdss
								astrom = astrom[where( astrom.gmag le 24)]
								astrom = astrom[where( astrom.gmag ge 16)]
							END
				'USNO' 	:	BEGIN
								;catalog is already USNO
								catflag = 0
							END
				'2MASS' :	BEGIN
								astrom = query2mass(info.ra, info.dec, radius)
								if (type(astrom,/string) eq -1) then begin
									info.catalog = 'USNO'
									break
								endif else catflag = 0
							END
			endcase
		endwhile
	endif else begin
		hdrnew = hdr
		newimage = imageName
	endelse
	info = struct_addtags(info, {image:newImage, hdr:hdrnew})								;add the new image name and header to the info structure
	
	adxy, hdrnew, astrom.raj2000, astrom.dej2000, xnew, ynew

	list = arr_struct({x:xnew, y:ynew, ra:astrom.raj2000, dec:astrom.dej2000});, era:astrom.e_raj2000, edec:astrom.e_dej2000})
	list = list[where( (xnew ge 0) and (xnew le info.xSize) and (ynew ge 0) and (ynew le info.ySize))]

	openw, lun, 'list.reg', /get_lun
		for i=0,(n_elements(list)-1) do printf, lun, list[i].x, list[i].y, format='(2(F-20.9))'
	free_lun, lun

	return, list
END

;-------------------------------------------------------------------------------------------

FUNCTION centroidObjects, list, info, astropar
	print
	;print, list
	print, 'Centroiding ' + strtrim(n_elements(list),1) + ' objects on image '+strtrim(info.image)


	d = {ra:list.ra, dec:list.dec, raj2000:list.ra, dej2000:list.dec}

	;print, d.ra, d.dec

	objects = runse(info.image, zp=29., exptime=1d, /batch, /cleanup, objectPars=objectPars, config=config)
	;plot, objects.mag_aper1, psym=1

	dd = seSelection(objects, {magcut:[10d, 28d], CLASS_STAR:0d, fwhmcut:[95d, 95d], apercut:[95d, 95d], sncut:[1d, 1d4]}, /simple, /batch)
	dd = dd[where(dd.ellipticity lt 0.5)]
	dd = dd[where(dd.flags lt 2)]	;gets rid of objects with possible bad centroiding
	sn = 1d/dd.magerr_aper1
	dd = dd[where( (sn gt 20d) )]	;signal to noise cut

	;match centroided objects to astrometric objects
	mDist = 10d
	;print, list
	img_DEC = median(list.DEC)*!dtor
	match = match_2d(dd.x_world*cos(img_DEC), dd.y_world, list.ra*cos(img_DEC), list.dec, mDist*2d/3600d, match_distance=minDist)
	iMatch = where( (match ne -1) and (minDist*3600d le mDist), nMatch)
	minDist = minDist[iMatch]
	ddd = struct_addtags(dd[iMatch], list[match[iMatch]])
	;plot, (ddd.x_world - ddd.ra)*3600d, (ddd.y_world - ddd.dec)*3600d, psym=1, /iso
	
	RA_offset = median(ddd.x_world - ddd.ra)
	DEC_offset = median(ddd.y_world - ddd.dec)

	;match the centroided objects to astrometric objects allowing for a bulk offset
	mDist = 1d	;config.SEEING_FWHM
	if (n_elements(astropar) gt 0) then mDist = (sqrt(robust_sigma(astropar.radev)^2d + robust_sigma(astropar.decdev)^2d) > mDist)	;search radius is 3x the greater of
	mDist = mDist*3d < 4d														;the rough astrometry fitting residuals
																												;or the seeing, but no greater than 4 arcsec
	img_DEC = median(list.DEC)*!dtor
	match = match_2d((dd.x_world - RA_offset)*cos(img_DEC), dd.y_world - DEC_offset, list.ra*cos(img_DEC), list.dec, 5d/3600d, match_distance=minDist)
	iMatch = where( (match ne -1) and (minDist*3600d le mDist), nMatch)
	minDist = minDist[iMatch]
	ddd = struct_addtags(dd[iMatch], list[match[iMatch]])
	;plot, (ddd.x_world - ddd.ra)*3600d, (ddd.y_world - ddd.dec)*3600d, psym=1, /iso
	;help, ddd
	
	;plot, minDist*3600d, psym=1
	
	d = ddd
	d.x = d.x_image
	d.y = d.y_image
	d = struct_trimtags(d, select=['X','Y','RA','DEC'])
	d = struct_addtags(d, arr_struct({raj2000:d.ra, dej2000:d.dec}))
	
	data = d
	print, strtrim(n_elements(data),1) + ' objects centroided!'

	return, data




	openw, lun, '0test.reg', /get_lun
	printf, lun, [transpose(d.x), transpose(d.y)]
	free_lun, lun
	
	pause
	
	


	;old way	
	stop

	img0 = mrdfits(info.image,/silent)							;read in the image array
	img = img0
	iBad = where( (img le -32768), nBad)						;bad pixels will be set to -32768
	if (nBad gt 0) then img[iBad] = !values.f_nan				;replace bad pixels with NaNs
	xcen = list.x												;input x-coordinates to be centroided
	ycen = list.y												;input y-coordinates to be centroided
	cBox = 20													;centering box size
	objrad = 5													;aperture size
	sAlg = 'median'												;sky estimation method
	skyrad = [15,25]											;inner and outer radius of sky annulus
	
	d = {xin:list.x, yin:list.y, ra:list.ra, dec:list.dec, raj2000:list.ra, dej2000:list.dec}
	
	flux = djs_phot(xcen, ycen, objrad, skyrad, img, cbox=cbox, salg=salg, flerr=flerr, $	;centroid the objects
		skyval=skyval, peakval=peakval)
		
	d = arr_struct(struct_addtags(d, {x:xcen, y:ycen, flux:flux, flerr:flerr, skyval:skyval, peakval:peakval}))

	iBad = where( $
		~finite(d.flux) $									;the flux is bad
		or ~finite(d.flerr) $								;the error on the flux is bad
		or ~finite(d.skyval) $								;the sky is bad
		or (d.peakval ge info.datamax) $					;there are saturated pixels in the aperture
		or (sqrt((d.x - d.xin)^2d + (d.y - d.yin)^2d) ge 20d) $	;the centroided position is too far away
		, nBad, complement=iGood, nComplement=nGood)		;

	if (nGood lt 10) then stop	;need more good objects than this to get a good solution
	d = d[iGood]
	
	d.xin = (xcen = d.x)
	d.yin = (ycen = d.y)
	cbox = 7												;now use a much smaller centroiding box
	
		
	flux = djs_phot(xcen, ycen, objrad, skyrad, img, cbox=cbox, salg=salg, flerr=flerr, $	;centroid the objects
		skyval=skyval, peakval=peakval)
	d.x = xcen
	d.y = ycen
	d.flux = flux
	d.flerr = flerr
	d.skyval = skyval
	d.peakval = peakval
	
	iBad = where( $
		~finite(d.flux) $									;the flux is bad
		or ~finite(d.flerr) $								;the error on the flux is bad
		or ~finite(skyval) $								;the sky is bad
		or (d.peakval ge info.datamax) $					;there are saturated pixels in the aperture
		or (sqrt((d.x - d.xin)^2d + (d.y - d.yin)^2d) ge 20d) $	;the centroided position is too far away
		, nBad, complement=iGood, nComplement=nGood)		;
	
	if (nGood lt 10) then stop	;need more good objects than this to get a good solution
	d = d[iGood]				;keep just the objects that were most likely to be centroided well
	++d.x		;change back to 1 indexing
	++d.y		;change back to 1 indexing


	openw, lun, '0test.reg', /get_lun
	printf, lun, [transpose(d.x), transpose(d.y)]
	free_lun, lun
	
	
	;plot, d.flux, -2.5*alog10(d.flux), psym=1
	
	ind = where(-2.5*alog10(d.flux) + 29. lt 24,num)
	dd = d[ind]
	openw, lun, '1test.reg', /get_lun
	printf, lun, [transpose(dd.x), transpose(dd.y)]
	free_lun, lun
		
	
	
	
	stop







	data = d
	print, strtrim(n_elements(data),1) + ' objects centroided!'

	return, data
END

;-------------------------------------------------------------------------------------------

PRO astrometry, image, object=object, pre=pre, overwrite=overwrite, batch=batch, nsig=nsig, $
		maxiter=maxiter, astropar=astropar, catalog=catalog, datamax=datamax, nthreads=nthreads, $
		rough=rough, bridge=bridge, nodestroy=nodestroy, centroid=centroid
	COMPILE_OPT idl2, HIDDEN 
	
	if (n_elements(image) eq 0) then image = choose('Specify image:', 'mos.fits')
	if (n_elements(pre) eq 0) then pre = 'astrom_'
	if (n_elements(overwrite) eq 1) then pre = ''
	if (n_elements(nsig) eq 0) then nsig = 2.
	if (n_elements(maxiter) eq 0) then maxiter = 10
	if (n_elements(catalog) eq 0) then catalog = 'sdss'
	if (n_elements(datamax) eq 0) then datamax = 180000.
	
	info = {radius:[50,40], rows:2., cols:3., catalog:catalog, datamax:datamax[0], pre:pre}			; ask radius  (longitude, latitude)
	
	make = 'yes'
	if (n_elements(object) gt 0) then begin
		hdu = 1
		fits_info, image, /silent, n_ext=n_ext
		for hdu=1,n_ext do begin
			dummy = mrdfits(image, hdu, /silent)
			if (type(dummy, /string) eq 'structure') then begin
				if tag_exist(dummy, 'raj2000') then begin
					astrom = dummy
					make = 'no'
				endif
				dhdu = create_struct('num'+strtrim(hdu,2), dummy)
				if ((hdu eq 1) and ~tag_exist(dummy, 'raj2000')) then hdus = dhdu else hdus = struct_addtags(hdus, dhdu)
			endif else break
		endfor
	endif
	;if usno then make = 'yes'
	
; Trying to keep the rough astrometry step. - Alexa	
	list = doRoughAstrom(image, make, info, rough=rough, astrom=astrom, astropar=astropar)	;make a rough astrometric solution

	data = centroidObjects(list, info, astropar)											;centroid the catalog objects on the image

	jaa_astrom1, data.x, data.y, data.ra, data.dec, info.image, error, $				;calibrate the astrometry
		nsig=nsig, maxiter=10, result=astropar, /quiet
	
	
	;makeSumPlot, data, info, astropar													;summary plot showing how good the astrometry came out
	

	if (n_elements(object) gt 0) then begin
		img = mrdfits(info.image, 0, info.hdr, /silent)					;read in image plus header
		file_delete, info.image+'.old', /allow_nonexistent
		dummy = file_search(info.image,count=nMatch)
		if (nMatch gt 0) then spawn, 'mv ' + info.image + ' ' + info.image + '.old'
		mwrfits, img, info.image, info.hdr					;write the image file and updated header
		if (n_elements(hdus) gt 0) then for i=0,(n_elements(hdus) - 1) do  mwrfits, hdus.(i), info.image, /silent	;write each of the hdus

		for tag=0,(n_tags(astrom)-1) do begin
			tagType = strupcase(type(astrom[0].(tag),/string))
			if (tagType ne 'STRING') then begin
				iBad = where(~finite(astrom.(tag)),nBad)
				if (nBad gt 0) then astrom[iBad].(tag) = -1
			endif
		endfor

		mwrfits, astrom, info.image, /silent								;write the unculled object list
		hdr = headfits(info.image)
		sxaddpar, hdr, 'NEXTENS', n_elements(hdus)+1, 'Number of extensions.', /savecomment
		modFits, info.image, 0, hdr
		file_delete, info.image+'.old', /allow_nonexistent
	endif
	centroid = data						;all the centroided points

end

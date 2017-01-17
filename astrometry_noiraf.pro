
function queryusno1, ra, dec, radius

	;----------get/trim initial catalog of objects---------------------------------------------
	print
	print, 'Querying USNO-B 1.0 ...'
	
	;ra, dec, radius
	
	astro = queryvizier('USNO-B1', [ra, dec], radius, /AllColumns)				;query USNO B 1.0
	print, strtrim(n_elements(astro),1) + ' objects returned.'
	
	maxErrPos	= 320d			;mas
	maxErrPM	= 15			;mas/yr
	maxMupr		= 8				;60%
	maxErrFit	= 5d			;4 = 400 mas = 0.4 arcseconds
	
	;---make preliminary cuts
	errPos = sqrt( (astro.E_RAJ2000)^2d +  (astro.E_DEJ2000)^2d)	;Error on Position of star (mas)
	errPM  = sqrt( (astro.E_PMRA)^2d + (astro.E_PMDE)^2d)			;Error on proper motion of star
	;mupr															;Probability of proper motion being correct
	errFit = sqrt( (astro.FIT_RA)^2d + (astro.FIT_DE)^2d)			;Error in the fitting of the RA and DEC
	
	use = where((finite(astro.B1mag) ne 0) and (finite(astro.B2mag ne 0)) and $ 	;Require at least 2 epochs of observation 
				(errPos le maxErrPos) and (errPM le maxErrPM) and $	
				(astro.mupr ge maxMupr) and (errFit le maxErrFit)	)
	
	astrom = astro[use]
	;------------------------------------------------------------------------------------------
	return, astrom

end



function queryusno, ra, dec, radius
	epoch = 2008.91
	;----------get/trim initial catalog of objects---------------------------------------------
	print
	print, 'Querying USNO-B 1.0 ...'
	
	;ra, dec, radius
	
	astro = queryvizier('USNO-B1', [ra, dec], radius, /AllColumns)				;query USNO B 1.0
	print, strtrim(n_elements(astro),1) + ' objects returned.'
	
	astro.raj2000 = astro.raj2000 + astro.pmra/1000d/3600d*(epoch - 2000.)										;precess the RA
	astro.dej2000 = astro.dej2000 + astro.pmde/1000d/3600d*(epoch - 2000.)										;	and DEC
	
	
	maxErrPos	= 500d	;320d			;mas
	maxPM	= 15			;mas/yr
	maxMupr		= 8				;60%
	maxErrFit	= 5d			;4 = 400 mas = 0.4 arcseconds
	
	;---make preliminary cuts
	astro = struct_addtags(astro, replicate({errpos:0d, pm:0d, errfit:0d}, n_elements(astro)))
	astro.errPos = sqrt( (astro.E_RAJ2000)^2d +  (astro.E_DEJ2000)^2d)	;Error on Position of star (mas)
	astro.PM  = sqrt( (astro.PMRA)^2d + (astro.PMDE)^2d)			;proper motion of star
	;mupr															;Probability of proper motion being correct
	astro.errFit = sqrt( (astro.FIT_RA)^2d + (astro.FIT_DE)^2d)			;Error in the fitting of the RA and DEC
	
	
	astrom = astro[where( ((finite(astro.B1mag) ne 0) or (finite(astro.R1mag ne 0)) and (finite(astro.B2mag) ne 0) or (finite(astro.R2mag ne 0)) ) )]		;require 2 epochs of observation
	astrom = astrom[where( (abs(astrom.pm) lt maxPM) )]																									;no more than 15 mas/yr
	
	
	; and $ 	;Require at least 2 epochs of observation 
	;			(errPos le maxErrPos) and (errPM le maxErrPM) and $	
	;			(astro.mupr ge maxMupr) and (errFit le maxErrFit)	)
	
	;astrom = astro[use]
	;------------------------------------------------------------------------------------------
	return, astrom

end












function queryusno3, ra, dec, radius, epoch=epoch
	epoch = 2008.91
	if (n_elements(epoch) eq 0) then epoch = 2000.0
	maxErrPos	= 500	;0.5 arcseconds				;320d			;mas
	maxPM	= 15			;mas/yr
	maxMupr		= 8				;60%
	maxErrFit	= 5d			;4 = 400 mas = 0.4 arcseconds
	
	;----------get/trim initial catalog of objects---------------------------------------------
	print
	print, 'Querying USNO-B 1.0 ...'
	
	;ra, dec, radius
	
	astro = queryvizier('USNO-B1', [ra, dec], radius, /AllColumns, /canada)				;query USNO B 1.0
	print, strtrim(n_elements(astro),1) + ' objects returned.'
	
	astro.raj2000 = astro.raj2000 + astro.pmra/1000d/3600d*(epoch - 2000.)										;precess the RA
	astro.dej2000 = astro.dej2000 + astro.pmde/1000d/3600d*(epoch - 2000.)										;	and DEC
	astro.e_raj2000 = round(sqrt((astro.e_raj2000)^2d + (astro.e_pmra*(epoch - astro.epoch))^2d))			;add the positional and proper motion errors in quadrature, convert to 
	astro.e_dej2000 = round(sqrt((astro.e_dej2000)^2d + (astro.e_pmde*(epoch - astro.epoch))^2d))
	
	
	astro = struct_addtags(astro, replicate({errpos:0d, pm:0d, errfit:0d}, n_elements(astro)))
	;---make preliminary cuts
	astro.errPos = sqrt( (astro.E_RAJ2000)^2d +  (astro.E_DEJ2000)^2d)	;Error on Position of star (mas)
	astro.errFit = sqrt( (astro.FIT_RA)^2d + (astro.FIT_DE)^2d)			;Error in the fitting of the RA and DEC, not sure exactly what this is
	
	
	astrom = astro[where( ((finite(astro.B1mag) ne 0) or (finite(astro.R1mag) ne 0)) and ((finite(astro.B2mag) ne 0) or (finite(astro.R2mag) ne 0) or (finite(astro.IMAG) ne 0)) )]		;require 2 epochs of observation
	;astrom = astrom[where( (abs(astrom.pm) lt maxPM) )]																									;no more than 15 mas/yr
	
	astrom = astrom[where( (astrom.mupr ge 6))]				;greater than 60% probability of proper mtion being correct
	astrom = astrom[where(astrom.errPos le maxErrPos)]		;position is known to better than 0.5", this includes proper motion uncertainty
	astrom = astrom[where(astrom.errfit le maxErrFit)]
	
	;	openw, lun, 'astrom.reg', /get_lun
	;	for i=0,(n_elements(astrom)-1) do printf, lun, astrom[i].raj2000, astrom[i].dej2000, format='(2(F-20.9))'
	;	free_lun, lun
	;	openw, lun, 'astrom1.reg', /get_lun
	;	for i=0,(n_elements(astrom)-1) do printf, lun, astrom[i].raj2000, astrom[i].dej2000, format='(2(F-20.9))'
	;	free_lun, lun
	
	;stop
	
	; and $ 	;Require at least 2 epochs of observation 
	;			(errPos le maxErrPos) and (errPM le maxErrPM) and $	
	;			(astro.mupr ge maxMupr) and (errFit le maxErrFit)	)
	
	;astrom = astro[use]
	;------------------------------------------------------------------------------------------
	return, astrom

end








function queryusno4, ra, dec, radius, epoch=epoch
	if (n_elements(epoch) eq 0) then epoch = 2000.0
	maxErrPos	= 500	;0.5 arcseconds				;320d			;mas
	maxPM	= 40			;mas/yr
	maxMupr		= 8				;60%
	maxErrFit	= 5d			;4 = 400 mas = 0.4 arcseconds
	
	;----------get/trim initial catalog of objects---------------------------------------------
	print
	print, 'Querying USNO-B 1.0 ...'
		
	astro = queryvizier('USNO-B1', [ra, dec], radius, /AllColumns, /canada)				;query USNO B 1.0
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
	
	;------------------------------------------------------------------------------------------
	return, astrom

end








function query2mass, ra, dec, radius
	
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

end

























;want no astrom and usno to mean different things
;usno just forces usno, or 2mass forces 2mass, while 




pro astrometry, image, object=object, pre=pre, overwrite=overwrite, batch=batch, nsig=nsig, maxiter=maxiter, astropar=astropar, $
				catalog=catalog, datamax=datamax, nthreads=nthreads, rough=rough, bridge=bridge, nodestroy=nodestroy, centroid=centroid

	COMPILE_OPT idl2, HIDDEN 
	
	if (n_elements(image) eq 0) then image = choose('Specify image:', 'mos.fits')
	if (n_elements(pre) eq 0) then pre = 'astrom_'
	if (n_elements(overwrite) eq 1) then pre = ''
	if (n_elements(nsig) eq 0) then nsig = 2.
	if (n_elements(maxiter) eq 0) then maxiter = 10
	if (n_elements(catalog) eq 0) then catalog = 'sdss'
	if (n_elements(datamax) eq 0) then datamax = 180000.
	
	
	radius = [50,40]							; ask radius  (longitude, latitude)
	
	rows = 2*1.
	cols = 3*1.
	
	
	
	
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
	
	
	spawn, 'ls ' + image, result
	if ~(result) then stop	;print, 'Image found.' else stop
	img = mrdfits(image, 0, hdr, /silent)					;read in image plus header
	
	xsize = sxpar(hdr, 'NAXIS1') - 1				;width of image
	ysize = sxpar(hdr, 'NAXIS2') - 1				;height of image
	
	if (xsize lt 3000) then begin
		rows = 6
		cols = 1
	endif
	
	if (make eq 'yes') then begin
	
		ra = sxpar(hdr, 'CRVAL1')
		dec = sxpar(hdr, 'CRVAL2')
	
		if (n_elements(batch) eq 0) then begin
	
			astrom = queryusno4(ra, dec, radius)
			;print, N_elements(astrom)
			;astrom = query2mass(ra, dec, radius)
		
			print
			print, 'You will cycle through ' + strtrim(round(rows*cols),1) + ' frames to do a quick fix to the WCS.'
			print, 'First click on the green square, then click on the star that it should belong to.'
	
			xInc = round(xSize/cols)										;number of steps in x for the crude astrometry correction
			yInc = round(ySize/rows)										;number of steps in y for the crude astrometry correction
	
	
			list = replicate({x:0d, y:0d, ra:0d, dec:0d, era:0d, edec:0d}, round(rows*cols))	;
	
	
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
			for j=0,(rows - 1) do begin
	
				minX = 0
				maxX = xInc
	
				for i=0,(cols - 1) do begin
		
					subImg = img[minX:maxX, minY:maxY]						;subImage to display
	
					;display, subImg, top=254, /aspect, min=median(img)-100, max=median(img)+1000, /silent		;display subImage
					disp_med = median(subimg[(idisp = where((subimg ne -32768) and (finite(subimg) eq 1)))])
					disp_sig = robust_sigma(subimg[idisp])
					display, subimg, min=disp_med - 10.*disp_sig, max=disp_med + 50.*disp_sig, top=254, /aspect, /silent
		
					
	
					adxy, hdr, astrom.raj2000, astrom.dej2000, xcat, ycat
					xcat = xcat - minx + xoff
					ycat = ycat - miny + yoff
			
			
					;display, subImg, top=254, /aspect, min=median(img)-100, max=median(img)+1000, /silent		;display subImage
					;disp_med = median(subimg[(idisp = where((subimg ne -32768) and (finite(subimg) eq 1)))])
					;disp_sig = robust_sigma(subimg[idisp])
					;display, subimg, min=disp_med - 10.*disp_sig, max=disp_med + 50.*disp_sig, top=254, /aspect, /silent
		
					oplot, xCat, yCat, psym=6, color=fsc_color('green')						;overplot catalog objects
	
					cursor, xCCat, yCCat, /up													;select catalog object choice
	
					ind = where( ((xCat - xCCat)^2d + (yCat - yCCat)^2d) eq min((xCat - xCCat)^2d + (yCat - yCCat)^2d))	;closest catalog object to cursor location
	
					cursor, xPix, yPix, /up								;select location of corresponding star

					gcntrd, img, xPix, yPix, xCen, yCen, 5			;centroid star
	
					if (xcen eq -1) then xcen = xPix				;if centroid fails, adopt initial guess
					if (ycen eq -1) then ycen = yPix				
	
					xoff = xoff + (xcen - xcat[ind[0]])
					yoff = yoff + (ycen - ycat[ind[0]])
	
					xPix = xCen + minX								;convert back to mosaic coordinates
					yPix = yCen + minY
			
					minX = minX + xInc								;increment minimum and maximum x-pixels for next subImage
					maxX = maxX + xInc
		
					maxX = min([maxX, xSize])						;restrict to actual width of mosaic
			
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
		
				maxY = min([maxY, ySize])									;restrict to actual height of mosaic
	
			endfor
			undefine, subImg
			undefine, img
	
	
			iniFile = 'iniWCS.reg'									;name of initial WCS file containing x,y,ra,dec
	
			file_delete, iniFile, /allow_nonexistent
	
	
			openw, lun, iniFile, /get_lun							;save objects to file
				printf, lun, [transpose(list.x), transpose(list.y), transpose(list.ra), transpose(list.dec)]
			free_lun, lun
		
			str = strsplit(image, '.fits', /extract, /regex)
			newImage = pre + str[0] + '.fits'
			if (n_elements(str) gt 1) then newImage = newImage + str[-1]
		
			if (n_elements(overwrite) eq 0) then begin
				if (image ne NewImage) then begin
					file_delete, newImage, /allow_nonexistent
					spawn, 'cp ' + image + ' ' + newImage
				endif
			endif

			;jaa_astrom, list.x, list.y, list.ra, list.dec, newImage, error, nsig=2.5, maxiter=10	;[[list.era],[list.edec]]/1000d/3600d		
			jaa_astrom1, list.x, list.y, list.ra, list.dec, newImage, error, nsig=2.5, maxiter=10	;[[list.era],[list.edec]]/1000d/3600d
	
			hdrNew = headfits(newImage)
	
			;astrom = astro[use]
	
		endif else begin
			;newImage = pre + (strsplit(image, '.fits', /extract, /regex))[0] + '.fits'		
			str = strsplit(image, '.fits', /extract, /regex)
			newImage = pre + str[0] + '.fits'
			if (n_elements(str) gt 1) then newImage = newImage + str[-1]
			if (n_elements(overwrite) eq 0) then begin
				if (image ne NewImage) then begin
					file_delete, newImage, /allow_nonexistent
					spawn, 'cp ' + image + ' ' + newImage
				endif
			endif
			hdrNew = hdr
		endelse
	
		if (n_elements(rough) ne 0) then return
	
		catflag = 1
		while catflag do begin
			case strupcase(catalog) of
				'SDSS' 	:	BEGIN
								print, 'Querying SDSS...'
								sdss = queryvizier('II/294', [ra, dec], radius, /AllColumns)
								if (type(sdss,/string) eq 'longint') then begin
									catalog = 'USNO'
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
								astrom = query2mass(ra, dec, radius)
								if (type(astrom,/string) eq -1) then begin
									catalog = 'USNO'
									break
								endif else catflag = 0
							END
			endcase
		endwhile
	endif else begin
		hdrnew = hdr
		newimage = image
	endelse
	
		
	
	adxy, hdrnew, astrom.raj2000, astrom.dej2000, xnew, ynew

	list = arr_struct({x:xnew, y:ynew, ra:astrom.raj2000, dec:astrom.dej2000});, era:astrom.e_raj2000, edec:astrom.e_dej2000})
	list = list[where( (xnew ge 0) and (xnew le xSize) and (ynew ge 0) and (ynew le ySize))]

	openw, lun, 'list.reg', /get_lun
		for i=0,(n_elements(list)-1) do printf, lun, list[i].x, list[i].y, format='(2(F-20.9))'
	free_lun, lun
	
	

	print
	print, 'Centroiding ' + strtrim(n_elements(list),1) + ' objects.'


	pr = phot_parallel(arr_struct({x:list.x, y:list.y}), {IMAGE:newImage, FWHMPSF:5, APERTURES:5, $
			WCSIN:'physical', ZMAG:0, ITIME:100, DATAMIN:-100., DATAMAX:datamax, MAXSHIFT:20, CALGORITHM:'centroid', $
			cbox:20, scale:1., SALGORITHM:'mode', annulus:15, dannulu:10}, $
			'g', nthreads=nthreads*abs(n_elements(rough)-1), bridge=bridge, nodestroy=nodestroy)


	ind = where( (pr.cier eq 0) and (pr.sier eq 0) and (pr.pier eq 0) and (pr.merr ne 'INDEF'), count, complement=bind)

	if (count eq 0) then stop

	openw, lun, '0_baddata.reg', /get_lun
		for i=0,(n_elements(bind)-1) do printf, lun, pr[bind[i]].xcenter, pr[bind[i]].ycenter, format='(2(F-20.9))'
	free_lun, lun



	
	data = struct_addtags(struct_addtags(list, pr), arr_struct({raj2000:list.ra, dej2000:list.dec}))
	data.x = data.xcenter
	data.y = data.ycenter
	data = data[ind]

	openw, lun, 'data.reg', /get_lun
		for i=0,(n_elements(data)-1) do printf, lun, data[i].x, data[i].y, format='(2(F-20.9))'
	free_lun, lun

	openw, lun, 'astrom.reg', /get_lun
		for i=0,(n_elements(astrom)-1) do printf, lun, astrom[i].raj2000, astrom[i].dej2000, format='(2(F-20.9))'
	free_lun, lun



	;jaa_astrom, data.x, data.y, data.ra, data.dec, newImage, error, nsig=nsig, maxiter=10, result=astropar, /quiet	;[[data.era],[data.edec]]/1000d/3600d
	jaa_astrom1, data.x, data.y, data.ra, data.dec, newImage, error, nsig=nsig, maxiter=10, result=astropar, /quiet	;[[data.era],[data.edec]]/1000d/3600d

	
	use = astropar.use
	hdr = headfits(newImage)
	adxy, hdr, data[use].ra, data[use].dec, newx, newy


	if 0 then begin
		erase
		cleanplot, /silent
		psopen, 'diag_astrometry_0'+catalog, /encapsulated, /color, /heiles, xs=8, ys=6, /inches
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
	endif
	
	
	if (n_elements(object) gt 0) then begin
		img = mrdfits(newimage, 0, hdr, /silent)					;read in image plus header
		;hdr = headfits(newImage)
		file_delete, newImage+'.old', /allow_nonexistent
		dummy = file_search(newImage,count=nMatch)
		if (nMatch gt 0) then spawn, 'mv ' + newImage + ' ' + newImage + '.old'
		mwrfits, img, newImage, hdr					;write the image file and updated header
		if (n_elements(hdus) gt 0) then for i=0,(n_elements(hdus) - 1) do  mwrfits, hdus.(i), newImage	;write each of the hdus

		for tag=0,(n_tags(astrom)-1) do begin
			tagType = strupcase(type(astrom[0].(tag),/string))
			if (tagType ne 'STRING') then begin
				iBad = where(~finite(astrom.(tag)),nBad)
				if (nBad gt 0) then astrom[iBad].(tag) = -1
			endif
		endfor

		mwrfits, astrom, newimage								;write the unculled object list
		hdr = headfits(newImage)
		sxaddpar, hdr, 'NEXTENS', n_elements(hdus)+1, 'Number of extensions.', /savecomment
		modFits, newImage, 0, hdr
		file_delete, newImage+'.old', /allow_nonexistent
	endif
	centroid = data						;all the centroided points

end

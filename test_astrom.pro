FUNCTION read_files, files, filenames
	filenames = file_search(files)
	return, filenames
END

PRO supCamData::astrom_store, centroid, shift, cats, shifts
	COMPILE_OPT, idl2, HIDDEN

	for i=0l, (n_elements(shfit) - 1l) do begib
		ind = where((centroid.x gt shift[i].delx) and 				$ 
			    (centroid.x lt (shift[i].delx + shift[i].xsize)) and	$
			    (centroid.ygt shift[i].dely) and 				$
			    (centroid.y lt (shift[i].dely + shift[i].ysize)), nobj)	$
		
		if (nobj eq 0) then stop

		dum = create_struct('num'+strtrim(i,2),					$	 
		      arr_struct({file:replicate(strtrim(shift[i].file,2),nobj), 	$
		      x:centroid[ind].x - shift[i].delx, y:centroid[ind].y - 		$
		      shift[i].dely, ra:centroid[ind].ra, dec:centroid[ind].dec}))	$
		      
		cat = struct_addtags(cat, dum)			    
	endfor

	if (n_elements(cats) gt 0) then ntags = strtrim(n_tags(cats),2) else ntags = '0'

	dcats  = create_struct('num'+ntags, cat)
        cats   = struct_addtags(cats, dcats)
        dshift = create_struct('num'+ntags, shift)
        shifts = struct_addtags(shifts, dshift)
END

FUNCTION getParams, data, date, filenames=filenames, start=star
	if (n_elements(start) eq 0) then begin
		CDELT1 = -5.6093104d-05     ;initial guess at RA scale (deg/pixel)
                CDELT2 = 5.6095829d-05

		for i=0,(n_tags(data)-1) do dhold = struct_append(dhold, data.(i))
		
		CRVAL1 = median(dhold.RA)
	        CRVAL2 = median(dhold.DEC)
        	undefine, dhold
		
		if (date lt 2008.58) then begin
                        if (strmid(filenames[0],0,1) eq 'n') then begin
                                CRPIX1 = [5372.53, 3285.94, 1162.94, -956.959, -3044.84,
					  5354.53, 3269.39, 1156.98, -966.578, -3053.81]
                                CRPIX2 = [-49.7044, -77.3883, -83.3438, -76.6404, -91.2516,
					  4005.36, 4019.26, 4018.71, 4004.06, 3979.33]
                                CROTA2 = -1*[-6.0914919e-03, -2.0648508, -4.9069205, 4.2234785,
					     1.3125425, 1.5521885, 4.2652558, 4.6701419, 
					    -1.6867016e-04, 0.71141694]
                        endif else begin
                                CRPIX1 = [5246.5142, 3165.0548, 1054.4712, -1066.5537, -3161.3642,
					  5244.9220, 3163.8577, 1051.1501, -1065.2115, -3159.0162]
                                CRPIX2 = [35.490705, 31.809698, 33.976075, 32.933892, 30.762235,
					  4105.4282, 4120.2479, 4118.9465, 4118.8147, 4102.6354]
                                CROTA2 = -1*[4.5489500, 4.8252375, 4.3339790, 4.5679359, 4.3634360,
    					     4.0064291, 4.0187368, 3.9573857, 4.4091170, 4.6163363]*1d-3
                                ; This method doesn't work if there is no rotation.
                                getrot, headfits(filenames[0]), imrot, cdelt, DEBUG=1
                                ANGLE = imRot*!dtor
                        endelse
                endif else begin
                        if (strmid(filenames[0],0,1) eq 'n') then begin
                        endif else begin
                        endelse
                        CRPIX1 = [5425.8691, 3351.0452, 1241.0316, -871.95972, -2964.5720,
				  5423.8774, 3351.3420, 1242.9849, -877.50110, -2966.1489]
                        CRPIX2 = [111.32944, 111.39816, 111.15529, 113.26266, 113.70063,
				  4315.5620, 4330.3096, 4338.3032, 4330.0840, 4318.5273]
                        CROTA2 = [-0.0055500342, -0.0053783910, -0.0054863331, -0.0053572311,
			  	  -0.0051588219, -0.0050752364, -0.0053784340, -0.0049788447,
				  -0.0052406706, -0.0055207304]
                        getrot, headfits(filenames[0]), imrot, cdelt, DEBUG=1
                        ANGLE = imRot*!dtor
                        CROTA2 = CROTA2 * 0d
                endelse
        endif else begin
                CDELT1 = start.CDELT1
                CDELT2 = start.CDELT2
                CRVAL1 = start.CRVAL1
                CRVAL2 = start.CRVAL2
                CRPIX1 = start.CRPIX1
                CRPIX2 = start.CRPIX2
                CROTA2 = start.CROTA2
                ANGLE = start.ANGLE
        endelse

        param = {CDELT1:CDELT1, CDELT2:CDELT2, CRVAL1:CRVAL1, CRVAL2:CRVAL2, 
		 CRPIX1:CRPIX1, CRPIX2:CRPIX2, CROTA2:CROTA2, ANGLE:ANGLE}

       ; Why is this defined twice?
        angleRange = !pi/(sqrt(findgen(10)+1))^3d
        ; I don't really understand what's going with these two lines...
        angleRange = dblarr(3) + 1d-1   ;angleRange[0]

	a = 0
	while (a lt n_elements(angleRange)) do begin
        	for i=0,10 do begin
	                param = [CDELT1, CDELT2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CROTA2, ANGLE]

	                ; What is the point of this line? What is res?
        	        z = model(data, param, res=res)
                	if (i eq 0) then res0 = res
	                ;window, 2
        	        ;plot, res[*,0], res[*,2], psym=1, title=roundx(likelihood(data, param),2)
                	;window, 0
                	;plot, res[*,0], res[*,3], psym=1, title=roundx(likelihood(data, param),2)
                	;oplot, res0[*,0], res0[*,3], psym=1, color=fsc_color('red')

	                ;want to check if the declination offset can be fixed with rotation
        	        CRVAL1 -= median(res[*,2])*cos(CRVAL2)/3600d
                	CRVAL2 -= median(res[*,3])/3600d
        	endfor
        	nAngles = 500
        	angles = (findgen(nAngles)/(nAngles - 1d) - 0.5)*2d*angleRange[a] + ANGLE
        	like = dblarr(nAngles)

        	for k=0,(nAngles - 1) do begin
                	ANGLE = ANGLE*0d + angles[k]
                	param = [CDELT1, CDELT2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CROTA2, ANGLE]
                	;z = model(data, param, res=res)
                	like[k] = likelihood(data, param)
        	endfor
		
		;pause
		;window, 1
	        ;plot, angles, like
       		;wset, 0

        	dummy = max(like,imax)
		; I don't really get why the index of the max would matter.
		; Is it just a test on the shape of the likelihood distribution?
        	if ((imax lt 200) or (imax gt 300)) then stop
        	;print, dummy
        	; Get the "final" Angle?
        	ANGLE = (angles[imax])[0]
        	param = [CDELT1, CDELT2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CROTA2, ANGLE]

	++a
	endwhile
	
	return, param
END

PRO test_astrom
	filenames = read_files('sA*')
	for i=1, (n_elements(mosfile) - 1) do begin
		if (i eq (n_elements(mosfile) - 1)) then 
			if (n_elements(nodestroy) ne 0) then 
				undefine, nodestroy
		
		if (n_elements(nobatch) eq 0) then begin
			hdr1 = headfits(mosfile[i])    
                	shdr1 = headfits(strtrim((mrdfits(mosfile[i], 1, /silent))[0].file))
                	xyxy, shdr1, shdr0, 0, 0, xdif, ydif    
                	extast, hdr0, astr              
                	astr.crval = astr.crval - [hms2dec(sxpar(shdr0, 'RA2000'))*15d, 	$
				     hms2dec(sxpar(shdr0, 'DEC2000'))] + 			$
                             	     [hms2dec(sxpar(shdr1, 'RA2000'))*15d, 			$
			    	     hms2dec(sxpar(shdr1, 'DEC2000'))]      			 
                 	putast, hdr1, astr      ; put this astrometry into the mosaics header
                 	modfits, mosfile[i], 0, hdr1    ; put the header back in the file
		endif
		mwrfits, astrom, mosfile[i], /silent 
            	 ; Run astrometry in batch mode now that the astrometric solution is close to correct
		 astrometry, mosfile[i], /object, /batch, /overwrite, catalog=catalog,		$ 
			     datamax=datamax[i], nthreads=nthreads, bridge=bridge, 		$
			     nodestroy=nodestroy, centroid=centroid 				$
		shifts = mrdfits(mosfile[i], 1, /silent)
		self->astrom_store, centroid, shift, cats, shifts
	endfor
	params = getParams(data, date, filenames)
END

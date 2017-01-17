



function supalign2, files, datamax, suffix=suffix, shift=shift, force=force
	
	numChips = 10
	numFiles = n_elements(files)
	numExps = numFiles/numChips
	
	date = fix(strtrim(strsplit(sxpar(headfits(files[0]), 'DATE-OBS'),'-',/extract)))								;extract the date from the domeflat
	date = date[0] + (date[1]-1)/12. + date[2]/365.															;turn this into a decimal
	;if ((date lt 20009) and (date gt 0)) then begin
	if (date gt 2008.58) then begin
		delxchip = round([0.0, 2076.6, 4187.7, 6303., 8397.3, 2.77, 2076.98, 4186.27, 6306.60, 8395.98])			;initial x-offset between chips which should be roughly correct
		delychip = round([0.0, 0.71, 1.85, 0.017, -1.40, -4200.88, -4215.95, -4224.48, -4218.27, -4204.58])		;initial y-offset
	endif else begin
		if (strmid(files[0],0,1) eq 'n') then begin
			delxchip = round([0, 2087, 4207, 6327, 8414, 19, 2100, 4203, 6341, 8425])
			delychip = round([0, 20, 36, 31, 48, -4060, -4067, -4070, -4052, -4030])
		delxchip = round([0, 2087, 4207, 6327, 8414, 19, 2100, 4203, 6341, 8425])
		delychip = round([0, 20, 36, 31, 48, -4096, -4076, -4060, -4065, -4048])
		endif else begin
			delxchip = round([	35., 	2115., 	4225., 	6350., 	8440., 50, 	2135, 	4250., 	6365., 	8450.])			;initial x-offset between chips which should be roughly correct
			delychip = round([	-75, 	-65., 	-55, 	-50., 	-35., -4140., -4150., -4140., -4130., -4105.])		;initial y-offset
		endelse
	endelse
	
	;endif else begin
	;	delxchip = round([0.0, 2077, 4185, 6302, 8391, 5, 2078, 4185, 6306, 8395])
	;	delychip = round([0.0, 1, 3, -3, -2, -4204, -4220, -4226, -4219, -4208])
	;endelse
	
	
	;if (date gt 2008.58) 	then chip = ['chihiro', 'clarisse', 'fio', 'kiki', 'nausicaa', 'ponyo', 'san', 'satsuki', 'sheeta', 'sophie'] $
	;						else chip = ['w67c1', 'w6c1', 'si006s'	, 'si002s'	, 'w7c3', 'w93c2'	, 'w9c2', 'si005s'	, 'si001s'	, 'w4c5']
										;[8, 5, 2, 1, 6, 9, 7, 3, 0, 4]
	
	
	outfiles = strarr(numExps)
	if (date le 2008.58) then order = [5,6,3,1,7,8,9,2,0,4] else order = indgen(10)	;[8, 5, 2, 1, 6, 9, 7, 3, 0, 4]	;indgen(10)
	file = ''
	for i=0,(numExps - 1) do file = [file, (files[i*10:i*10+9])[order]]
	file = file[1:*]
		
	for i=0,(numExps - 1) do begin
		outfile = 'tmos_' + (strsplit(file[i*10], '_', /extract, /regex))[1]
		if (n_elements(suffix) gt 0) then outfile = outfile + suffix
		outfiles[i] = outfile + '.fits'
	endfor
	
	xsize = fltarr(numChips)
	ysize = fltarr(numChips)
	gain = fltarr(numchips)
	;scale = fltarr(numchips)
	datamax = dblarr(numExps) + 1d30
	
	msg = 'Building mosaic'
	if (numExps gt 1) then msg = msg + 's'
	for i=0,(numExps - 1) do begin
		
		dummy = file_search(outfiles[i]+'.gz', count=exist)	;check to see if the mosaic file already exists
		if (exist gt 0) then begin				;if so then 
			print, "MOSAIC ALREADY EXISTS
			; Uncommenting this out to try to skip this step if I've run it already. - Alexa
			fits_info, outFiles[i], n_ext=n_ext, /silent	;count the number of extensions
			n_ext = sxpar(headfits(outFiles[i]+'.gz'), 'NEXTENS')
			if (n_ext eq 0) then exist = 0 $		;if there are none, then the mosaic needs to be remade
			else begin					;otherwise
				shift = mrdfits(outfiles[i]+'.gz',1,/silent)	;get the shift structure
				if (type(shift,/string) eq 'structure') then begin
					dummy = where(file[i*10] eq strtrim(shift.file), nmatch) ;see if the files making up the mosaic match the first file
					if (nmatch eq 0) then exist = 0				;if not, then the mosaic needs to be remade
				endif else nmatch = 0
			endelse
		endif
	
	
		for c=0,(numChips - 1) do begin			;cycle through each chip and ...
				xsize[c] = sxpar((hdr = headfits(file[i*10+c])), 'NAXIS1')
				ysize[c] = sxpar(hdr, 'NAXIS2')		
				if (c eq 0) then hdr0 = hdr
				;gain[c] = sxpar(hdr, 'GAIN')
				;scale[c] = sxpar(hdr, 'CHPSCALE')
				;if (scale[c] eq 0) then scale[c] = (mrdfits('supcamchipscale.fits',1,/silent))[c].scale
				datamax[i] = min([datamax[i], sxpar(hdr,'SATURATE')])
				if (datamax[i] eq 0) then stop
				if (exist and (n_elements(force) eq 0)) then statusbar, i*numChips + c, numchips*numExps, msg=msg
		endfor
		shift = arr_struct({file:file[i*10:(i+1)*10-1], xsize:xsize, ysize:ysize, delx:delxchip, dely:delychip + abs(min(delychip))})
		
		if (~exist or (n_elements(force) eq 1)) then begin										;if the file doesn't exist or the force keyword is set, then remake mosaic
			mos = fltarr(max(delxchip + xsize), max(abs(delychip)) + max(delychip[0:4] + ysize[0:4])) - 32768		
	
			for c=0,(numChips - 1) do begin
				mos[delxchip[c], delychip[c] + abs(min(delychip))] = mrdfits(file[i*10+c], 0, /silent)
				statusbar, i*numChips + c, numchips*numExps, msg=msg
			endfor
			
			bPix = where(~finite(mos),nBPix)
			if (nBPix gt 0) then mos[bPix] = -32768
	
			;file_delete, outfiles[i] + '.gz', /allow_nonexistent
		
			mwrfits, mos, outfiles[i], /create
	
			addhead, outFiles[i], file[i*numChips:(i + 1)*numChips - 1], hdr=hdr0, xoff=widbuf, yoff=abs(min(delychip))		;create a header
			hdr = headfits(outFiles[i])
			;sxaddpar, hdr, 'SATURATE', 50000.*gain[c]*scale[c]/exptime, /saveComment
	
			mwrfits, shift, outFiles[i]
			
			sxaddpar, hdr, 'NEXTENS', 1, 'Number of extensions.', /savecomment
			modfits, outFiles[i], 0, hdr
			;spawn, 'gzip ' + outFiles[i]
		endif else begin
			;fits_info, outFiles[i], n_ext=n_ext, /silent
			n_ext = sxpar(hdr, 'NEXTENS')
			if (n_ext eq 0) then begin
				mwrfits, shift, outFiles[i]
				hdr = headfits(outFiles[i])
				sxaddpar, hdr, 'NEXTENS', 1, 'Number of extensions.', /savecomment
				modFits, outFiles[i], 0, hdr
				;file_delete, outfiles[i] + '.gz', /allow_nonexistent
				;spawn, 'gzip ' + outFiles[i]
			endif
		endelse
		
		if (i eq 0) then mosfile = outfiles[i] else mosfile = [mosfile, outfiles[i]]
		;mosfile = mosfile + '.gz'
		
	endfor
	
	return, mosfile

end




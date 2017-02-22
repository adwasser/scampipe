;+
;
;-
;-------------------------------------------------------------------------------------------

FUNCTION supCamData::init, filter=filter, object=object, flatType=flatType, outFile=outFile, $
		pre=pre, batch=batch, coaddType=coaddType
	compile_opt idl2, hidden
	self->getChipsDate
	self.badPixValue = !values.f_nan
	if (n_elements(filter) gt 0) then self.filter = filter else self.filter = 'NULL'
	if (n_elements(object) gt 0) then self.object = object else self.object = 'NULL'
	if (n_elements(flatType) gt 0) then self.flatType = flatType else self.flatType = 'NULL'
	if (n_elements(outFile) gt 0) then self.outFile = outFile else self.outFile = 'NULL'
	if (n_elements(pre) gt 0) then self.pre = pre else self.pre = 'NULL'
	if (n_elements(batch) gt 0) then self.batch = 1 else self.batch = 0
	if (self.date gt 2008.58) then self.satLevel = 2^16 else self.satLevel = 2^15			;set the saturation level depending on the date
	if (n_elements(coaddType) eq 0) then self.coaddType = 'median' else self.coaddType = coaddType
	self.nonLinFrac = 0.9																	;assume non-linearity sets in around 90% of full well depth
	self.satMaskGrow = [2,2]
	self.crGrow = 1
	return, 1
END

;-------------------------------------------------------------------------------------------

PRO supCamData::cleanup
	compile_opt idl2, hidden
	;-- free memory allocated to pointer when destroying object
END

;-------------------------------------------------------------------------------------------

FUNCTION supCamData::hdr_search, commands
;--Search file headers for a given keyword

	list = ['SUPA','H','To_R','f','z','n','q','g','A','c','as','x']
	listing = ''
	for i=0,(n_elements(list)-1) do listing = [listing,list[i]+'*.fits',list[i]+'*.fits.gz','dir'+list[i]+'/'+list[i]+'*.fits','dir'+list[i]+'/'+list[i]+'*.fits.gz']
 	listing = listing[1:*]

	flag = 0
	foreach filetype, listing do begin
		file = file_search(filetype, count=nFiles)
		for i=0,(nfiles - 1) do begin
			for c=0,(n_elements(commands)-1) do dummy = execute(commands[c])
			if flag then break
		endfor
		if flag then break
	endforeach
	if (flag eq 0) then value = 'NULL'
	return, value
END

;-------------------------------------------------------------------------------------------

PRO supCamData::checkInput, ndFilt=ndFilt, ndFlat=ndFlat, ndObj=ndObj, ndOut=ndOut
;--Check for necessary keywords
	COMPILE_OPT idl2, HIDDEN

	if ((self.filter eq 'NULL') and (n_elements(ndFilt) ne 0)) then begin
		self.filter = self->hdr_search(["object = strtrim(sxpar((hdr = headfits(file[i])),'OBJECT'),2)", $
									"if (object ne 'BIAS') and (strpos(object,'FLAT') eq -1) and (object ne '0') then flag = 1", $
									"if flag then value = strtrim(sxpar(hdr,'FILTER01'),2)"])
		if (self.filter eq 'NULL') then begin
			print
			self.filter = choose('Which filter? (Must be exact header syntax)', 'W-S-G+')
		endif
	endif

	if ((self.object eq 'NULL') and (n_elements(ndObj) ne 0)) then begin
		self.object = self->hdr_search(["value = strtrim(sxpar(headfits(file[i]),'OBJECT'),2)", $
									"if (value ne 'BIAS') and (strpos(value,'FLAT') eq -1) and (value ne '0') then flag = 1"])
		if (self.object eq 'NULL') then begin
			print
			self.object = choose('Name of science object: (Must be exact header syntax)', 'NGC 3377')
		endif
	endif

	if (n_elements(ndFlat) ne 0) then begin
		if (self.flatType ne 'NULL') then begin
			if ((self.flattype ne 'DOMEFLAT') and (self.flattype ne 'SKYFLAT') and (self.flattype ne 'SUPERSKYFLAT')) then begin
				print
				self.flatType = choose('What type of flats? (Must be exact header syntax)', 'DOMEFLAT')
			endif
		endif else begin
			print
			self.flatType = choose('What type of flats? (Must be exact header syntax)', 'DOMEFLAT')
		endelse
	endif

	if ((self.outFile eq 'NULL') and (n_elements(ndOut) ne 0)) then	self.outfile = choose('Name of final mosaic:', 'mos_x.fits')
END

;-------------------------------------------------------------------------------------------

PRO supCamData::write2file, list, outFile
;--Write list to file
	COMPILE_OPT idl2, HIDDEN
	;spawn, 'ls ' + outFile, result																;check for outfile
	;if (result) then spawn, 'rm ' + outFile													; remove if found
	;result = file_search(outfile, count=nres)
	file_delete, outfile, /allow_nonexistent
	openw, lun, outFile, /get_lun																;save list of flats to file
		printf, lun, transpose(list)
	free_lun, lun
END

;-------------------------------------------------------------------------------------------

PRO supCamData::mkDir, str
	dir = (file_search(str, count=ndir))[0]
	if (ndir lt 1) then FILE_MKDIR, str
END

;-------------------------------------------------------------------------------------------

FUNCTION supCamData::file_check, count
;--Create a directory or retrieve files from it given a prefix
	self->mkDir, 'dir'+self.pre
	files = file_search('dir'+self.pre+'/*fits*', count=nFiles)
	if (nFiles gt 0) then file_move, files, './'
	return, file_search(self.pre+'*fits*', count=count)											;search for the files in the current directory
END

;-------------------------------------------------------------------------------------------

FUNCTION supCamData::getFiles, count, pre=pre, dirPre=dirPre
;--Create a directory or retrieve files from it given a prefix
	if (n_elements(dirPre) eq 0) then dirPre = 'dir'
	if (n_elements(pre) eq 0) then pre = self.pre
	dir = (file_search(dirPre+pre, count=ndir))[0]
	if (nDir gt 0) then begin
		files = file_search(dirPre+pre+'/*fits*', count=nFiles)
		if (nFiles gt 0) then file_move, files, './'
	endif
	return, file_search(pre+'*fits*', count=count)												;search for the files in the current directory
END

;-------------------------------------------------------------------------------------------

FUNCTION supCamData::getFiles1, dirPre=dirPre, pre=pre, count=count
;--Create a directory or retrieve files from it given a prefix
	if (n_elements(pre) eq 0) then pre = ''
	if (n_elements(dirPre) gt 0) then begin
		dir = (file_search(dirPre, count=ndir))[0]
		if (nDir gt 0) then begin
			files = file_search(dirPre+'/*fits*', count=nFiles)
			if (nFiles gt 0) then file_move, files, './'
		endif
	endif
	return, file_search(pre+'*fits*', count=count)												;search for the files in the current directory
END

;-------------------------------------------------------------------------------------------

PRO supCamData::gZip, str, msg=msg0
	if (n_elements(msg0) gt 0) then msg = msg0
	files = file_search(str, count=nFiles)
	for i=0,(nFiles - 1) do begin
		spawn, 'gzip '+files[i]
		if (nFiles gt 1) then statusBar, i, nFiles, msg=msg
	endfor
END

;-------------------------------------------------------------------------------------------

PRO supCamData::putInDir, str, dir
	files = file_search(str, count=nFiles)
	if (nFiles gt 0) then begin
		self->mkDir, dir
		spawn, 'mv '+str+' '+dir+'/'
	endif
END

;-------------------------------------------------------------------------------------------

PRO supCamData::fileMove, oldName, newName
	files = file_search(oldName,count=nFiles)													;
	if (nFiles gt 0) then file_move, oldName, newName
END

;-------------------------------------------------------------------------------------------

PRO supCamData::fileDelete, str
	files = file_search(str,count=nFiles)														;
	if (nFiles gt 0) then file_delete, files
END

;-------------------------------------------------------------------------------------------

FUNCTION supCamData::rmGZ, fileNames
	newNames = fileNames
	for i=0,(n_elements(fileNames) - 1) do begin
		ipos = strpos(fileNames[i],'.gz')
		if (iPos gt -1) then newNames[i] = strtrim(strmid(fileNames[i],0,iPos),2)
	endfor
	return, newNames
END

;-------------------------------------------------------------------------------------------

PRO supCamData::gunzip, inFileNames, outFileNames
	nFiles = n_elements(inFileNames)
	outFileNames = self->rmGZ(inFileNames)
	for i=0,(nFiles - 1) do begin
		if (inFileNames[i] ne outfileNames[i]) then begin
			file = file_Search(inFileNames[i],count=nMatch)
			if nMatch then spawn, 'gunzip -f '+inFileNames[i]
		endif
	endfor
END

;-------------------------------------------------------------------------------------------

PRO supCamData::getChipsDate
;--Return appropriate list of chips given a hdr containing the date

	str = 'Determining date and ccd names...'
	print & print, str, format='(A-'+strtrim(strlen(str),2)+',$)'
	value = self->hdr_search(["object = strtrim(sxpar((hdr = headfits(file[i])),'OBJECT'),2)", $
								"if (object ne 'BIAS') and (strpos(object,'FLAT') eq -1) and (object ne '0') then flag = 1", $
								"if flag then date = fix(strtrim(strsplit(sxpar(hdr, 'DATE-OBS'),'-',/extract)))", $
								"if flag then value = strtrim(date[0] + (date[1]-1)/12. + date[2]/365.,2)", $
								"if flag then if (value gt 2008.58) then value = [value, 'chihiro', 'clarisse', 'fio', 'kiki', 'nausicaa', 'ponyo', 'san', 'satsuki', 'sheeta', 'sophie']" + $
								" else value = [value, 'w67c1', 'w6c1', 'si006s', 'si002s', 'w7c3', 'w93c2', 'w9c2', 'si005s', 'si001s', 'w4c5']"])
	if (value[0] eq 'NULL') then stop
	date = float(value[0])
	chips = value[1:*]
	if 0 then begin
		;old way
		foreach filetype, listing do begin
			files = file_search(filetype, count=nFiles)
			if (nFiles gt 0) then begin
				foreach tag, ['DOMEFLAT','SKYFLAT','SUPERSKYFLAT','BIAS'] do begin
					files = self->supCheckConstraints(Files, !null, {OBJECT:tag}, count=count)
					if (count eq 0) then break
				endforeach
				if (count gt 0) then begin
					hdr = headfits(files[0])
					date = fix(strtrim(strsplit(sxpar(hdr, 'DATE-OBS'),'-',/extract)))			;extract the date from the domeflat
					date = date[0] + (date[1]-1)/12. + date[2]/365.

					if (date gt 2008.58) then $
						chips = ['chihiro', 'clarisse', 'fio', 'kiki', 'nausicaa', 'ponyo', 'san', 'satsuki', 'sheeta', 'sophie'] $
						else chips = ['w67c1', 'w6c1', 'si006s', 'si002s', 'w7c3', 'w93c2', 'w9c2', 'si005s', 'si001s', 'w4c5']
					break
				endif
			endif
		endforeach
	endif
	print, 'done!'
	print
	print, date
	self.chips = chips
	self.date = date
END

;-------------------------------------------------------------------------------------------

PRO supCamData::rename_files
;--Rename Files (change everything from SUPA* to H*)
	COMPILE_OPT idl2, HIDDEN
	self.pre = 'SUPA'

	self->checkInput, /ndfilt
	rfiles = self->getFiles(nFiles)																;retrieve files
	print
	if (nFiles gt 0) then begin
		print, 'Renaming Files...'

		files = [self->supCheckConstraints(rFiles, {FILTER01:self.filter}, {OBJECT:'BIAS'}), self->supCheckConstraints(rFiles, {OBJECT:'BIAS'})]
		if (n_elements(files) eq 0) then begin
			message, 'Their are no files with this filter.', /continue, /informational
			return
		endif


		self->write2file, files, 'namechange.lis'


; Commented out because of problems with sdfred paths.
; Run namechange for within the shell and restart program
		spawn, 'namechange.csh namechange.lis'													;change the names



		self->putInDir, self.pre+'*.fits', 'dir'+self.pre

		self.pre = 'H'
		self->mkDir, 'dir'+self.pre
		;add the keywords here



		if (self.date le 2008.58) then begin
			files = file_search(self.pre+'*.fits', count=count)									;identify all the renamed files
			for i=0,(count - 1) do begin														;cycle through each one
				hdr = headfits(files[i])														;	read in its header
				for c=0,3 do begin																;	add in these keywords to make the pipeline behave in the gaincalc portion
					SXADDPAR, hdr, 'S_EFMN' + strtrim(c+1,2) + '1', c*512 + 1, /SaveComment
					SXADDPAR, hdr, 'S_EFMN' + strtrim(c+1,2) + '2', 1, /SaveComment
					SXADDPAR, hdr, 'S_EFMX' + strtrim(c+1,2) + '1', (c+1)*512, /SaveComment
					SXADDPAR, hdr, 'S_EFMX' + strtrim(c+1,2) + '2', 4096, /SaveComment
				endfor
				modfits, files[i], 0, hdr														;store the header back to the file
			endfor
		endif



		files = file_search(self.pre+'*.fits', count=nFiles)
		for i=0,(nFiles - 1) do begin
			self->gZip, files[i]																;gZip the new files
			self->putInDir, files[i]+'.gz', 'dir'+self.pre									;move them to a directory
			statusBar, i, nFiles, msg='GZip-ing files ...'
		endfor
		self->fileDelete, 'tmp*'																;clean up an temporary files
		self->fileDelete, 'namechange.lis'														;ditto

		;print, 'Done!'
	endif else begin
		print
		print, 'No '+self.pre+'*.fits files to rename.'
		print
	endelse
	self.pre = 'H'
	self->putInDir, self.pre+'*.fits.gz', 'dir'+self.pre										;move them to a directory
	;print
END

;-------------------------------------------------------------------------------------------

PRO supCamData::get_gain, gain=gain
;--Calculate Gain and Multiply by it (change To_RH* to gTo_RH*)
	COMPILE_OPT idl2, HIDDEN
	if (self.pre eq 'NULL') then self.pre = 'f'

	self->checkInput, /ndfilt, /ndobj

	gainFiles = self->getFiles1(pre='supCamGain',dirPre='GAIN', count=nGainfiles)				;retrieve existing superbiases

	if (nGainFiles gt 0) then begin
		rGainFiles = file_search('supCamGain_'+self.filter+'.fits',count=nMatch)
		if nMatch then redo = choose('Gain file already exists, recalculate?', 'no') else redo = 'no'
	endif else redo = 'YES'
	if (strupcase(redo) eq 'NO') then begin
		gain = create_struct(self.chips[0], mrdfits(rgainfiles[-1],1,/silent))						;calculate gain
		for i=2,10 do gain = struct_addtags(gain, create_struct(self.chips[i-1], mrdfits(rgainfiles[-1],i,/silent)))
		;spawn, 'mv supcamgain_'+filter+'.fits GAIN/'
	endif else begin
		;print
		;print, 'Calculating Gain & Read Noise...'

		rFiles = self->getFiles(nFiles, pre='H')												;make a list of all the H* files
		if (nFiles eq 0) then stop
		biases 	= self->supCheckConstraints(rFiles, {OBJECT:'BIAS'}, count=nBiases)
		if (nBiases eq 0) then begin
			print, 'There are no biases!'
			stop
		endif
		biases = biases[0:min([50,n_elements(biases)])-1]

		rFiles = self->getFiles(nFiles, pre='To_R')												;make a list of all the To_R* files
		if (nFiles eq 0) then stop

		flats 	= self->supCheckConstraints(rFiles, {FILTER01:self.filter, OBJECT:'DOMEFLAT'},count=nFlats)	;make sure to use dome or skyflats for this
		if (nFlats eq 0) then begin
			print, 'There are no available flats for filter '+self.filter+'!'
			allFlats = self->supCheckConstraints(rFiles, {OBJECT:'DOMEFLAT'},count=nAllFlats)	;make sure to use dome or skyflats for this
			if (nAllFlats eq 0) then print, 'There are non available flats for any filter!' $
			else begin
				for i=0,(nAllFlats - 1) do begin
					hdr = headfits(allFlats[i])														;get the header
					print, allFlats[i], sxpar(hdr,'FILTER01')
				endfor
			endelse
			stop
		endif
		flats = flats[0:min([50,n_elements(flats)])-1]
		gain = self->gainCalc(flats, biases, outfile='supCamGain_'+self.filter+'.fits')
	endelse

	self->putInDir, 'supCamGain_'+self.filter+'.fits', 'GAIN'
	self->putInDir, 'H*.fits.gz', 'dirH'														;return the H* files to their directory
END

;-------------------------------------------------------------------------------------------

FUNCTION supCamData::gainCalc, flats, biases, outFile=outFile, monitor=monitor

	flat = flats
	bias = biases

	if (self.chips[0] eq 'chihiro') then pmax = 50000 else pmax = 25000

	exptimes = fltarr(n_elements(flat))
	for i=0,(n_elements(flat)-1) do exptimes[i] = float(sxpar(headfits(flat[i]),'EXPTIME'))		;get the exposure time for each flat

	if (difference(minmax(exptimes)) ne 0) then begin											;if the exposure times are not all the same
		uniqtimes = exptimes[uniq(exptimes, sort(exptimes))]									;identify the unique exposure times
		times = fltarr(n_elements(uniqtimes))
		for i=0,(n_elements(times)-1) do times = n_elements(where(exptimes eq uniqtimes[i]))	;count the number of exposures with each unique exposure time
		flat = flat[where(exptimes eq uniqtimes[(where(times eq max(times)))[0]])]				;keep only the exposures with the most common exposure time
	endif

	k = 0
	for c=0,(n_elements(self.chips)-1) do begin														;for each chip
		dumgain = create_struct(self.chips[c], {mean:fltarr(4), sigma:fltarr(4), ron:fltarr(4), ronsig:fltarr(4), mni:fltarr(4)})
		if (c eq 0) then gains = dumgain else gains = struct_addtags(gains, dumgain)			;create or append to a gains structure to hold the results

		detector = self.chips[c]																		;the chip under consideration
		;detector = strtrim(sxpar(headfits(flat[c]), 'DETECTOR'),2)
		listflat = self->supcheckconstraints(flat, {DETECTOR:detector}, count=nListF)			;keep only the flats for this chip
		listbias = self->supcheckconstraints(bias, {DETECTOR:detector}, count=nListB)			;keep only the biases for this chip

		if (n_elements(listflat)*n_elements(listbias) gt 0) then begin							;if there are files then continue
			for type=0,1 do begin																;for each type (flat and bias)
				case type of
					0:	list = listflat															;choose the appropriate list
					1:	list = listbias
				endcase
				nList = n_elements(list)

				med = fltarr(4, n_elements(list))												;create an array to hold the median values
				var = fltarr(4, n_elements(list), n_elements(list))								;create an array to hold the variance values

				for i=0,(nList-1) do begin														;for each file in this list
					subImg0 = getsubimgs(list[i])												;retrieve subchip images from input file
					for sc=0,3 do begin																	;for each subchip
						gpix = where( (subImg0.(sc) gt 0) and (subImg0.(sc) lt pmax), count)			;identify the good pixels
						if (count gt 0) then dmed = median(subImg0.(sc)[gpix]) else stop			;calculate a median value
						med[sc,i] = dmed																;store to the array of median values
					endfor

					if (i lt (nList - 1)) then begin
						for f=(i+1),(nList-1) do begin								;for each subsequent image on the same detector
							subImg1 = getsubimgs(list[f])										;retrive subchip images
							for sc=0,3 do begin													;for each subchip
								timg = (subImg0.(sc) - subImg1.(sc))[170:340,2000:2170]								;create a difference image
								gpix = where( (subImg0.(sc)[170:340,2000:2170] gt 0) and (subImg0.(sc)[170:340,2000:2170] lt pmax) and (subImg1.(sc)[170:340,2000:2170] gt 0) and (subImg1.(sc)[170:340,2000:2170] lt pmax), count)	;identify the good pixels
								if (count gt 0) then dvar = robust_sigma(timg[gpix])^2d else stop	;calculate the variance of the difference image
								var[sc,i,f] = dvar												;save to the array of variances

								if (n_elements(monitor) gt 0) then begin

									mid = biweight_mean(timg[gpix])
									sig0 = robust_sigma(timg[gpix])

									range = mid + [-5,5]*sig0
									plothist, timg[gpix], xrange=range, xhist, yhist, title=chip + ' | chiplet = ' + strtrim(sc+1,2) + ' type = ' + (['FLAT', 'BIAS'])[type]
									npixel = 100
									xpsf = findgen(nPixel)/float(nPixel)*difference(range) + min(range)
									sig = (sig0/difference(range)*nPixel)
									xdummy = min(abs(mid - xpsf), mind)
									ipsf = psf_gaussian(npixel=nPixel, st_dev=sig, centroid=mind, /double, ndimen=1)*max(yhist)
									oplot, xpsf, ipsf, color=fsc_color('red')

									gpix = where( (subImg0.(sc)[170:340,2000:2170] gt 0) and (subImg0.(sc)[170:340,2000:2170] lt pmax), count)			;identify the good pixels
									f1 = median((subImg0.(sc)[170:340,2000:2170])[gpix], /even)

									gpix = where( (subImg1.(sc)[170:340,2000:2170] gt 0) and (subImg1.(sc)[170:340,2000:2170] lt pmax), count)			;identify the good pixels
									f2 = median((subImg1.(sc)[170:340,2000:2170])[gpix], /even)

									print, (f1 + f2)/dvar
									pause
								endif
								statusbar, k, 4*((nListB*(nListB-1)/2)+(nListF*(nListF-1)/2))*n_elements(self.chips), msg='Calculating Gain & Read Noise ...'
								++k
							endfor
						endfor
					endif
				endfor

				case type of
					0:	BEGIN
							fmed = med
							fvar = var
						END
					1:	BEGIN
							bmed = med
							bvar = var
						END
				endcase
			endfor

			;at this point I have the median of every flat and bias on this chip and the variance of the difference image between the main chip and each successive one
			;print, 'and'

			kk = 0
			for f=0,(n_elements(fmed[0,*])-2) do begin											;for each flat except the last one

				f1 = fmed[*,f]																	;the median value of the first flat
				for ff=(f+1),(n_elements(fmed[0,*])-1) do begin
					f2 = fmed[*,f+1]															;the median value of the second flat
					varf = fvar[*,f,ff]															;the variance of their difference
					for b=0,(n_elements(bmed[0,*])-2) do begin									;for each bias except the last one
						b1 = bmed[*,b]															;the median value of the first bias
						for bb=(b+1),(n_elements(bmed[0,*])-1) do begin
							b2 = bmed[*,b+1]													;the median value of the second bias
							varb = bvar[*,b,b+1]												;the variance of their difference
							dgain = ((f1 + f2) - 0*(b1 + b2)) / (varf - varb)			;calculate the gain from the above values, don't subtract off the bias if these are overscan subtracted
							if (kk eq 0) then gain = dgain else gain = [[gain], [dgain]]			;record each gain value
							dron = dgain*sqrt(varb)/sqrt(2.)									;calculate the read out noise
							if (kk eq 0) then ron = dron else ron = [[ron], [dron]]				;record each RON value
							++kk
						endfor
					endfor
				endfor
			endfor
			;gains.(c).mean 	= mean(gain, dimension=2)													;store the average gain values to the gains structure
			;gains.(c).sigma = stddev(gain, dimension=2)
			gains.(c).mean = wtd_mean(mean(gain, dimension=2), 1d/stddev(gain, dimension=2)^2d)			;the mean of all values across chiplets for each chip
			gains.(c).sigma = stddev(mean(gain, dimension=2))
			gains.(c).ron = wtd_mean(mean(ron, dimension=2), 1d/stddev(ron, dimension=2)^2d)
			gains.(c).ronsig = stddev(mean(ron, dimension=2))
		endif
	endfor
	;gains = create_struct('num0', mrdfits('supcamgain.fits',1))
	;for i=2,(10) do gains = struct_addtags(gains, create_struct('num'+strtrim(i,2), mrdfits('supcamgain.fits',i)))
	gain = gains.(0)
	for i=1,(n_tags(gains)-1) do gain = struct_append(gain, gains.(i))

	psopen, 'gain_readnoise', /encapsulated, xs=8, ys=6, /inches;, /heiles
		!p.multi=[0,1,2]
		!p.charsize = 1.3
		!p.charthick = 3
		multiplot
		ploterror, indgen(10), gain.mean[0], gain.sigma[0], ytitle='Gain [e-/ADU]', xrange=[-1,10], thick=3, errthick=3, title=sxpar(headfits(flat[0]),'DATE-OBS'), psym=10
		sharpcorners, thick=3
		multiplot
		ploterror, indgen(10), gain.ron[0], gain.ronsig[0], ytitle='Read Noise [e-]', xrange=[-1,10], thick=3, xtitle='Chips [' + self.chips[0] + ' - ' + self.chips[-1] + ']', errthick=3, psym=10
		sharpcorners, thick=3
	psclose
	self->putInDir, 'gain_readnoise.eps', 'dir_DIAG'

	;gain[i-1] = ((mean(f1) + mean(f2)) - (mean(b1) + mean(b2))) / (variance(f1 - f2) - variance(b1 - b2))

	if (self.date lt 2008.58) then begin											;multiply chips by median gain so as to not introduce offsets again
		dgain = gains.(0).mean
		for i=1,(n_tags(gains)-1) do dgain = [dgain, gains.(i).mean]
		mdgain = median(dgain, /even)
		for i=0,(n_tags(gains)-1) do begin
			gains.(i).mni = gains.(i).mean
			gains.(i).mean = mdgain
		endfor
	endif

	if (n_elements(outFile) eq 0) then outFile = 'supCamGain.fits'
	doutfile = file_search(outFile,count=noutfile)
	if (noutfile gt 0) then begin								;if this file already exists
		file_delete, outfile+'.old', /allow_nonexistent			;	delete any existing duplicates
		file_move, outfile, outfile+'.old'						;	append .old to the existing outfile
	endif
	for i=0,(n_tags(gains)-1) do mwrfits, gains.(i), outFile	;write the gains to the outfile

	return, gains
END

;-------------------------------------------------------------------------------------------

function supCamData::supCheckConstraints, files0, constraint, remove, count=count
	COMPILE_OPT idl2, HIDDEN
	if (n_elements(files0) eq 0) then stop
	files = files0
	remonly = 0
	if (n_elements(constraint) eq 0) then remonly = 1 else value = tag_names(constraint)
	if (n_elements(remove) eq 0) then remove = {OBJECT:'FAKE'}
	remvalue = tag_names(remove)
	count = 0
	for l=0,(n_elements(remvalue)-1) do begin													;for each removal constraint
		for k=0,max([0,(n_elements(value)-1)]) do begin											;for each retention constraint
			j = 0																				;counter variable
			for i=0,(n_elements(Files) - 1) do begin											;for each file
				hdr = headfits(Files[i])														;get the header
				if remonly then constraint = {NAXIS1:strcompress(sxpar(hdr, (value = 'NAXIS1')),/remove_all)}
				if ((strcompress(sxpar(hdr, value[k]),/remove_all) eq strcompress(constraint.(k), /remove_all)) and (strcompress(sxpar(hdr, remvalue[l]),/remove_all) ne strcompress(remove.(l), /remove_all))) then begin								;if this file meets the constraint
					if (j eq 0) then tlist = Files[i] else tlist = [tlist, Files[i]]			;keep it
					j++																			;and count it
				endif
			endfor
			if (n_elements(tlist) gt 0) then begin												;if any files remain then
				files = tlist																	;	reset the file list to reflect the reduced list
				undefine, tlist																	;	erase the temporary list
			endif else return, !NULL															;otherwise return nothing
		endfor
	endfor
	count = n_elements(files)
	return, files																				;return the list of files that meet all the constraints
end

;-------------------------------------------------------------------------------------------

PRO supCamData::overscan_subtract, monitor=monitor
;--Overscan Subtraction (change H* to To_RH*)
	COMPILE_OPT idl2, HIDDEN
	if (self.pre eq 'NULL') then self.pre = 'H'
	self->checkinput, /ndfilt

	post = 'To_R'

	rFiles = self->getFiles(count)																;get the self.pre files
	dummy = self->getFiles(pre=post)															;move any already-overscan-subtracted images to the working directory

	if (count gt 0) then begin
		files = self->supCheckConstraints(rFiles, {FILTER01:self.filter}, {OBJECT:'BIAS'}, count=nFiles)	;isolate files that are not biases taken with a specific filter
		if (nFiles eq 0) then stop

		dummy = file_search(post+files, count=nSuperBiasSub)

		if (nSuperBiasSub eq nFiles) then begin													;check whether the files are already bias corrected
			if ~self.batch then begin															;if the batch keyword is not set
				print
				redo = choose('Files are already bias corrected (overscan + superbias), and trimmed, redo?', 'No') ;ask if we should redo things
			endif else redo = 'NO'																;if the batch keyword is set, then don't redo
		endif else redo = 'YES'

		if ((strupcase(redo) eq 'YES') or (strupcase(redo) eq 'Y')) then begin
			filesNoGz = strarr(nFiles)															;string array to hold the above filenames w/o a .gz extension
			for i=0,(nFiles - 1) do filesNoGz[i] = strtrim(strmid(files[i],0,strpos(files[i],'.gz')),2) ;trim off the .gz extension from each filename
			dummy = file_search(post+filesNoGz, count=nOverScanSub)

			if (nOverScanSub eq nFiles) then begin												;check whether the files have just been overscan subtracted
				if ~self.batch then begin														;if the batch keyword is not set
					print
					redo = choose('Files are already bias corrected (overscan), redo?', 'No')	;then ask if we should redo the overscan subtraction
				endif else redo = 'NO'															;if the batch keyword is set, then don't redo the overscan subtraction
			endif else redo = 'YES'

			if ((strupcase(redo) eq 'YES') or (strupcase(redo) eq 'Y')) then begin				;if we need to redo the overscan subtraction
				;print
				;print, 'Subtracting Overscan...'
				self->overScanSub, files, /nopattern, monitor=monitor, statusMsg='Subtracting overscan from science frames ...'	;subtract the bias from the overscan region
			endif

			rSuperBias = self->getFiles(nSuperBias, pre='SUPERBIAS',dirPre='')					;retrieve existing superbiases

			needSuperBias = 0
			chipNames = strarr(nFiles)
			for i=0,(nFiles - 1) do chipNames[i] = strtrim(sxpar(headfits(files[i]),'DETECTOR'));record the chip names of all the science files
			chipNames = chipNames[sort(chipNames)]												;sort the chip names
			chipNames = chipNames[uniq(chipNames)]												;keep just the unique names

			if (nSuperBias gt 0) then begin														;if there are some superbiases
				SBchipNames = strarr(nSuperBias)												;make a string array to hold their names
				for i=0,(nSuperBias - 1) do SBchipNames[i] = strtrim(sxpar(headfits(rSuperBias[i]),'DETECTOR'))	;record the chip names of all the superbiases
				SBchipNames = SBchipNames[sort(SBchipNames)]									;sort the chip names
				SBchipNames = SBchipNames[uniq(SBchipNames)]									;keep just the unique names.

				foreach chip, chipNames do begin												;cycle through the science files
					ind = where(chip eq SBchipNames,nMatch)										;check whether there is a corresponding superbias
					if (nMatch eq 0) then begin													;if not, then
						needSuperBias = 1														;set the needSuperBias flag
						break																	;and break out of the loop
					endif
				endforeach
			endif else needSuperBias = 1														;set the needSuperBias flag if there are no superbiases

			if ~needSuperBias then begin
				if ~self.batch then redo = choose('Superbiases already exist, remake?', 'no') else redo = 'no'
			endif else redo = 'yes'

			if ((strupcase(redo) eq 'YES') or (strupcase(redo) eq 'Y')) then begin
				biases = self->supCheckConstraints(rfiles, {OBJECT:'BIAS'})						;identify the biases
				nBiases = n_elements(biases)
				if (nBiases gt 0) then begin													;if there are any
					biasesNoGz = strarr(nBiases)			;string array to hold the above filenames w/o a .gz extension
					for i=0,(nBiases - 1) do biasesNoGz[i] = strtrim(strmid(biases[i],0,strpos(biases[i],'.gz')),2) ;trim off the .gz extension from each filename
					dummy = file_search(post+biasesNoGz,count=nBiasesNoGz)
					if (nBiasesNoGz lt nBiases) then begin
						self->overScanSub, biases, /nopattern, monitor=monitor, statusMsg='Subtracting overscan from bias frames ...'	;subtract the bias from the overscan region
						for i=0,(nBiases - 1) do self->gZip, post+biasesNoGz[i]
					endif
				endif else stop

				if (self.date gt 2008.58) then begin
					bFiles = file_search('To_R*fits.gz',count=nBFiles)
;					if (nBFiles eq 0) then stop
					obiases = self->supCheckConstraints(bfiles, {OBJECT:'BIAS'}, count=nBias)	;isolate the biases
					if (nbias eq 0) then begin
						message, 'There are no biases here to make the superbias...', /continue, /informational
						stop
					endif

					k = 0
					foreach chip, chipNames do begin											;for each chip
						cbiases = self->supCheckConstraints(obiases, {DETECTOR:chip}, count=nCBiases)	;isolate just those biases
						if (nCBiases eq 0) then stop
						for i=0,(n_elements(cbiases)-1) do begin								;for each on
							img = mrdfits(cbiases[i],0, hdr,/silent)							;read in the file
							if (i eq 0) then hold = fltarr(n_elements(img[*,0]),n_elements(img[0,*]),n_elements(cbiases))	;first time through make a holding array
							hold[*,*,i] = img													;store file to holding array
							statusbar, k, n_elements(chipNames)*n_elements(cBiases), msg='Creating SuperBiases...'
							++k
						endfor
						;superbias = median(hold,dimension=2,/even)								;take the median value at each pixel
						;superbias = median(hold,dimension=3,/even)								;take the median value at each pixel
						superbias = median(hold)								;take the median value at each pixel
						bpix = where(finite(superbias) eq 0, nbpix)								;check for bad pixels
						if (nbpix gt 0) then superbias[bpix] = self.badPixValue					;mark any that exist
						;spawn, 'mv SUPERBIAS_'+chip+'.fits SUPERBIAS_'+chip+'.fits.old'		;rename previous versions
						self->fileDelete, 'SUPERBIAS_'+chip+'.fits'								;delete any old versions
						SXADDPAR, hdr, 'OBJECT', 'SUPERBIAS', /SaveComment						;change object type to SUPERBIAS
						mwrfits, superbias, 'SUPERBIAS_'+chip+'.fits',hdr						;save the superbias
					endforeach
				endif
			endif else spawn, 'ls To_R*fits', bfiles											;find the resulting files from both subtractions

			if (self.date gt 2008.58) then begin
				k = 0
				foreach chip, chipNames do begin												;for each chip
					bFiles = file_search('To_R*fits', count=nBFiles)							;find the resulting files from both subtractions
					ofiles = self->supCheckConstraints(bfiles, {FILTER01:self.filter}, {OBJECT:'BIAS'}, count=nOFiles)	;isolate non-biases
					if (nOFiles eq 0) then stop
					superbias = mrdfits('SUPERBIAS_'+chip+'.fits',0,/silent)					;read in the superbias
					cfiles = self->supCheckConstraints(ofiles, {DETECTOR:chip})					;identify the non-biases
					for i=0,(n_elements(cFiles)-1) do begin										;for each non-bias
						img = mrdfits(cfiles[i],0,hdr,/silent)									;read in file
						;bpix = where( (img eq -32768) or (superbias eq -32768), nbpix)			;identify bad pixels
						bPix = where( ~finite(img), nBPix)										;identify bad pixels
						img = img - superbias													;subtract superbias
						if (nbpix gt 0) then img[bpix] = self.badPixValue						;mark bad pixels
						sxaddhist, 'supcampipe.pro: SUPERBIAS subtracted.', hdr					;add history
						;spawn, 'mv ' + cfiles[i] + ' ' + cfiles[i] + '.old'					;move original file
						self->fileMove, cFiles[i], cFiles[i]+'.old'
						mwrfits, img, cfiles[i], hdr											;save new file
						self->fileDelete, cFiles[i]+'.old'										;delete original file
						self->fileDelete, cFiles[i]+'.gz'										;delete original file
						self->gZip, cFiles[i]													;gZip the new files
						statusbar, k, n_elements(chipNames)*n_elements(cFiles), msg='Subtracting SuperBiases ...'
						++k
					endfor
					self->putInDir, 'SUPERBIAS_'+chip+'.fits', 'SUPERBIAS'						;move new superbiases to SUPERBIAS directory
				endforeach
			endif
		endif
	endif else print, 'No '+self.pre+'*.fits files to overscan subtract.'
	print
	self->putInDir, self.pre+'*.fits.gz', 'dir'+self.pre										;move the input files to a directory
	self->putInDir, post+'*.fits.gz', 'dir'+post												;move the output files to a directory
	self.pre = post
END

;-------------------------------------------------------------------------------------------

PRO supCamData::overScanSub, files, monitor=monitor, nopattern=nopattern, pre=pre, sm=sm, statusMsg=statusMsg
	COMPILE_OPT idl2, HIDDEN
	if (n_elements(sm) eq 0) then sm = 5
	;monitor = 1
	porder = 7
	nFiles = n_elements(files)
	if (n_elements(statusMsg) ne 0) then msg = statusMsg else msg = 'Subtracting overscan ...'

	for f=0,(nFiles-1) do begin
		img = mrdfits(files[f], 0, hdr, /silent)
		img = img + sxpar(hdr, 'BZERO')
		date = fix(strtrim(strsplit(sxpar(hdr, 'DATE-OBS'),'-',/extract)))
		date = date[0] + (date[1]-1)/12. + date[2]/365.

		if (((date gt 2008.58) and (median(img) lt 50000.)) or ((date le 2008.58) and (median(img) lt 35000.))) then begin					;check to see if the flat is near the non-linearity regime
			bias = float(img)
			;overscan
			if (date gt 2008.58) then begin
				nchiplet = 4
				noverscan = 2
			endif else begin
				nchiplet = 1
				noverscan = 1
			endelse

			xmin = intarr(nchiplet)
			xmax = intarr(nchiplet)
			ymin = intarr(nchiplet)
			ymax = intarr(nchiplet)

			for c=1,nchiplet do begin
				if (date gt 2008.58) then begin
					xemin = sxpar(hdr, 'S_EFMN' + strtrim(c,2) + '1') - 1
					xemax = sxpar(hdr, 'S_EFMX' + strtrim(c,2) + '1') - 1
					yemin = sxpar(hdr, 'S_EFMN' + strtrim(c,2) + '2') - 1
					yemax = sxpar(hdr, 'S_EFMX' + strtrim(c,2) + '2') - 1

					xomin = sxpar(hdr, 'S_OSMN' + strtrim(c,2) + '1') - 1
					xomax = sxpar(hdr, 'S_OSMX' + strtrim(c,2) + '1') - 1
					yomin = sxpar(hdr, 'S_OSMN' + strtrim(c,2) + '2') - 1
					yomax = sxpar(hdr, 'S_OSMX' + strtrim(c,2) + '2') - 1
				endif else begin
					xemin = sxpar(hdr, 'EFP-MIN1') - 1
					xemax = xemin + sxpar(hdr, 'EFP-RNG1') - 1
					yemin = sxpar(hdr, 'EFP-MIN2') - 1
					yemax = yemin + sxpar(hdr, 'EFP-RNG2') - 1

					if (xemin eq 0) then begin
						xomin = xemax + 1
						xomax = sxpar(hdr, 'PRD-RNG1') - 1
					endif else begin
						xomin = 0
						xomax = xemin - 1
					endelse
					if (yemin eq 0) then begin
						yomin = yemax + 1
						yomax = sxpar(hdr, 'PRD-RNG2') - 1
					endif else begin
						yomin = 0
						yomax = yemin - 1
					endelse
                                     endelse


				xmin[c-1] = xemin
				xmax[c-1] = xemax
				ymin[c-1] = yemin
				ymax[c-1] = yemax

				for i=0,(noverscan - 1) do begin

					if (i eq 0) then reg = img[xomin:xomax, *] else begin		;isolate the x-overscan region
						reg = img[xemin:xemax,yomin:yomax]						;isolate the y-overscan region
						reg = transpose(reg)									;transpose this in order to use the same code as the x-overscan
					endelse

					mn = mean(reg, dimension=1)									;calculate the mean of each row
					mn = transpose(rebin(mn, n_elements(reg[0,*]),n_elements(reg[*,0])))		;create 2d image with appropriate mean value placed at each pixel
					sig = stddev(reg, dimension=1)								;calculate the standard deviation of each row
					sig = transpose(rebin(sig, n_elements(reg[0,*]),n_elements(reg[*,0])))		;create 2d image with appropriate sigma value placed at each pixel
					sigma = (reg - mn) / sig									;record the number of sigma away each pixel is
					bPix = where((sigma gt 3.0) or (finite(sigma) eq 0), count)	;identify outliers
					if (count gt 0) then reg[bPix] = !values.f_nan				;mark those pixels with a NaN

					num = total(reg*0+1., 1, /nan)								;calculate the number of viable rows
					mn = mean(reg, dimension=1, /nan)							;recalculate the mean of each row
					sig = stddev(reg, dimension=1, /nan)/sqrt(num)				;calculate the standard deviation of the mean for each row

					bpix = where((finite(mn) eq 0) or (finite(sig) eq 0) or (sig eq 0.0), nbpix)						;identify bad pixels
					if (nbpix gt 0) then begin
						mn[bPix] = mean(mn, /nan)
						sig[bPix] = 9999.
					endif

					ind = indgen(n_elements(mn))
					if (i eq 0) then ind = ind[min([yemin, yomin]):max([yemax, yomax])]

					if (N_elements(monitor) ne 0) then begin
						erase
						!p.multi=[0,1,2]
						multiplot
						plot, ind, mn, yrange=median(mn)+[-1,1]*3.*robust_sigma(mn);minmax(mn)
					endif
					mn = mn[ind]
					sig = sig[ind]

					nsig = abs(mn - median(mn))/sig

					;plot, ind, sig, yrange=[0,1]
					bPix = where( (nsig gt 5.*robust_sigma(nsig)) , count)
					if (count gt 0) then sig[bPix] = sig[bpix]*nsig[bpix]
					;oplot, ind, sig, color=fsc_color('red')

					dummy = poly_fit(ind, mn, porder, measure_errors=sig[ind], yfit=yfit, yband=yband, yerror=yerror)
					if (N_elements(monitor) ne 0) then oplot, ind, yfit, color=fsc_color('red')

					;if (N_elements(monitor) ne 0) then print, yerror
					for k=0,10 do begin
						dummy = poly_fit(ind, mn, porder, measure_errors=sqrt(sig[ind]^2d + (yfit - mn)^2d), yfit=yfit, yband=yband, yerror=yerror)
						if (N_elements(monitor) ne 0) then oplot, ind, yfit, color=fsc_color('red')
						++k

					endfor
					if (N_elements(monitor) ne 0) then oplot, ind, yfit, color=fsc_color('green')
					;if (N_elements(monitor) ne 0) then print, yerror
					;if (N_elements(monitor) ne 0) then pause

	;				if (yerror gt 2) then stop		;check for a bad fit
	;need a diagnostic plot for these
	;maybe a histogram of residuals for the 99% percentile

					if (n_elements(nopattern) eq 0) then begin
						tsig = sig[ind]
						bpix = where(tsig gt (5.*robust_sigma(tsig) + biweight_mean(tsig)), nbpix)

						tmn = mn
						if (nbpix gt 0) then tmn[bpix] = !values.f_nan
						yfit0 = yfit
						yfit = smooth(tmn,sm, /edge_truncate, /nan)
						bpix = where(finite(yfit) eq 0, nbpix)
						if (nbpix gt 0) then begin
							yfit[bpix] = interpol(yfit, ind, bpix, /nan)
							;stop
						endif

						if (N_elements(monitor) ne 0) then begin
							multiplot
							plot, ind, yfit0 - mn, yrange=[-1,1]*50.*robust_sigma(sig), psym=-3
							oplot, ind, yfit - mn, color=fsc_color('red'), psym=-3, thick=2
							;pause
							;erase
						endif

					endif

					if (n_elements(hold) eq 0) then begin
						ntags = '0'
						dummy = struct_addtags({a:0d}, {b:0d})
						undefine, dummy
					endif else ntags = strtrim(n_tags(hold),2)
					hold = struct_addtags(hold, create_struct('num'+ntags, {file:files[f], ind:ind, yfit:yfit, mn:mn}))

					;yfit = reverse(yfit)
					bias = finite(img)*0.
					if (i eq 0) then bias[min([xemin,xomin]),min([yemin, yomin])] = transpose(rebin(yfit, n_elements(yfit), (max([xemax,xomax]) - min([xemin,xomin]) + 1) ) )	$
								else bias[xemin,yemin] = rebin(yfit, n_elements(yfit), (yemax - yemin + 1))

;					bpix = where( finite(img) eq 0, nbad)
;					if (nbad gt 0) then stop

					bpix = where( finite(bias) eq 0, nbad)
					if (nbad gt 0) then stop
					img -= bias
					bpix = where( finite(img) eq 0, nbad)
					if (n_elements(nbad) gt 0) then img[bpix] = self.badPixValue
					;pause
				endfor
			endfor

			sxaddhist, 'overscansub.pro : Overscan subtracted', hdr	;record history

			if (date gt 2008.58) then begin
				chipOrder = strtrim((indgen(4)+1)[sort(xmin)],2)
				xmax = xmax[sort(xmin)]
				ymin = ymin[sort(xmin)]
				ymax = ymax[sort(xmin)]
				xmin = xmin[sort(xmin)]

				newImg = fltarr(total(xmax - xmin + 1), max(ymax[0] - ymin[0]) + 1)

				for c=0,3 do begin

					newxmin = (total((xmax - xmin)[0:c-1]) + c)*min([c,1])
					newymin = ymin[c]-min(ymin)
					newImg[newxmin,newymin] = img[xmin[c]:xmax[c], ymin[c]:ymax[c]]
					SXADDPAR, hdr, 'BZERO', 0., /SaveComment
					SXADDPAR, hdr, 'S_EFMN' + chiporder[c] + '1', newxmin + 1, /SaveComment
					SXADDPAR, hdr, 'S_EFMN' + chiporder[c] + '2', newymin + 1, /SaveComment
					SXADDPAR, hdr, 'S_EFMX' + chiporder[c] + '1', newxmin + 1 + (xmax[c] - xmin[c]), /SaveComment
					SXADDPAR, hdr, 'S_EFMX' + chiporder[c] + '2', newymin + 1 + (ymax[c] - ymin[c]), /SaveComment
				endfor
			endif else begin

				newImg = float(img[xmin:xmin+2048-1, ymin:ymin+4096-1])
				SXADDPAR, hdr, 'BZERO', 0., /SaveComment
				SXADDPAR, hdr, 'EFP-MIN1', 1, /SaveComment
				SXADDPAR, hdr, 'EFP-MIN2', 1, /SaveComment
				SXADDPAR, hdr, 'EFP-RNG1', 2048, /SaveComment
				SXADDPAR, hdr, 'EFP-RNG2', 4096, /SaveComment

				;for c=0,3 do begin
				;	SXADDPAR, hdr, 'S_EFMN' + strtrim(c+1,2) + '1', c*512 + 1, /SaveComment
				;	SXADDPAR, hdr, 'S_EFMN' + strtrim(c+1,2) + '2', 1, /SaveComment
				;	SXADDPAR, hdr, 'S_EFMX' + strtrim(c+1,2) + '1', (c+1)*512, /SaveComment
				;	SXADDPAR, hdr, 'S_EFMX' + strtrim(c+1,2) + '2', 4096, /SaveComment
				;endfor
			endelse

			sxaddhist, 'overscansub.pro : Trimmed', hdr	;record history

			if (n_elements(pre) eq 0) then pre = ''

			outfile = pre + 'To_R' + files[f]
			dfile = file_search(outFile, count=ndfile)
			;if (ndfile gt 0) then file_delete, outfile
			;if (ndfile gt 0) then file_move, outfile, outfile+'.old', /overwrite
			;spawn, 'mv ' + outFile + ' ' + outFile + '.old'

			if ((ipos = strpos(outFile,'.gz')) gt 0) then outFile = strmid(outFile,0,ipos)

			mwrfits, newImg, outFile, hdr, /create

		endif
		statusbar, f, nfiles, msg=msg	;msg='Subtracting overscan...'
		;print, outfile
		;pause
             endfor


	ntags = n_tags(hold)
	sq = (sqrt(float(ntags)))
	cleanplot, /silent
	erase
	file = file_search('diag_overscan.eps', count=count)
	if (count gt 0) then file_move, 'diag_overscan.eps', 'diag_overscan_old.eps', /overwrite
	psopen, 'diag_overscan', /encapsulated, xs=10, ys=10, /inches
	!p.multi=[0,ceil(sq),floor(sq)]
	file = strarr(ntags)
	for i=0,(n_tags(hold)-1) do begin
		file[i] = hold.(i).file
		d = arr_struct(struct_trimtags(hold.(i), except='FILE'))
		d = d[sort(abs(d.yfit - d.mn))]
		numd = n_elements(d.yfit)
		ind = indgen(numd)
		ind = ind[0:round(numd*0.999)-1]
		multiplot
		plothist, (d.yfit - d.mn)[ind], bin=0.1, xrange=[-5,5], xstyle=1
;		legend, strtrim(i,2), /top, /left, box=0
	endfor
	psclose
	if (count gt 0) then file_move, 'diag_overscan.dat', 'diag_overscan_old.dat', /overwrite
	openw, lun, 'diag_overscan.dat', /get_lun
	printf, lun, [transpose(strtrim(indgen(n_elements(file)),2)), transpose(file)]
	free_lun, lun
	self->putInDir, 'diag_overscan*', 'dir_DIAG'
END

;-------------------------------------------------------------------------------------------

PRO supCamData::make_flats, sm=sm
;--Make Flat Field Frames (create obj_mflat*.fits files)
	COMPILE_OPT idl2, HIDDEN

	if (self.pre eq 'NULL') then self.pre = 'To_R'
	if (self.pre ne 'To_R') then self.pre = choose('Are you sure you want this prefix?', self.pre)

	self->checkinput, /ndfilt, /ndflat

	;calculate gain
	self->get_gain, gain=gainInfo

	gains = dblarr(n_tags(gainInfo))
	rdNzs = dblarr(n_tags(gainInfo))
	for i=0,(n_tags(gainInfo)-1) do begin
		gains[i] = mean(gainInfo.(i).mean)
		rdNzs[i] = mean(gainInfo.(i).ron)
	endfor
	gain = mean(gains)
	readNoise = mean(rdNzs)

	print
	print, 'Creating Flat Fields...'

	rfiles = self->file_check(nFiles)															;create a 'pre' directory or retrieve files from it

	;file_mkdir, 'dir'+pre
	;file_move, 'dir'+pre+'/'+pre+'*.fits', './'
	;rFiles = file_search(pre+'*.fits')				;make a list of all the gTo_RH* files

	case self.flattype of
		'DOMEFLAT'		:	flattypes = 'DOMEFLAT'												;if flattype is DOMEFLAT then make those
		'SKYFLAT'		:	flattypes = ['DOMEFLAT', 'SKYFLAT']									;if flattype is SKYFLAT then make those and DOMEFLATS
		'SUPERSKYFLAT'	:	BEGIN																;if flattype is SUPERSKYFLAT then check for existence of the other flattypes
								flats = self->supCheckConstraints(rFiles, {FILTER01:self.filter, OBJECT:'SKYFLAT'})		;identify domeflats
								if (n_elements(flats) eq 0) then begin
									flats = self->supCheckConstraints(rFiles, {FILTER01:self.filter, OBJECT:'DOMEFLAT'})	;or skyflats in the absence of domeflats
									if (n_elements(flats) eq 0) then stop else flattypes = 'DOMEFLAT'
								endif else flattypes = ['DOMEFLAT', 'SKYFLAT']
							END
	ENDCASE

	foreach flatType, flatTypes do begin														;make each flatfield type

		rFlats = file_search(flatType + '_' + self.filter + '/obj_mflat_'+self.chips[0]+'.fits', count=nMatch)

		if nMatch then begin
			if ~self.batch then remake = choose(flatType+'s already exist. Remake them?', 'no') else remake = 'no'
		endif else remake = 'yes'

		if (strupcase(remake) eq 'YES') then begin

			self->mkDir, self.flatType + '_' + self.filter

			flats = self->supCheckConstraints(rFiles, {FILTER01:self.filter, OBJECT:flatType}, count=nFlats)
			nChips = n_elements(self.chips)
			nExp = nFlats/nChips
			medVals = dblarr(nExp)
			nVert = 20
			nHorz = 2
			locMed = dblarr(nChips,nExp,4,nVert,nHorz)
			sigL = locMed
			sigH = locMed

			;find the median value of each mosaic and chip, and the local
;			if 1 then begin
;				k = 0l
;				for ixp=0,(nExp - 1) do begin
;					index = 0l
;					for ic=0,(nchips - 1) do begin
;						img = mrdfits(flats[ixp*nchips+ic], 0, header, /silent)

;						ysize = sxpar(header,'NAXIS2')
;						if (ic eq 0) then holdpix = fltarr(long(sxpar(header,'NAXIS1'))*ysize*nchips) + !values.f_nan			;make an array to hold all the pixels

;						xMin = [sxpar(header, 'S_EFMN11'), sxpar(header, 'S_EFMN21'), sxpar(header, 'S_EFMN31'), sxpar(header, 'S_EFMN41')] - 1l	;min x pix of each subchip
;						xMax = [sxpar(header, 'S_EFMX11'), sxpar(header, 'S_EFMX21'), sxpar(header, 'S_EFMX31'), sxpar(header, 'S_EFMX41')] - 1l	;max x pix of each subchip
;						for ia=0,3 do begin														;for each subchip
;							if (ia eq 0) then index_loc = index
;							subImg = img[xmin[ia]:xmax[ia],*]
;							nsubImg = n_elements(subImg)*1l
;							holdpix[index] = subImg[lindgen(nsubImg)]
;							index = index + nsubImg
;							vert = round(findgen(nvert + 1)/nvert*ysize)
;							for iy=0,(nvert - 1) do begin
;								subImg_y = subImg[*,vert[iy]:(vert[iy+1]-1l)]
;								horz = round(findgen(nHorz + 1)/nHorz*n_elements(subImg_y[*,0]))
;								for ix=0,(nHorz - 1) do begin
;									subImg_x = subImg_y[horz[ix]:(horz[ix+1]-1l),*]
;									subImg_x = subImg_x[where(finite(subImg_x) and (subImg_x gt 500) and (subImg_x lt 60000),count)]
;									locMed[ic,ixp,ia,iy,ix] = median(subImg_x,/even)
;									sigL[ic,ixp,ia,iy,ix] = robust_sigma(subImg_x[where(subImg_x lt locMed[ic,ixp,ia,iy,ix])] - locMed[ic,ixp,ia,iy,ix],/zero)
;									sigH[ic,ixp,ia,iy,ix] = robust_sigma(subImg_x[where(subImg_x ge locMed[ic,ixp,ia,iy,ix])] - locMed[ic,ixp,ia,iy,ix],/zero)
;									undefine, subImg_x
;									statusbar, k, nChips*nExp*4*nVert*nHorz, msg='Calculating noise characteristics for flat-field frames ...'
;									++k
;								endfor
;								undefine, subImg_y
;							endfor
;							undefine, subImg
;						endfor
;						undefine, img
;					endfor
;					holdpix = holdpix[where(finite(holdpix) and (holdpix gt 500) and (holdpix lt 60000),count)]
;					medVals[ixp] = median(holdPix,/even)
;					undefine, holdpix
;				endfor
;				mwrfits, medVals, 'diag_medVals.fits', /create
;				mwrfits, locMed, 'diag_locMed.fits', /create
;				mwrfits, sigL, 'diag_sigL.fits', /create
;				mwrfits, sigH, 'diag_sigH.fits', /create
;			endif else begin
;				medVals = mrdfits('diag_medVals.fits',0,/silent)
;				locMed = mrdfits('diag_locMed.fits',0,/silent)
;				sigL = mrdfits('diag_sigL.fits',0,/silent)
;				sigH = mrdfits('diag_sigH.fits',0,/silent)
;			endelse

			k = 0l
			;sigma clip each exposure, then median combine
			for ic=0,(nchips - 1) do begin
				chip = self.chips[ic]
				fchip = self->supCheckConstraints(flats, {DETECTOR:chip}, count=nExp)
				header = headfits(fchip[0])
				xMin = [sxpar(header, 'S_EFMN11'), sxpar(header, 'S_EFMN21'), sxpar(header, 'S_EFMN31'), sxpar(header, 'S_EFMN41')] - 1l	;min x pix of each subchip
				xMax = [sxpar(header, 'S_EFMX11'), sxpar(header, 'S_EFMX21'), sxpar(header, 'S_EFMX31'), sxpar(header, 'S_EFMX41')] - 1l	;max x pix of each subchip
				xsize = sxpar(header, 'NAXIS1')
				ysize = sxpar(header, 'NAXIS2')
                                medVals = dblarr(nExp)
                                nVert = 20
                                nHorz = 2
                                locMed = dblarr(nChips,nExp,4,nVert,nHorz)
                                sigL = locMed
                                sigH = locMed
				nsig = [4.,4.]
				grow = 0

				if (n_elements(rows) gt 2) then if (rows[-2] ge rows[-1]) then rows = rows[0:-2]

				flatpix = fltarr(xsize,ysize)
				hold = fltarr(xsize, ysize,nexp)+!values.f_nan
				mask = intarr(xsize, ysize,nexp)

				for ixp=0,(nexp - 1) do begin
					img = mrdfits(fchip[ixp],0,hdr,/silent)
					vert = round(findgen(nvert + 1)/nvert*ysize)
					for ia=0,3 do begin
						subImg = img[xmin[ia]:xmax[ia],*]
						for iy=0,(nvert - 1) do begin
							subImg_y = subImg[*,vert[iy]:(vert[iy+1]-1l)]
							horz = round(findgen(nHorz + 1)/nHorz*n_elements(subImg_y[*,0]))
							for ix=0,(nHorz - 1) do begin
								subImg_x = subImg_y[horz[ix]:(horz[ix+1]-1l),*]
                                                                locMed[ic,ixp,ia,iy,ix] = median(subImg_x,/even)
								sigL[ic,ixp,ia,iy,ix] = robust_sigma(subImg_x[where(subImg_x lt locMed[ic,ixp,ia,iy,ix])] - locMed[ic,ixp,ia,iy,ix],/zero)
								sigH[ic,ixp,ia,iy,ix] = robust_sigma(subImg_x[where(subImg_x ge locMed[ic,ixp,ia,iy,ix])] - locMed[ic,ixp,ia,iy,ix],/zero)

								bPix = where( (subImg_x lt (locMed[ic,ixp,ia,iy,ix] - nsig[0]*sigL[ic,ixp,ia,iy,ix])) or (subImg_x ge (locMed[ic,ixp,ia,iy,ix] + nsig[1]*sigH[ic,ixp,ia,iy,ix])), nbad)
                                                                ;if (iy eq 10) then stop
								if (nbad gt 0) then begin

									subImg_x = maskgrow(subImg_x, bPix, grow, !values.f_nan)

								endif
								subImg_y[horz[ix],0] = subImg_x
								undefine, subImg_x
								statusbar, k, nChips*nExp*4*nVert*nHorz, msg='Masking bad pixels & median combining flat-field frames...'
								++k
							endfor
							subImg[0,vert[iy]] = subImg_y
							undefine, subImg_y
						endfor
						img[xmin[ia],0] = subImg
						undefine, subImg
                                             endfor
					medVals[ixp] = median(img[where(finite(img) and img gt 500 and img lt 60000,count)],/even)
					hold[0,0,ixp] = img/medVals[ixp]
                                     endfor
				flatpix[0,0] = median(hold,dim=3,/even)/gain

				bPix = where(~finite(flatPix),nBad)
				if (nBad gt 0) then flatPix[bPix] = self.badPixValue
				flatPix[0:9,*] = self.badPixValue
				flatPix[*,0:9] = self.badPixValue
				flatPix[-10:*,*] = self.badPixValue

				;add info to the header?
				sxaddhist, 'supcampipe.pro: '+flattype+' constructed', header	;add some history
				sxaddhist, 'supcampipe.pro: Flat divided by gain value: '+strtrim(gain,2), header	;add some history
				SXADDPAR, header, 'GAIN', float(gain), /savecomment
				SXADDPAR, header, 'RDNOISE', float(readNoise), /savecomment


				;mask the AG probe
				AGX = sxpar(header, 'S_AG-X')
	   			j_limit = 5000
			    if (self.date lt 2008.58) then begin
			    	if ((chip eq 'si006s') or (chip eq 'si002s') or (chip eq 'w6c1')) then j_limit = AGX*60 - 2300
			    	if ((chip eq 'w67c1') or (chip eq 'w7c3')) then j_limit = AGX*60 - 2000
			    endif else begin
			    	if ((chip eq 'clarisse') or (chip eq 'fio') or (chip eq 'kiki')) then j_limit = AGX*60 - 2300
			    	if ((chip eq 'chihiro') or (chip eq 'nausicaa')) then j_limit = AGX*60 - 2000
			    endelse
				if (j_limit le (ySize-1)) then flatPix[*,j_limit:*] = self.badPixValue
				outFile = self.flatType + '_' + self.filter + '/obj_mflat_'+chip+'.fits'
				mwrfits, flatPix, outFile, header, /create
                                ;spawn, 'mask_for_AGX '+outfile+'
                                ;'+outfile+' '+strtrim(j_limit,2)+'
                                ;-32768'

			endfor
		endif
	endforeach	;will end with flattype = SKYFLAT if it exists, otherwise it will be DOMEFLAT

	self->putInDir, self.pre+'*.fits', 'dir'+self.pre											;move all the To_RH* files to a dir

	if (strupcase(self.flatType) eq 'SUPERSKYFLAT') then begin									;if superskyflats are specified
		stop
		;i'm not sure the code below is correct..., should flatType be self.flatType?
		;supFinalSuperFlats, self.filter, flatType, sm=sm										;make the final super sky flats
		spawn, 'mkdir ' + self.flatType + '_' + self.filter										;make flats folder for this filter
		spawn, 'mv obj_mflat*.fits ' + self.flatType + '_' + self.filter + '/'					;put the flats in the folder
	endif

	;cleanup some extraneous files
	;delList = ['mnah*fits', 'ssbtmp*fits', 'tmp*', 'blankmap*', 'check*','mask_mkflat.cat','flatList']
	;for i=0,(n_elements(delList)-1) do begin
	;	files = file_search(delList[i], count=nFiles)
	;	if (nFiles gt 0) then file_delete, files
	;endfor
	foreach file, 'diag_'+['sigL','sigH','medVals','locMed']+'.fits' do self->putInDir, file, 'dir_DIAG';move all the To_RH* files to a dir
	print
END

;-------------------------------------------------------------------------------------------

PRO supCamData::flatfield_data
;--Flat Field Object Data (change To_RH* to fTo_RH*)
	COMPILE_OPT idl2, HIDDEN
	if (self.pre eq 'NULL') then self.pre = 'To_R'

	self->checkinput, /ndfilt, /ndflat, /ndobj

	print
	print, 'Flat Fielding Data...'

	;if (date lt 2008.58) then begin
	;	print
	;	print, 'These data have photometric flatfields.'
	;	if (choose('Use improved flatfields?','yes') eq 'yes') then begin
	;		flats = file_search('capak_flats/' + filter + '/obj_mflat*', count=nfiles)
	;		if (nfiles eq 0) then stop $
	;						 else spawn, 'ls -1 capak_flats/' + filter + '/obj_mflat* > finalFlats'	;make a list of the combined flat field frames
	;	endif else spawn, 'ls -1 ' + flatType + '_' + filter + '/obj_mflat* > finalFlats'	;make a list of the combined flat field frames
	;endif else spawn, 'ls -1 ' + flatType + '_' + filter + '/obj_mflat* > finalFlats'	;make a list of the combined flat field frames


	;if (date lt 2008.58) 	then flats = file_search('capak_flats/' + filter + '/obj_mflat*', count=nfiles) $
	;						else flats = file_search(flatType + '_' + filter + '/obj_mflat*', count=nfiles)
	flats = file_search(self.flatType + '_' + self.filter + '/obj_mflat*', count=nFiles)

	rfiles = self->getFiles(nFiles)																;make a list of all the gTo_RH* files

	files = self->supCheckConstraints(rFiles, {FILTER01:self.filter, OBJECT:self.object}, count=count)
	if (count eq 0) then stop

	k = 0
	foreach chip, self.chips do begin
		fchip = self->supCheckConstraints(Files, {FILTER01:self.filter, OBJECT:self.object, DETECTOR:chip}, count=nFiles)	;identify the science frames for a particular chip
		fflat = self->supCheckConstraints(flats, {DETECTOR:chip}, count=nFlats)					;identify the flat frames for a particular chip
		if (nFiles eq 0) then stop
		if (nFlats eq 0) then stop
		flat = mrdfits(fflat,0,fhdr,/silent)													;get the flat
		readNoise = sxpar(fhdr,'RDNOISE')
		if (readNoise eq 0) then stop

		for i=0,(nFiles - 1) do begin															;for each science frame
			img = mrdfits(fchip[i],0,hdr,/silent)												;read it in
			;bpix = where( (img eq -32768) or (flat eq -32768), nbpix)							;find bad pixels in the science and flatfield frames
			bPix = where( ~finite(img) or ~finite(flat), nBPix)									;find bad pixels in the science and flatfield frames
			if (nbpix eq 0) then stop															;
			img = img / flat																	;divide the science frame by the flat

	;		flatmod = file_search('flatmod_'+chip+'.fits',count=nflatmod)
	;		if (nflatmod eq 0) then stop
	;		img = img / mrdfits(flatmod,0,/silent)

			if (nbpix gt 0) then img[bpix] = self.badPixValue									;replace bad pixels by -32768

			sxaddpar, hdr, 'GAIN', 1., /savecomment
			sxaddpar, hdr, 'RDNOISE', readNoise, /savecomment
			sxaddhist, 'supcampipe.pro: Flat-fielded.', hdr
			sxaddpar, hdr, 'SATURATE', self.satLevel/median(flat), 'Saturation/non-linearity level', before='HISTORY'	;add in saturation level keyword (~where non-linearity sets in)

			outfile = 'f'+fchip[i]
			if ((ipos = strpos(outFile,'.gz')) gt 0) then outFile = strmid(outFile,0,ipos)

			mwrfits, img, outFile, hdr, /create													;write the new image to file
			self->gZip, outFile

			statusbar, k*nFiles + i, nFiles*10, msg='Dividing by flatfield frames...'
		endfor
		++k
	endforeach

	;supcam_write2file, files, 'sciList'
	;spawn, 'ffield.csh finalFlats sciList'														;apply the flat field correction

	self->putInDir, self.pre+'*.fits.gz', 'dir'+self.pre										;move all the To_RH* files to a dir
	dummy = file_search('flatmod*.fits',count=nMatch)
	if (nMatch gt 0) then self->putInDir, 'flatmod*.fits', 'dir_flatmod'						;move all the To_RH* files to a dir
	rfiles = self->getFiles1(pre='flatmod', dirPre='dir_flatmod')								;
	;add in the code to use flatmods

	self.pre = 'f'
	self->mkDir, 'dir'+self.pre
	self->putInDir, self.pre+'*.fits.gz', 'dir'+self.pre										;move all the To_RH* files to a dir
END

;-------------------------------------------------------------------------------------------

PRO supCamData::clean_cosmics
;--Run L.A. Cosmic on the flat-fielded images to remove cosmic-rays.
	self->checkinput, /ndfilt, /ndobj
	exceptState = !EXCEPT
	!EXCEPT = 0
	compile_opt idl2, hidden
	if (self.pre eq 'NULL') then self.pre = 'f'
	post = 'z'
	rfiles = self->getfiles(nFiles)
	if (nFiles eq 0) then stop
	files = self->supCheckConstraints(rFiles, {FILTER01:self.filter, OBJECT:self.object}, count=nFiles)
	if (nFiles eq 0) then stop
	for i=0,(nFiles - 1) do begin
		hdr = headfits(files[i])
		readNoise = sxpar(hdr,'RDNOISE')
		gain = sxpar(hdr,'GAIN')
		outFile = post + strtrim(strmid(files[i],0,(ipos = strpos(files[i],'.gz'))),2)
		if (ipos eq -1) then stop
		;newIm = jaa_la_cosmic(files[i], outlist=outFile, gain=gain, readn=readNoise, bPix=bPix, /verbose) ;replace cosmic ray pixels with NaNs
		newIm = jaa_la_cosmic(files[i], outlist=outFile, gain=gain, readn=readNoise, bPix=bPix, /quiet) ;replace cosmic ray pixels with NaNs
		if (n_elements(bPix) gt 0) then begin
			saturate = sxpar(hdr, 'SATURATE')
			if (saturate eq 0) then stop
			img = mrdfits(files[i],0,/silent)					;read in the original image
			notSatPix = where(img[bPix] lt saturate*self.nonLinFrac, nNotSatPix, complement=satPix, nComplement=nSatPix)
												;find CR pixels that are not saturated
			if (nSatPix gt 0) then newIm[bPix[satPix]] = img[bPix[satPix]]		;restore any saturated pixels
			if (nNotSatPix gt 0) then bPix = bPix[notSatPix] else undefine, bPix	;remove saturated pixels from the bad pixel array
			undefine, img
		endif
		if ((self.crGrow gt 0) and (n_elements(bPix) gt 0)) then $			;if there are cosmic rays, and the grow parameter is set
			newIm = maskgrow(newIm, bPix, self.crGrow, !values.f_nan)		;then grow the masked pixel mask by crGrow
		sxaddhist, 'supcamdata__define.pro: Pixels affected by cosmic rays replated with NaNs.', hdr	;add history
		mwrfits, newIm, outFile, hdr, /create						;write the new file
		statusbar, i, nFiles, msg='Removing cosmic rays ...'
	endfor
	self->putInDir, self.pre+'*.fits.gz', 'dir'+self.pre					;move all the To_RH* files to a dir
	self->gZip, post+'*.fits', msg='GZipping output files ...'
	self->putInDir, post+'*.fits.gz', 'dir'+post						;move all the To_RH* files to a dir
	self.pre = post
	!EXCEPT = exceptState
END

;-------------------------------------------------------------------------------------------

PRO supCamData::scale_chips
;--Multiplicatively scale chips (change from gqfTo_RH* to AgqfTo_RH*)
	compile_opt idl2, hidden
	if (self.pre eq 'NULL') then self.pre = 'z'
	print
	print, 'Correcting for chip-to-chip QE/Gain differences ...'
	files = self->getfiles(nFiles)
        if (self.date lt 2008.58) then upgrade = 0 else upgrade = 1
	if (nFiles gt 0) then chip2chipqe_fast, files, post, upgrade else stop
	self->putInDir, self.pre+'*.fits.gz', 'dir'+self.pre										;move all the To_RH* files to a dir
	self->gZip, post+'*.fits'
	self->putInDir, post+'*.fits.gz', 'dir'+post												;move all the To_RH* files to a dir
	self->putInDir, 'diag_chipscaling.eps', 'dir_DIAG'
	self->putInDir, 'diag_supcamchipscale.fits', 'dir_DIAG'
	self.pre = post
END

;-------------------------------------------------------------------------------------------

PRO supCamData::distortion_correction
;--Distortion/Atmospheric Dispersion Correction (change qfTo_RH* to gqfTo_RH*)
	self->fileDelete, 'distCorrList'
	COMPILE_OPT idl2, HIDDEN
	if (self.pre eq 'NULL') then self.pre = 'c'
	post = 'g'
	self->checkinput, /ndfilt, /ndobj

	print
	print, 'Distortion/Atmospheric Dispersion Correction...'

	rfiles = self->getfiles(nFiles)
	if (nFiles eq 0) then stop

	files = self->supCheckConstraints(rfiles, {FILTER01:self.filter, OBJECT:self.object})
	self->write2file, files, 'distCorrList'

	spawn, 'distcorr.csh distCorrList'															;apply the distortion/dispersion correction

	self->putInDir, self.pre+'*.fits.gz', 'dir'+self.pre

	for i=0,(nFiles - 1) do begin
		ipos = strpos(files[i],'.gz')
		if (ipos eq -1) then outFile = post+files[i] $
			else outFile = post + strtrim(strmid(files[i],0,ipos),2)
		if (iPos eq -1) then stop
		img = mrdfits(post+files[i],0,hdr,/silent)
		bPix = where(img eq -32768, nBPix)
		if (nBPix gt 0) then img[bPix] = !values.f_nan
		sxaddhist, 'supcamdata__define.pro: Distortion/astmospheric dispersion corrected for.', hdr
		mwrfits, img, outFile, hdr
		self->fileDelete, post+files[i]
	endfor

	spawn, 'gzip '+post+'*.fits'
	self->putInDir, post+'*.fits.gz', 'dir'+post
	self.pre = post
	self->fileDelete, 'distCorrList'
	print
END

;-------------------------------------------------------------------------------------------

PRO supCamData::mask_AGprobe
;--Masking the AG Probe (change from gqfTo_RH* to AgqfTo_RH*)
	compile_opt idl2, hidden
	if (self.pre eq 'NULL') then self.pre = 'g'
	post = 'A'
	self->checkinput, /ndfilt

	rfiles = self->getfiles(nFiles)
	if (nFiles eq 0) then stop

	files = self->supCheckConstraints(rfiles, {FILTER01:self.filter})
	;self->write2File, files, 'AGProbe'
	;spawn, 'mask_AGX.csh AGProbe'																;Mask the AG probe

	flag = 0
	for i=0,(nFiles - 1) do begin
		img = mrdfits(files[i],0,hdr,/silent)
		ipos = strpos(files[i],'.gz')
		if (ipos eq -1) then outFile = post+files[i] $
			else outFile = post + strtrim(strmid(files[i],0,ipos),2)
		chip = strtrim(sxpar(hdr, 'DETECTOR'))
		AGX = sxpar(hdr, 'S_AG-X')
  		j_limit = 5000
	    ySize = n_elements(img[0,*])
	    if (self.date lt 2008.58) then begin
	    	if ((chip eq 'si006s') or (chip eq 'si002s') or (chip eq 'w6c1')) then j_limit = AGX*60 - 2300
	    	if ((chip eq 'w67c1') or (chip eq 'w7c3')) then j_limit = floor(AGX*60 - 2000)
	    endif else begin
	    	if ((chip eq 'clarisse') or (chip eq 'fio') or (chip eq 'kiki')) then j_limit = AGX*60 - 2300
	    	if ((chip eq 'chihiro') or (chip eq 'nausicaa')) then j_limit = floor(AGX*60 - 2000)
	    endelse
		if (j_limit le (ySize-1)) then begin
			img[*,j_limit:*] = self.badPixValue
			flag = 1
		endif
		sxaddhist, 'supcamdata__define.pro: AG probe masked.', hdr
		mwrfits, img, outFile, hdr
		statusbar, i, nFiles, msg='Masking the AG Probe ...'
	endfor
	self->putInDir, self.pre+'*.fits.gz', 'dir'+self.pre

	self->gZip, post+'*.fits', msg='GZipping output files ...'
	self->putInDir, post+'*.fits.gz', 'dir'+post
	self.pre = post
	self->fileDelete, post+'*.fits'
END

;-------------------------------------------------------------------------------------------

PRO supCamData::subtract_sky
;--Subtract the sky
	compile_opt idl2, hidden

	if (self.pre eq 'NULL') then self.pre = 'A'
	post = 's'

	self->checkInput, /ndfilt, /ndobj

	;str = 'Determining the Sky ...'
	;print
	;print, str, format='(A-'+strtrim(strlen(str),2)+',$)'

	rfiles = self->getFiles(nFiles)
	if (nFiles eq 0) then stop

	nChips = n_elements(self.chips)
	files = self->supCheckConstraints(rFiles, {FILTER01:self.filter, OBJECT:self.object}, count=nFiles)
	if (nFiles eq 0) then stop

	outFiles = files
	for i=0,(nFiles - 1) do begin
		ipos = strpos(files[i],'.gz')
		if (ipos gt -1) then outFiles[i] = post + strtrim(strmid(files[i],0,ipos),2)
	endfor

	nExp = nFiles / nChips
	skymap = fltarr(5,2,nexp)
	k = 0
	for ixp=0,(nExp - 1) do begin
		for ic=0,(nChips - 1) do begin
			img = mrdfits(files[ixp*nChips+ic],0,/silent)
			hold = img[where(finite(img) and (img gt 0))]
			mmm, hold, mode
			undefine, img
			skymap[ic mod (nChips/2),ic/(nChips/2),ixp] = mode
			statusbar, k, nChips*nExp, msg='Determining the sky ...'
			++k
		endfor
	endfor

	k = 0
	for ixp=0,(nExp - 1) do begin
		for ic=0,(nChips - 1) do begin
			img = mrdfits(files[ixp*nChips+ic],0,hdr,/silent)
			bPix = where( ~finite(img), nBPix)
			img -= (medVal = median(skymap[*,*,ixp],/even))
			if (nbpix gt 0) then img[bPix] = self.badPixValue
			sxaddpar, hdr, 'SKYLEVEL', medVal, '[e-/s]'
			sxaddhist, 'supCamData__define.pro: Sky subtracted, median(mode(chips)) = '+strtrim(medVal,2), hdr
			mwrfits, img, outFiles[ixp*nChips+ic], hdr, /create
			statusbar, k, nChips*nExp, msg='Subtracting the sky ...'
			++k
		endfor
	endfor

	self->putInDir, self.pre+'*.fits*', 'dir'+self.pre
	self->gZip, post+'*.fits', msg='GZipping output files ...'
	self->putInDir, post+'*.fits.gz', 'dir'+post
	self.pre = post
END

;-------------------------------------------------------------------------------------------

PRO supCamData::astrom, redoastrom=redoastrom, catalog=catalog, nthreads=nthreads, nobatch=nobatch, noinherit=noinherit, distcorr=distcorr
;--Calibrate astrometry
	compile_opt idl2, hidden

	;post = 'As'
	;supcam_astrom_x, pre, redoastrom=redoastrom, catalog=catalog, nthreads=nthreads, nobatch=nobatch, noinherit=noinherit

	if (self.pre eq 'NULL') then self.pre = 's'

	files = self->file_check(nFiles) ;create a 'pre' directory or retrieve files from it
	if (nFiles eq 0) then stop

	;files = [files[170:179],files[-10:*]]

	post = 'aS'
	self->calibAstrom, files, post, self.date, redoastrom=redoastrom, catalog=catalog, $
		nthreads=nthreads, nobatch=nobatch, noinherit=noinherit, distcorr=distcorr	;perform astrometry
	file_delete, 'dir'+post, /allow_nonexistent, /recursive
	self->gZip, post+'*.fits', msg='GZipping output files ...'
	;spawn, 'gzip '+post+'*.fits'
	self->putInDir, post+'*.fits.gz', 'dir'+post
	self->putInDir, self.pre+'*.fits.gz', 'dir'+self.pre
	foreach str, ['*.reg','*.cl','dumpFile*','xoutPhot*','xphotcoords*'] do self->putInDir, str, 'dir_JUNK'
	self.pre = post
END

;-------------------------------------------------------------------------------------------

PRO	supCamData::astrom_store, centroid, shift, cats, shifts
	COMPILE_OPT idl2, HIDDEN
	for i=0l,(n_elements(shift)-1l) do begin
		ind = where( (centroid.x gt shift[i].delx) and (centroid.x lt (shift[i].delx + shift[i].xsize)) and (centroid.y gt shift[i].dely) and (centroid.y lt (shift[i].dely + shift[i].ysize)), nobj)
		if (nobj eq 0) then stop
		dum = create_struct('num'+strtrim(i,2), arr_struct({file:replicate(strtrim(shift[i].file,2),nobj), x:centroid[ind].x - shift[i].delx, $
			y:centroid[ind].y - shift[i].dely, ra:centroid[ind].ra, dec:centroid[ind].dec}))
		;dum = create_struct('num'+strtrim(i,2), arr_struct({file:replicate(strtrim(shift[i].file,2),nobj), x:centroid[ind].x - shift[i].delx, $
		;	y:centroid[ind].y - shift[i].dely, xerr:float(centroid[ind].xerr), yerr:float(centroid[ind].yerr), ra:centroid[ind].ra, dec:centroid[ind].dec}))
		cat = struct_addtags(cat, dum)
	endfor
	if (n_elements(cats) gt 0) then ntags = strtrim(n_tags(cats),2) else ntags = '0'
	dcats = create_struct('num'+ntags, cat)
	cats = struct_addtags(cats, dcats)
	dshift = create_struct('num'+ntags, shift)
	shifts = struct_addtags(shifts, dshift)
END

;-------------------------------------------------------------------------------------------

;use nobatch when the astrometry is coming out badly
PRO supCamData::calibAstrom, file, post, date, redoastrom=redo, catalog=catalog, nthreads=nthreads, nobatch=nobatch, noinherit=noinherit, distcorr=distcorr
	COMPILE_OPT idl2, HIDDEN

	nodestroy = 1
	savfile = 'supcamastrom_save.sav'	; save file name
	dummy = file_search(savfile, count=savexist)	; check for existence
	if ((savexist eq 0) and (n_elements(redo) gt 0)) then undefine, redo	; if it doesn't exist and redo is set, undefine redo

	if (n_elements(redo) eq 0) then begin
		dummy = file_search(savfile, count=savexist)

		exptime = fltarr(n_elements(file))
		for i=0,(n_elements(file) - 1) do exptime[i] = sxpar(headfits(file[i]), 'OEXPTIME')

		nind = reverse(sort(exptime[indgen(n_elements(file)/10)*10]))
		newind = 0
		for j=0,(n_elements(nind)-1) do newind = [newind, indgen(10)+nind[j]*10] ; reorder so that the rough astrometry is calibrated to a mosaic with a longer exposure time
		file = file[newind[1:*]]

		mosfile = supalign2(file[0:9], datamax, shift=shift) ; make rough mosaics

		datamax0 = datamax
		astrometry, mosfile, /object, /overwrite, catalog=catalog, datamax=datamax, nthreads=nthreads, bridge=bridge, nodestroy=nodestroy, centroid=centroid

		self->astrom_store, centroid, shift, cats, shifts

		if (n_elements(file) gt 10) then begin
			mosfile = [mosfile, supalign2(file[10:*], datamax)] ; make rough mosaics
			datamax = [datamax0,datamax]
			if (n_elements(nobatch) ne 0) then begin
				for i=1,(n_elements(mosfile) - 1) do begin
					astrometry, mosfile[i], /overwrite, catalog=catalog, datamax=datamax[i], /rough
					stop
				endfor
			endif
		endif

		fits_info, mosfile[0], /silent, n_ext=n_ext
		for hdu=1,n_ext do begin	; cycle through the hdus of the first mosaic to find the astrometry object list
			dummy = mrdfits(mosfile[0], hdu, /silent)
			if (type(dummy, /string) eq 'structure') then begin
				if tag_exist(dummy, 'raj2000') then begin
					astrom = dummy
					break
				endif
			endif else stop
		endfor

		;donate updated astrometry to other mosaics
		hdr0 = headfits(mosfile[0])	; get the header of the first mosaic
		shdr0 = headfits(strtrim((mrdfits(mosfile[0], 1, /silent))[0].file)) ; get the header of its first subimage

		for i=1,(n_elements(mosfile) - 1) do begin	; for each subsequent mosaic
			if (i eq (n_elements(mosfile) - 1)) then if (n_elements(nodestroy) ne 0) then undefine, nodestroy

			if (n_elements(nobatch) eq 0) then begin
				hdr1 = headfits(mosfile[i])	; get the mosaic header
				shdr1 = headfits(strtrim((mrdfits(mosfile[i], 1, /silent))[0].file)) ; get the first subimage header
				xyxy, shdr1, shdr0, 0, 0, xdif, ydif	; calculate offset between exposures
				extast, hdr0, astr		;get the astrometry from the first mosaic
				astr.crval = astr.crval - [hms2dec(sxpar(shdr0, 'RA2000'))*15d, hms2dec(sxpar(shdr0, 'DEC2000'))] + $
					 [hms2dec(sxpar(shdr1, 'RA2000'))*15d, hms2dec(sxpar(shdr1, 'DEC2000'))]	;offset reference coords to the current mosaic
				putast, hdr1, astr	; put this astrometry into the mosaics header
				modfits, mosfile[i], 0, hdr1	; put the header back in the file
			endif

			mwrfits, astrom, mosfile[i], /silent	; write the object list to the file as well
			astrometry, mosfile[i], /object, /batch, /overwrite, catalog=catalog, datamax=datamax[i], nthreads=nthreads, bridge=bridge, nodestroy=nodestroy, centroid=centroid 			;run astrometry in batch mode now that the astrometric solution is close to correct
			shift = mrdfits(mosfile[i],1,/silent)
			self->astrom_store, centroid, shift, cats, shifts
		endfor
		save, mosfile, shifts, cats, filename=savfile
	endif else restore, savfile
	if (n_elements(astr) gt 0) then undefine, astr
	if (n_elements(noinherit) eq 0) then begin
		for i=0,(n_elements(mosfile) - 1) do begin
			;dastr = supcam_mos_astromMC(cats.(i), date, filenames=strtrim(shifts.(i).file,2), maxiter=2, start=start, distcorr=distcorr, post=post)
			;astr = struct_append(astr, dastr)
			;start = dastr

                                ;call write_cats to write out catalog
                                ;thematpl files to .csv
                                ;will also write out list of files to
                                ;perform astrometry on
                   write_cat,cats.(i)

		endfor

                   ;make the list of files to be calibrated in python
                   spawn,'ls -1 sAgczfTo* > astrom_file_list.lst'
				print,"In seperate window, run: python -m full_solver 'astrom_file_list.lst'"
				print,"After full_solver is finished running, .continue or begin from step 11."
                   ;spawn,"python -m full_solver 'astrom_file_list.lst'"
				stop
		;for i=0,(n_elements(mosfile) - 1) do astrom_inherit, mosfile[i], pre=post				;cycle through each mosaic and inherit astrometry into subimages
	endif
END

;-------------------------------------------------------------------------------------------

PRO supCamData::scale_exps
;--Multiplicatively scale each exposure to the median throughput level of all the exposures
;--Compare the flux levels of stars between exposures to determine the relative scaling
	compile_opt idl2, hidden
	if (self.pre eq 'NULL') then self.pre = 'aS'
	post = 'x'
	print, 'supCamData::scale_exps - need to add logic to skip re-calculation of the scaling factors...'
	scaleFits = 'diag_ScaleExposures.fits'
	files = self->file_check(nFiles)
	if (nFiles eq 0) then stop					;create a 'pre' directory or retrieve files from it

	mosFiles = file_search('tmos*fits*',count=nMosFiles)

	if (nMosFiles eq 0) then stop

	if (nMosFiles gt 1) then begin

		dummy = file_search(scaleFits,count=nMatch)

		if (nMatch gt 0) then begin
			redo = 'NO'
			scFits = mrdfits(scaleFits,1,/silent)
			if (n_elements(scFits) eq nMosFiles) then begin					;make sure that there are as many entries as mosaic files
				for k=0,(nMosFiles - 1) do begin					;make sure there is an entry for each mosaic file
					dummy = where(strtrim(scFits[k].fileName) eq mosFiles,nMatch)
					if ~nMatch then begin
						redo = 'YES'						;if an entry is missing, then we must remake scaleFits
						break
					endif
				endfor
			endif else redo = 'YES'								;remake if there is a different number of entries than mosaic files
		endif else redo = 'YES'																	;if scaleFits does not exist, then we must make it
		if (strupcase(redo) eq 'NO') then redo = strupcase(choose('Scaling already determined, recalculate?', 'no')) ;Make sure it should not be remade
		if (strupcase(redo) eq 'YES') then begin
			catList = list()
			for i=0,(nMosFiles - 1) do begin
				objects = runse(mosFiles[i], zp=29., exptime=1d, /batch, /cleanup)
				dCat = seSelection(objects, {magcut:[18d, 22d], CLASS_STAR:0.95d, fwhmcut:[90d, 70d], apercut:[90d, 90d], sncut:[40d, 80d]}, /simple, /batch)	;want to remove the explicit keywords once they are propagated in the header
				if (n_elements(dCat) lt 50) then stop
				catList.add, dCat, i																	;add the catalog to the list
			endfor
			;
			;self->fileDelete, 'use_wt_tmos*fits'
			;self->fileDelete, 'use_tmos*fits'
			;self->fileDelete, 'check.fits'

			ratio = dblarr(nMosFiles-1,nMosFiles) + !values.d_nan
			ratio[0,0] = 1d
			for i=0,(nMosFiles - 2) do begin
				print, mosFiles[i]
				d0 = catList[i]
				if (i gt 0) then ratio[i,i] = 1d
				d0.flux_aper2 *= ratio[i,i]
				for j=(i+1),(nMosFiles - 1) do begin
					d1 = catList[j]
					d1.flux_aper2 *= ratio[i,i]
					match = match_2d(d0.x_world, d0.y_world, d1.x_world, d1.y_world, 5d/3600d, match_distance=mindist)
					ind = where( (match ne -1) and (mindist*3600d lt 4d), ngood)

					if (nGood eq 0) then stop
					d0 = d0[ind]
					d1 = d1[match[ind]]
					print, d1.flux_aper2/d0.flux_aper2
					ratio[i,j] = median(d1.flux_aper2/d0.flux_aper2,/even)								;
				endfor
				ratio[i,*] *= ratio[0,i]	;make all these ratio relative to the first exposure
			endfor
			scale = median(ratio,dim=1,/even)
			scale /= median(scale,/even)
		endif else scale = 1d/scFits.scale

		k = 0
		for i=0,(nMosFiles - 1) do begin

			pattern = (strsplit(strsplit(mosFiles[i],'tmos_',/extract,/regex),'.fits',/extract,/regex))[0]
			files = file_search(self.pre+'*'+pattern+'*fits*',count=nFiles)
			if (nFiles eq 0) then stop
			for j=0,(nFiles - 1) do begin
				img = mrdfits(files[j],0,hdr,/silent)
				;bPix = where( ~finite(img), nBPix)
				img /= scale[i]			;normalize the image to the median flux throughput of all the exposures
				;if (nBPix gt 0) then img[bPix] = !values.f_nan
				outFile = post + strtrim(strmid(files[j],0,(ipos = strpos(files[j],'.gz'))),2)
				sxaddpar, hdr, 'EXPSCALE', 1d/scale[i], 'Multiplicatively scaled to median throughput.'
				hSat = sxpar(hdr,'SATURATE')
				if (hSat eq 0) then stop
				sxaddpar, hdr, 'SATURATE', hSat/scale[i], /saveComment
				gain = sxpar(hdr,'GAIN')
				if (gain eq 0) then stop
				sxaddpar, hdr, 'GAIN', gain*scale[i], /saveComment
				mwrfits, float(img), outFile, hdr, /create
				statusBar, k, nFiles*nMosFiles, msg='Rescaling exposures to common throughput level ...'
				++k
			endfor
		endfor
		self->gZip, post+'*.fits', msg='GZipping output files ...'
		self->putInDir, post+'*.fits.gz', 'dir'+post
		mwrfits, arr_struct({filename:mosfiles, scale:1d/scale}), scaleFits, /create
	endif else begin
		post = self.pre
		print
		print, 'Only one mosaic found.'
		print, 'Skipping multiplicative scaling.'
		print
	endelse
	self->putInDir, self.pre+'*.fits*', 'dir'+self.pre
	self->putInDir, scaleFits, 'dir_DIAG'
	foreach str, ['secat','SExtractor*.sex','SEoutput','SEconfig_*uparm','*nnw'] do self->putInDir, str, 'dir_JUNK'
	self.pre = post
END

;-------------------------------------------------------------------------------------------

PRO supCamData::reproject
;--Reproject Images to a Common Astrometric Solution
	compile_opt idl2, hidden
	if (self.pre eq 'NULL') then self.pre = 'x'
	files = self->file_check(nFiles)														;create a 'pre' directory or retrieve files from it

	for i=0,(nFiles - 1) do begin															;for each file
		self->gunZip, files[i], newFileName
		img = mrdfits(newFileName, 0, hdr, /silent)											;read in the image and header
		bpix = where(img le -32768, nbpix)													;find bad pixels
		if (nbpix gt 0) then begin
			print, 'there should not be any pixels with -32768, only NaNs'
			stop
			img[bpix] = !values.f_nan														;make them NaNs
		endif
		mwrfits, float(img), newFileName, hdr, /create										;overwrite the original
		;spawn, 'mFixNan ' + files[i] + ' ' + files[i] + ' -1d-1000 -32768'					;turn any bad pixels into NaNs
		statusbar, i, nFiles, msg='Converting bad pixels to NaNs...'
	endfor
	files = self->file_check(nFiles)														;create a 'pre' directory or retrieve files from it
	if (nFiles eq 0) then stop

	file_mkdir, 'scratch_proj'																;make a scratch directory
	file_move, files, 'scratch_proj/'+files, /overwrite										;move the relevant files to this directory

	file_delete, 'projdir', /allow_nonexistent, /recursive									;get rid of any existing projdir files
	file_mkdir, 'projdir'																	;make a directory to hold the projected images

	spawn, 'mImgtbl scratch_proj images-rawdir.tbl'											;Generate image metadata table
	spawn, 'mMakeHdr -n images-rawdir.tbl template.hdr'										;Make template header

	nthreads = 1

	itime = systime(1)
	if (nthreads gt 1) then begin
		nrows = file_lines('images-rawdir.tbl')												;count the rows in the overlap.tbl
		openr, lun, 'images-rawdir.tbl', /get_lun											;open it
		d = replicate({val:' '}, nrows)														;create a structure to hold the data
		readf, lun, d																		;read the table into the structure
		free_lun, lun																		;close the file

		spawn, 'pwd', basedir
		for i=0,(nthreads - 1) do begin
			if (i eq 0) then bridge = create_struct('n'+strtrim(i,1), Obj_New('IDL_IDLBridge', output='')) $
						else bridge = struct_addtags(bridge, create_struct('n'+strtrim(i,1), Obj_New('IDL_IDLBridge', output='')))	;create the child processes
			bridge.(i) -> SetVar, 'basedir', basedir
			bridge.(i) -> Execute, "cd, '" + basedir
		endfor

		status = intarr(nthreads)															;array carrying status of each child process
		pending = intarr(nthreads)															;array carrying status of readout state

		itime = systime(1)
		idle = 0
		ncomp = 0
		noverlaps = n_elements(d) - 2
		itime = systime(1)
		for i=0,(noverlaps - 1) do begin
			openw, lun, 'images-rawdir'+strtrim(i,2)+'.tbl'									;open a new file
			foreach j, [0,1,i+2] do printf, lun, d[j].val									;write the header and a single entry of the table to this file
			free_lun, lun

			bridge.(idle) -> Execute, '.full'
			bridge.(idle) -> SetVar, 'basedir', basedir
			bridge.(idle) -> Execute, "cd, '" + basedir
			bridge.(idle) -> Setvar, 'i', strtrim(i,2)
			bridge.(idle) -> Execute, "spawn, 'mProjExec -p scratch_proj images-rawdir' + i + '.tbl template.hdr projdir stats' + i + '.tbl", /nowait		;initiate cross-correlation

			status[idle] = 1																;mark the status as busy
			pending[idle] = 1																;specify information is to be read out

			;itime = systime(1)
			while 1 do begin
				while 1 do begin
					for ind=0,(nthreads - 1) do status[ind] = bridge.(ind)->status()		;record the status of the child process
					idle = (where((status eq 0) or (status eq 2), count))[0]				;idle marks an idle child
					if (count gt 0) then break	;wait, 0.1 else break								;if no children are present then wait for a process to complete, otherwise exit
				endwhile
				ind = where(((status eq 0) or (status eq 2)) and (pending eq 1), count)		;ind marks an idle child ready to be read out
				for ready=0,(count - 1) do begin											;for each one of the above
					status[ind[ready]] = 0													;reset status and pending state
					pending[ind[ready]] = 0
					ii = bridge.(ind[ready])->GetVar('i')
					file_delete, 'images-rawdir' + ii + '.tbl'
					file_delete, 'stats' + ii + '.tbl'
					++ncomp																	;number of processes completed
					;if (ncomp eq 100) then PRINT, 'Time Elapsed: ' + strtrim(string(systime(1) - itime, format='(F9.2)'),1) + ' seconds'
					;if (ncomp eq 20) then PRINT, 'Time Elapsed: ' + strtrim(string(systime(1) - itime, format='(F9.2)'),1) + ' seconds'
				endfor
				if ((i lt (noverlaps - 1)) or (total(pending) eq 0)) then break				;if not the last iteration then continue, if last iteration then make sure everything is read out
			endwhile
			;PRINT, 'Time Elapsed: ' + strtrim(string(systime(1) - itime, format='(F9.2)'),1) + ' seconds'
		endfor
		for i=0,(nthreads - 1) do bridge.(i)->Cleanup										;destroy each child process
		undefine, bridge																	;get rid of the bridge structure
            endif else spawn, 'mProjExec -p scratch_proj images-rawdir.tbl template.hdr projdir stats.tbl' ;reproject each image
;       	endif      else spawn, 'mProjExec -p scratch_proj images-rawdir.tbl template.hdr projdir stats.tbl' ;reproject each image
;        print,'Due to montage path issues, run: '
;        print,'mProjExec -p scratch_proj images-rawdir.tbl template.hdr projdir stats.tbl'
;        print,'in a seperate window. This will take a while. Then, when finished, type:'
;        print,'> .c'
;        print,'in the IDL window. You can close the seperate window.'
;	stop
	spawn, 'mImgtbl -c projdir projected.tbl'												;Generate new image metadata table

	file_move, 'scratch_proj/'+files, 'dir'+self.pre+'/'									;move these files to their folder

	self->fileDelete, 'scratch_proj'

	self->mkDir, 'projdir/areaOrig'															;make an area directory inside of corrdir
	afiles = file_search('projdir/hdu*area.fits', count=count)								;find the area files
	for i=0,(count - 1) do begin
		dummy = file_search('projdir/areaOrig/'+aFiles[i]+'.gz',count=nMatch)
		if (nMatch gt 0) then begin
			print, 'A file with this name already eximImgtblsts in the projdir/areaOrig folder.'
			print, 'You might be overwriting the original area files.'
			stop
                    endif

 		if sxpar(headfits(afiles[i]),'BITPIX') eq -64 then spawn, 'mConvert -b -32 ' + afiles[i] + ' ' + afiles[i]	;convert the file back to float
		self->gZip, aFiles[i]
		;file_move, afiles[i]+'.gz', 'projdir/areaOrig/'+(strsplit(afiles[i], '/', /extract, /regex))[-1]	;if any, then move them to the area directory
		self->putInDir, aFiles[i]+'.gz', 'projdir/areaOrig'
		;statusbar, i, count, msg='Converting area images to float and GZip-ing ...'
	endfor
	files = file_search('projdir/hdu*fits', count=nfiles)									;find all the bg rectified files
	for i=0,(nfiles - 1) do begin
		if (sxpar(headfits(files[i]),'BITPIX') eq -64) then spawn, 'mConvert -b -32 ' + files[i] + ' ' + files[i]		;convert the file back to float
		;statusbar, i, nFiles, msg='Converting projected images to float and GZip-ing ...'
	endfor
END

;-------------------------------------------------------------------------------------------

PRO supCamData::weight_area
;--Determine Weighting of each Exposure
	compile_opt idl2, hidden
	self->checkInput, /ndfilt

	;get information about the images and save this to the diag_coadd data structure to be used when co-adding
	if (self.pre eq 'NULL') then self.pre = 'x'
	files = self->getFiles(nFiles)												;move
	if (nFiles eq 0) then stop
	skip = ['LONGPOLE', 'CUNIT1', 'CUNIT2', 'WCS-ORIG', 'RADECSYS', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'LONPOLE', 'LATPOLE', 'PV2_1', 'PV2_2']
	saturate = fltarr(nFiles)													;array to hold the saturation values for each file
	gain = fltarr(nFiles)														;array to hold the gain values for each file
	for i=0,(nfiles - 1) do gain[i] = sxpar(headfits(files[i]),'GAIN')			;grab the gain value from the file header
	for i=0,(nfiles - 1) do saturate[i] = sxpar(headfits(files[i]),'SATURATE')	;grab the saturation level from the file header
	nexp = nFiles/n_elements(self.chips)										;the number of exposures is the number of images divided by the number of CCDs
	diag_coadd = {COADDTYP:self.coaddType, UNITS: 1., EXPTIME:1., NUMEXP:nexp, $
		SATURATE:min(saturate), iSATURATE:saturate, iGAIN:gain, filter:self.filter}

	; This line only works for post-April 2008 data. Need to switch.
	;print, 'The following line  needs to be hardcoded for pre- and post- April 2008 datg
	;stop
	;pfiles = file_search('projdir/hdu*si001s.fits', count=nPFiles)				;get the names of the reprojected images
	;pfiles = file_search('projdir/hdu*chihiro.fits', count=nPFiles)				;get the names of the reprojected images

	if (self.date lt 2008.58) then begin
		pfiles = file_search('projdir/hdu*si001s.fits', count=nPFiles)
	endif else begin
	 	pfiles = file_search('projdir/hdu*chihiro.fits', count=nPFiles)				;get the names of the reprojected images
	endelse

	if (nPFiles eq 0) then stop

	pre = strmid((strsplit(pfiles[0],'hdu0_',/extract,/regex))[-1],0,1)			;make a string to search for pre-reprojected images
	dummy = self->getFiles(count, pre=pre)										;retrieve these files, put them in the working directory
	if (count eq 0) then stop

	wtFiles = strarr(npFiles)
	for i=0,(npFiles - 1) do begin
		str = (strsplit(pfiles[i],'hdu0_',/extract,/regex))[-1]					;trim the prefix from the file name
		iPos = strpos(str,'_',/reverse_search)									;search for the underscore in the image name
		if (ipos[0] eq -1) then stop
		wtFileName = strmid(str,0,ipos[-1])
		wtFile = file_search(wtFileName+'*.fits*',count=nWtFile)				;find the files that have the wtFileName prefix
		if (nWtFile eq 0) then stop
		wtFile = wtFile[0]														;just need one of these image names, otherwise headfits will fail
		hdr = headfits(wtFile)													;get the header
		skyLevel = float(sxpar(hdr, 'SKYLEVEL'))								;extract the count rate of the sky background
		if (skyLevel eq 0) then stop
		expTime = float(sxpar(hdr,'OEXPTIME'))									;get the original exposure time, before rescaling
		if (expTime eq 0) then stop
		expScale = float(sxpar(hdr,'EXPSCALE'))									;get the multiplicative scaling that this exposure needs to match up to the other exposures
		if (expScale eq 0) then stop
		wt = expTime / (skyLevel*expScale^2d)									;the weighting factor for this exposure
		wtInfo = struct_append(wtInfo, {pre:wtFileName, wt:wt})					;save this information to a structure
	endfor

	case strupcase(self.coaddType) of
		'MEDIAN':	finalGain = 2/3.*nExp*mean(gain)							;assume only 2/3 of the data are being used when median combining
		'MEAN'	:	finalGain = 1.*nExp*wtd_mean(gain, reform(rebin(wtInfo.wt/total(wtInfo.wt),nExp,nFiles/nExp),nFiles))
		else	:	stop
	endcase
	diag_coadd = struct_addtags(diag_coadd, {GAIN:finalGain, weight:wtInfo.wt/total(wtInfo.wt)})
	mwrFits, diag_coadd, 'diag_coadd.fits', /create

	aFiles = file_search('projdir/areaOrig/hdu*fits*', count=nAFiles)			;find all the bg rectified files
	if (nAFiles eq 0) then stop
	self->fileDelete, 'hdu*area.fits'											;get rid of any already modified area files
	self->fileDelete, 'projdir/hdu*area.fits'									;get rid of any already modified area files
	self->mkDir, 'projdir/area'													;make an area directory inside of projdir
	self->fileDelete, 'projdir/area/hdu*area.fits'								;make sure there are no already modified area files

	wt = dblarr(nAFiles)														;create a weight array that will contain the weighting for each file
	for i=0,(nAFiles - 1) do begin
		str = (strsplit(aFiles[i],'hdu0_',/extract,/regex))[-1]					;add the projected prefix to each area file name
		iPos = strpos(str,'_',/reverse_search)									;find the last underscore in the file name
		if (ipos[0] eq -1) then stop
		str = strmid(str,0,ipos[-1])											;trim off the underscore and everything after from the file name
		preFileName = str + '.fits'												;add a .fits extension to the file name string, and store to variable
		iPos = strpos(str,'_',/reverse_search)									;find the first underscore in the file name string
		if (ipos[0] eq -1) then stop
		str = strmid(str,0,ipos[-1])											;trim off the previously added prefix

		iwt = where(str eq wtinfo.pre,iMatch)									;grab the coadd weighting from the wtInfo structure
		if (iMatch eq 0) then stop
		wt[i] = wtinfo[iWt].wt/total(wtInfo.wt)									;normalize the weights to 1

		aimg = mrdfits(aFiles[i],0,ahdr,/silent)								;read in the area image array
		imgName = 'projdir/'+(strsplit((strsplit(afiles[i], '/', /extract, /regex))[-1],'_area.fits',/extract,/regex))[0] + '.fits' ;construct the science image file name from the area file name
		dummy = file_search(imgName,count=nMatch)								;make sure the science image exists
		if (nMatch eq 0) then stop
		img = mrdfits(imgName,0,hdr,/silent)									;read in the science image array
		dummy = file_search(preFileName,count=nMatch)							;read in the corresponding image before projection
		if (nMatch eq 0) then stop
		saturate = sxpar(headfits(preFileName),'SATURATE')						;extract the saturation level from the header
		if (saturate eq 0) then stop											;error check
		if (strupcase(self.coaddType) eq 'MEAN') then begin						;if the co-addition type is set to mean, this means weighted mean
			aimg *= wt[i]														;multiply the area image by the image weighting
			msg = 'Weighting images by exposure time and sky level ...'			;
		endif else msg = 'Masking out saturated pixels ...'
		bPix = where(img gt saturate*self.nonLinFrac, nBPix)					;find pixels with values above the non-linearity level
		if (nBPix gt 0) then aImg = maskgrow(aImg, bPix, self.satMaskGrow, 0.)	;mask out those pixels and surrounding pixels with a grow size set by self.satMaskGrow
		mwrfits, float(aImg), self->rmGZ('projdir/area/'+(strsplit(afiles[i], '/', /extract, /regex))[-1]), ahdr, /create	;write the area file in the appropriate directory
		statusbar, i, nAFiles, msg=msg
	endfor
	self->putInDir, self.pre+'*.fits', 'dir'+self.pre							;put the pre-reprojected images back in their directory
END

;-------------------------------------------------------------------------------------------

PRO supCamData::coAdd
	compile_opt idl2, hidden
	self->checkinput, /ndout
	;if (self.pre eq 'NULL') then self.pre = 'x'
	;self->putInDir, self.pre+'*fits*', 'dir'+self.pre

	coaddFile = 'diag_coadd.fits'
	dummy = file_search(coaddFile,count=nMatch)
	if (nMatch eq 0) then stop
	coaddInfo = mrdfits(coaddFile,1,/silent)

	;files = file_search('corrdir/hdu*fits', count=count)
	files = file_search('projdir/hdu*fits', count=nFiles)
	if (nFiles eq 0) then stop
	badPixVal = -32768
	for i=0,(nFiles-1) do begin
		img = mrdfits(files[i],0,hdr,/silent)
		bPix = where(~finite(img),nBad)
		if (nBad gt 0) then begin
			img[bPix] = badPixVal
			mwrfits, img, files[i], hdr, /create
		endif
		statusbar, i, nFiles, msg='Checking for NaNs ...'
	endfor
	print
	areaFiles = file_search('projdir/area/hdu*fits', count=nAreaFiles)
	if (nAreaFiles gt 0) then file_move, areaFiles, 'projdir/'

	file_delete, 'final', /allow_nonexistent, /recursive					;get rid of any existing final files
	file_mkdir, 'final'

	spawn, 'mAdd -p projdir -a ' + strlowcase(coaddInfo.coaddTyp) + ' projected.tbl template.hdr final/mosaic.fits'

	foreach suf, ['']+'.fits' do begin			;'_area'
		img = float(mrdfits('final/mosaic' + suf, 0, hdr, /silent))
		bpix = where(~finite(img), nbpix, ncomplement=nGPix)
		if (nGPix eq 0) then stop
		;if (nbpix gt 0) then img[bpix] = 0.
		doutfile = (strsplit(self.outfile, '.fits', /extract, /regex))[0] + suf
		dummy = file_search(doutfile,count=nMatch)
		if (nMatch gt 0) then begin
			self->fileDelete, doutfile+'.old'
			spawn, 'mv ' + doutfile + ' ' + doutfile + '.old'
		endif

		sxaddpar, hdr, 'FILTER', coaddInfo.filter, 'Imaging filter.'
		sxaddpar, hdr, 'COADDTYP', coaddInfo.coaddtyp, 'Type of co-addition performed.'
		sxaddpar, hdr, 'UNITS', coaddInfo.units, '[e-/s]'
		sxaddpar, hdr, 'EXPTIME', coaddInfo.exptime, '[s]'
		sxaddpar, hdr, 'NUMEXP', coaddInfo.numexp, 'Number of contributing exposures.'
		case strupcase(coaddInfo.coaddtyp) of
			'MEDIAN':	coaddComment = 'Gn*<exptime>*Nexp*2/3'
			'MEAN'	:	coaddComment = 'Gn*<exptime>*Nexp'
			else	:	stop
		endcase
		sxaddpar, hdr, 'GAIN', coaddInfo.gain, 'Effective Gain: '+coaddComment
		sxaddpar, hdr, 'SATURATE', coaddInfo.saturate, 'Saturation/non-linearity level'
		mwrfits, img, doutfile, hdr, /create
	endforeach
	self->putInDir, 'projdir/hdu*area.fits', 'projdir/area'
	print
END

;-------------------------------------------------------------------------------------------

PRO supCamData__define
;--Class Definition
	compile_opt idl2, hidden

	void = {supCamData, chips:strarr(10), date:0., filter:'', object:'', $
		flatType:'', outFile:'', pre:'', batch:0, badPixValue:0., $
		satLevel:0., nonLinFrac:0., satMaskGrow:[0,0], crGrow:0, $
		coAddType:''}

END











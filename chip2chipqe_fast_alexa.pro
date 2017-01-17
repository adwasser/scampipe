FUNCTION scoreit, plot=plot
	;compile_opt idl2
	common share_c2cqe, imgs, xsize, ysize, wid, modx, i, scores, xxx, clafio, kikfio, sansat, shesat, chicla, ponsan, naukik, sopshe, satfio, sancla, shekik, ponchi, sopnau

	claFio;, plot=plot	
	kikFio;, plot=plot
	sanSat;, plot=plot
	sheSat;, plot=plot
	chiCla;, plot=plot
	ponSan;, plot=plot
	nauKik;, plot=plot
	sopShe;, plot=plot
	satFio;, plot=plot
	sanCla;, plot=plot
	sheKik;, plot=plot
	ponChi;, plot=plot
	sopNau;, plot=plot

	dev = abs(scores - median(scores))
	out = where( dev gt 3.*robust_sigma(scores) , nout)			;find outliers, for instance, interfaces where a large galaxy sits
	if (nout gt 0) then scores[out[where(dev[out] eq max(dev[out]))]] = 0.		;only get rid of the biggest outlier
	
	scores[6] = 0.		;?
	return, total(abs(scores))
	
END

;-------------------------------------------------------------------------------------------

FUNCTION newtIFy, X, plot=plot
	plot = 1
	common share_c2cqe
	modx = X	
	ret = scoreit(plot=plot)	
	boundDist = abs(median(modx) - 1.)
	IF (boundDist gt 0.2) THEN ret = ret + boundDist 
	++xxx
	RETURN, ret
END

;-------------------------------------------------------------------------------------------

FUNCTION fluxScale, img0, img1, scale
	gPix0 = where(finite(img0), countGPix0)
	gPix1 = where(finite(img1), countGpix1)
	img0 = TEMPORARY(img0) * scale[0]	
	img1 = TEMPORARY(img1) * scale[1]		
	zero = biweight_mean(img0[gPix0])
	one = biweight_mean(img1[gPix1])
	RETURN, 1d - (one / zero)
END

;-------------------------------------------------------------------------------------------

FUNCTION subImgLR, img0, img1, wid, backwards=backwards
	
	if (n_elements(backwards) ne 0) then begin
		dum = img0
		img0 = img1
		img1 = dum
	endif

	nRows = n_elements(img1[0,*])
	dum1 = dblarr(wid+1,nRows)
	for z=0,(nRows-1) do begin																	;for each row
		ind = where(finite(img1[*,z]), count)													;find which pixels are finite
		if (count gt (wid + 1)) then dum1[0,z] = img1[(reverse(ind))[0:wid],z] $				;if there are enough, then store them
			else dum1[0,z] = dblarr(wid+1)+!values.f_nan										;otherwise, store NaNs
	endfor
	
	nRows = n_elements(img0[0,*])
	dum0 = dblarr(wid+1,nRows)
	for z=0,(nRows-1) do begin	
		ind = where(finite(img0[*,z]), count)
		if (count gt (wid + 1)) then dum0[0,z] = img0[ind[0:wid],z] $
			else dum0[0,z] = dblarr(wid+1)+!values.f_nan
	endfor

	len = min([n_elements(dum0[0,*]),n_elements(dum1[0,*])])-1

	if (n_elements(backwards) eq 0) then return, {reg0:dum0[*,0:len], reg1:dum1[*,0:len]} $
									else return, {reg0:dum1[*,0:len], reg1:dum0[*,0:len]}
END

;-------------------------------------------------------------------------------------------

FUNCTION fluxScaleLR, data0, scale, plot=plot

	dum0 = data0.(0)
	gPix = where(finite(dum0))
	dum0[gpix] = dum0[gpix]*scale[0]

	dum1 = data0.(1)
	gPix = where(finite(dum1))
	dum1[gpix] = dum1[gpix]*scale[1]

	dum = median(dum1,dimension=1)/median(dum0,dimension=1)
	dum = dum[where(finite(dum))]

	for xyz=0,10 do dum = dum[where( abs(dum - median(dum)) le 2.5*stddev(dum - median(dum)) )]
	
	mmm, dum, res

	if (res gt 0) then begin
		dum = median(dum0,dimension=1)/median(dum1,dimension=1)
		dum = dum[where(finite(dum))]
		for xyz=0,10 do dum = dum[where( abs(dum - median(dum)) le 2.5*stddev(dum - median(dum)) )]
		mmm, dum, res
		res = 1d/res
	endif
	
	if (n_elements(plot) gt 0) then begin
		plot, median(dum0,dimension=1), yrange=[-15,15]+mean([median(dum0),median(dum1)]), title=strtrim(1d - res,2)
		oplot, median(dum1,dimension=1), color=fsc_color('red')
		pause
	endif
	RETURN, 1d - res
END

;-------------------------------------------------------------------------------------------

FUNCTION subImgTB, img0, img1, wid


	nCols = n_elements(img1[*,0])
	dum1 = dblarr(wid+1,nCols)
	for z=0,(nCols-1) do begin
		ind = where(finite(img1[z,*]),count)
		if (count gt (wid + 1)) then dum1[0,z] = transpose(img1[z,(reverse(ind))[0:wid]]) $
			else dum1[0,z] = dblarr(wid+1)+!values.f_nan
	endfor

	nCols = n_elements(img0[*,0])
	dum0 = dblarr(wid+1,nCols)
	for z=0,(nCols-1) do begin
		ind = where(finite(img0[z,*]),count)
		if (count gt (wid + 1)) then dum0[0,z] = transpose(img0[z,ind[0:wid]]) $
			else dum0[0,z] = dblarr(wid+1)+!values.f_nan
	endfor

	len = min([n_elements(dum0[0,*]),n_elements(dum1[0,*])])-1

	return, {reg0:dum0[*,0:len], reg1:dum1[*,0:len]}
END

;-------------------------------------------------------------------------------------------

FUNCTION fluxScaletb, data0, scale, plot=plot

	dum0 = data0.(0)
	;gpix = where(dum0 ne -32768)
	gPix = where(finite(dum0))
	dum0[gpix] = dum0[gpix]*scale[0]

	dum1 = data0.(1)
	;gpix = where(dum1 ne -32768)
	gPix = where(finite(dum1))
	dum1[gpix] = dum1[gpix]*scale[1]

	dum = median(dum1,dimension=1)/median(dum0,dimension=1)
	;dum = dum[0:min([n_elements(dum0[0,*]),n_elements(dum1[0,*])])-1]
	;dum = dum[where( (dum gt 0) and (dum lt 2.) )]
		;dum = dum[where( (dum gt -32768)  )]
		dum = dum[where(finite(dum))]

	for xyz=0,10 do dum = dum[where( abs(dum - median(dum)) le 2.5*stddev(dum - median(dum)) )]

	mmm, dum, res

	if (res gt 0) then begin
		dum = median(dum0,dimension=1)/median(dum1,dimension=1)
		;dum = dum[0:min([n_elements(dum0[0,*]),n_elements(dum1[0,*])])-1]
		;dum = dum[where( (dum gt 0) and (dum lt 2.) )]
			;dum = dum[where( (dum gt -32768)  )]
			dum = dum[where(finite(dum))]
		
		for xyz=0,10 do dum = dum[where( abs(dum - median(dum)) le 2.5*stddev(dum - median(dum)) )]
		mmm, dum, res
		res = 1d/res
	endif
	
	if (n_elements(plot) gt 0) then begin
		plot, median(dum0,dimension=1), yrange=[-15,15]+mean([median(dum0),median(dum1)]), title=strtrim(1d - res,2)
		oplot, median(dum1,dimension=1), color=fsc_color('red')
		pause
	endif
	
	RETURN, 1d - res
END

;-------------------------------------------------------------------------------------------

PRO chiCla, plot=plot
	common share_c2cqe
	k = 1
	l = k - 1
	m = 0
	scale = [modx[k], modx[l]]
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		chicla = subImgLR(img0, img1, wid)
	endif 
	scores[m] = fluxScaleLR(chicla, scale, plot=plot)
END

;-------------------------------------------------------------------------------------------

PRO claFio, plot=plot
	common share_c2cqe
	k = 2
	l = k - 1
	m = 1
	scale = [modx[k], modx[l]]
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		claFio = subImgLR(img0, img1, wid)
	endif
	scores[m] = fluxscalelr(claFio, scale, plot=plot)
END

;-------------------------------------------------------------------------------------------

PRO kikFio, plot=plot
	common share_c2cqe
	k = 2
	l = k + 1
	m = 2
	scale = [modx[k], modx[l]]
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		kikFio = subimglr(img0, img1, wid, /backwards)
	endif
	scores[m] = fluxscalelr(kikFio, scale, plot=plot)
END

;-------------------------------------------------------------------------------------------

PRO nauKik, plot=plot
	common share_c2cqe
	k = 3
	l = k + 1
	m = 3
	scale = [modx[k], modx[l]]
	
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		nauKik = subimglr(img0, img1, wid, /backwards)
	endif 

	scores[m] = fluxscalelr(nauKik, scale, plot=plot)


END

;-------------------------------------------------------------------------------------------

PRO ponSan, plot=plot
	common share_c2cqe
	k = 6
	l = k - 1
	m = 4
	scale = [modx[k], modx[l]]
	
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		ponSan = subimglr(img0, img1, wid)
	endif 

	scores[m] = fluxscalelr(ponSan, scale, plot=plot)

END

;-------------------------------------------------------------------------------------------

PRO sanSat, plot=plot
	common share_c2cqe
	k = 7
	l = k - 1
	m = 5
	scale = [modx[k], modx[l]]
	
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		sanSat = subimglr(img0, img1, wid)
	endif 

	scores[m] = fluxscalelr(sanSat, scale, plot=plot)

END

;-------------------------------------------------------------------------------------------

PRO satFio, plot=plot
	common share_c2cqe
	k = 2
	l = k + 5
	m = 6
	scale = [modx[k], modx[l]]
	
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		satFio = subimgtb(img0, img1, wid)
	endif 

	scores[m] = fluxscaletb(satFio, scale, plot=plot)
	
END	

;-------------------------------------------------------------------------------------------

PRO sheSat, plot=plot
	common share_c2cqe
	k = 7
	l = k + 1	
	m = 7
	scale = [modx[k], modx[l]]
	
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		sheSat = subimglr(img0, img1, wid, /backwards)
	endif 

	scores[m] = fluxscalelr(sheSat, scale, plot=plot)
	
END

;-------------------------------------------------------------------------------------------

PRO sopShe, plot=plot
	common share_c2cqe
	k = 8
	l = k + 1
	m = 8
	scale = [modx[k], modx[l]]
	
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		sopShe = subimglr(img0, img1, wid, /backwards)
	endif 

	scores[m] = fluxscalelr(sopShe, scale, plot=plot)
	
END

;-------------------------------------------------------------------------------------------

PRO ponChi, plot=plot
	common share_c2cqe
	k = 0
	l = k + 5
	m = 9
	scale = [modx[k], modx[l]]
	
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		ponChi = subimgtb(img0, img1, wid)
	endif 

	scores[m] = fluxscaletb(ponChi, scale, plot=plot)
	
END	

;-------------------------------------------------------------------------------------------

PRO sanCla, plot=plot
	common share_c2cqe
	k = 1
	l = k + 5	
	m = 10
	scale = [modx[k], modx[l]]
	
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		sanCla = subimgtb(img0, img1, wid)
	endif 

	scores[m] = fluxscaletb(sanCla, scale, plot=plot)
	
END	

;-------------------------------------------------------------------------------------------

PRO sheKik, plot=plot
	common share_c2cqe
	k = 3
	l = k + 5		
	m = 11
	scale = [modx[k], modx[l]]
	
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		sheKik = subimgtb(img0, img1, wid)
	endif 

	scores[m] = fluxscaletb(sheKik, scale, plot=plot)
	
END	

;-------------------------------------------------------------------------------------------

PRO sopNau, plot=plot
	common share_c2cqe
	k = 4
	l = k + 5
	m = 12
	scale = [modx[k], modx[l]]
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		sopNau = subimgtb(img0, img1, wid)
	endif 
	scores[m] = fluxscaletb(sopNau, scale, plot=plot)
END	

;---------------------------------------------------------------------------------------------------

PRO chip2chipqe_fast, files, pre	
	common share_c2cqe	;common block
	
	starttime0 = SYSTIME(1)
	
	numFiles 	= N_ELEMENTS(files)											
	numChips 	= 10	
	numExps 	= numFiles/numChips		
	
	file = strarr(numChips, numExps)
	print, "stopiing to switchover chip names"
	stop

;	chip = ['chihiro', 'clarisse', 'fio', 'kiki', 'nausicaa', $
;		'ponyo', 'san', 'satsuki', 'sheeta', 'sophie'] ;names of each detector
	; I don't know if the order matters here at all
	chip = ['si001s', 'si002s', 'si005s', 'si006s', 'w4c5', $
		'w67c1', 'w6c1', 'w7c3', 'w93c2', 'w9c2'] ;names of each detector
	
	xsize = INTARR(numChips, numExps)
	ysize = INTARR(numChips, numExps)
	
	FOR k=0,(numChips - 1) DO BEGIN	;cycle through each chip
		j = 0
		FOR i=0,(numFiles - 1) DO BEGIN ;cycle through each image
			IF (N_ELEMENTS(STRSPLIT(files[i], chip[k], /regex)) gt 1) THEN begin	;see if it is on this chip
				file[k,j] = files[i]						;record the filename
				hdr = headfits(file[k,j])					;retrieve the header
				xsize[k,j] = sxpar(hdr, 'NAXIS1')				;record the width
				ysize[k,j] = sxpar(hdr, 'NAXIS2')				;record the height
				j++
			ENDIF
		ENDFOR
	ENDFOR
	
	scores = DBLARR(13)	;create a scores array to hold the scores for the minimization routine
	modx = DBLARR(10) + 1	;
	
	IF (N_ELEMENTS(wid) EQ 0) THEN wid = 50	;set the width of the region around the edges of the chips to analyze
	
	addit = FLTARR(numChips, numExps) + 1.	;
	date = fix(strtrim(strsplit(sxpar(headfits(file[0]), 'DATE-OBS'),'-',/extract))) ;extract the date information
	date = date[0] + date[1]/12. + date[2]/365.					;turn this into a decimal
	
		
	FOR i=0,(numExps - 1) DO BEGIN		;cycle through each exposure
		if (n_elements(imgs) gt 0) then undefine, imgs
		if (n_elements(hdrs) gt 0) then undefine, hdrs		
		FOR k=0,(numChips - 1) DO BEGIN	;cycle through each chip
			img = CREATE_STRUCT('num' + STRTRIM(k,1), float(mrdfits(file[k,i], 0, hdr, /silent)))	;read in the image to a structure
			imgs = struct_addtags(imgs, img)	;store it to a structure
			dhdr = CREATE_STRUCT('num' + STRTRIM(k,1), hdr)	;store the header to a structure
			hdrs = struct_addtags(hdrs, dhdr)	;store this to a structure
			xsize[k,i] = N_ELEMENTS(imgs.(k)[*,0])	;store the actual width
			ysize[k,i] = N_ELEMENTS(imgs.(k)[0,*])	;store the actual height
		ENDFOR
		
		IF (i gt 0) THEN BEGIN		;if this is beyond the first exposure
			x = addit[*,i-1]	;set addit to the results of the previous exposure
		ENDIF ELSE BEGIN											
			if (date le 2009.) then begin	;for older data with lingering non-linearity issues, allow a different scaling for different exposure times
				x = float([1.00000, 0.975279, 0.940946, 0.975610, 0.954359, 1.06361, 0.936070, 0.940543, 0.963282, 1.00205])
			endif else x = float([1.00000, 0.962169, 0.933242, 0.970565, 1.01779, 1.09830, 0.960591, 0.939223, 0.979983, 1.03439])	;initial guesses at the scale, use last images values IF possible
		ENDELSE

		xxx = 0
		;test = newtify(x)
		;stop
		;addit[*,i] = NEWTON(x, 'newtify', /double)		;NEWTON-raphson solver
		
		
		addit[*,i] = tnmin('newtify', x, bestmin=bestmin, /autoderivative, /quiet)					;determine the relative scaling between chips
		addit[*,i] = addit[*,i]/addit[0,i]															;normalize by the first chip
		
		
		
		window, 0
		device, decomposed=0
		color = rebin(findgen(220/numExps*numExps),numExps) + (255 - 220/numExps*numExps)/2
		loadct, 39, /silent
		plot, [0], [0], xrange=[-1,10], /xstyle, xticklen=1d-9, xtickinterval=1, yrange=minmax(addit)+[-1,1]*0.01, /ystyle, $
			xtickname=['chi','cla','fio','kiki','nau','pon','san','sat','she','sop']
		for j=0,(n_elements(addit[0,*])-1) do oplot, findgen(10), addit[*,j], psym=10, color=color[j], thick=2
		loadct, 0, /silent		
		statusbar, i, numExps, msg='Finding relative scaling between chips ...'
	ENDFOR
		
	if (date le 2009.) then begin																	;for older data with lingering non-linearity issues, allow a different scaling for different exposure times
		finaddit = addit
		etime = fltarr(numexps)
		for i=0,(numexps - 1) do etime[i] = sxpar(headfits(file[0,i]), 'EXPTIME')
		etimes = etime[sort(etime)]
		etimes = etimes[uniq(etimes)]
		for j=0,(n_elements(etimes)-1) do begin 			;for each unique exposure time do
			ind = where(etime eq etimes[j], count)
			if (count gt 1) then begin
				for k=0,(numchips - 1) do finaddit[k,[ind]] = median(addit[k,[ind]], /even)
				;for k=0,(numchips - 1) do finaddit[k,[ind]] = wtd_mean(addit[k,[ind]], etime[ind])
			endif																					;otherwise just leave finaddit column as is
		endfor	
	endif else finaddit = rebin(median(addit, dimension=2, /even), [numchips, numexps])				;for newer data just take the median
	
	
	
	cleanplot, /silent
	erase
	psopen, 'diag_chipscaling', /encapsulated, /color, /heiles, xs=8, ys=8, /inches
	!p.multi=[0,1,2]
	!p.charsize = 1.5
	!p.charthick = 3
	yrange = float(roundx(1 + [-1,1]*max(abs(finaddit - 1))*1.5,3))
	plot, indgen(10), finaddit[*,0], xrange=[-1,5], xstyle=1, /nodata, ytitle='Scale Factor', title='Chip Scaling relative to Chihiro', $
		xticks=6, yrange=yrange, ystyle=1, xtickname=[' ','chi','cla','fio','kik','nau',' ']
	oploterror, indgen(5), finaddit[0:4,0], (stddev(addit - finaddit,dimension=2)/sqrt(n_elements(addit[0,*])))[0:4], $
		psym=10, thick=5, errthick=3
	oplot, indgen(10)-1, intarr(10)+1, linestyle=2, thick=4
	sharpcorners, thick=4
	
	ploterror, indgen(5), finaddit[5:*,0], (stddev(addit - finaddit,dimension=2)/sqrt(n_elements(addit[0,*])))[5:*], $
		psym=10, xrange=[-1,5], xstyle=1, ytitle='Scale Factor', thick=5, errthick=3, $
		xtickname=[' ','pon','san','sat','she','sop',' '], xticks=6, yrange=yrange, ystyle=1
	sharpcorners, thick=4
	oplot, indgen(10)-1, intarr(10)+1, linestyle=2, thick=4
	cleanplot, /silent
	psclose
	
	
	;mwrfits, arr_struct({name:chip, scale:finaddit[*,0], err:stddev(addit - finaddit,dimension=2)/sqrt(n_elements(addit[0,*]))}), 'diag_supcamchipscale.fits', /create
	
	post = 'c'
	FOR i=0,(numExps - 1) DO BEGIN
	
		FOR k=0,(numChips - 1) DO BEGIN
			
			img = float(mrdfits(file[k,i], 0, hdr, /silent))				
			outFile = post + file[k,i]
			;SPAWN, 'rm ' + outFile
			file_delete, outFile, /allow
			badPix = where(~finite(img), countBPix)
			expTime = float(sxpar(hdr, 'EXPTIME'))
			if (expTime eq 0) then stop
			img = TEMPORARY(img) * finaddit[k,i] / expTime											;multiply by the scaling factor and divide by the exposure time
			if (countBPix gt 0) then img[badPix] = !values.f_nan
			sxaddpar, hdr, 'EXPTIME', 1., /saveComment, before='HISTORY'												;reset the exposure time to 1 second
			sxaddpar, hdr, 'OEXPTIME', exptime, 'Original exposure time, before rescaling. [s]', before='HISTORY'		;rename the original exposure time		
			sxaddpar, hdr, 'CHPSCALE', finaddit[k,i], 'Chip scaling factor [normalized to chihiro]', before='HISTORY'	;add in chip scaling keyword	
			saturate = sxpar(hdr, 'SATURATE')
			if (saturate eq 0) then stop
			gain = sxpar(hdr, 'GAIN')
			if (gain eq 0) then stop
			sxaddpar, hdr, 'SATURATE', saturate*finaddit[k,i]/exptime, 'Saturation/non-linearity level', before='HISTORY'	;add in saturation level keyword (~where non-linearity sets in)			
			sxaddhist, 'chipToChipQE.pro: Chip to chip QE differences corrected for.', hdr
			sxaddpar, hdr, 'GAIN', gain*exptime, /saveComment, before='HISTORY'
			sxaddhist, 'chipToChipQE.pro: Image divided by exposure time of ' + strtrim(exptime,2) + ' s', hdr
			mwrfits, FLOAT(img), outFile, hdr
		ENDFOR
	
	ENDFOR
	
	pre = post

END





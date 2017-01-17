FUNCTION scoreit, plot=plot
	;compile_opt idl2
	common share_c2cqe, imgs, xsize, ysize, wid, modx, i, scores, xxx, clafio, kikfio, sansat, shesat, chicla, ponsan, naukik, sopshe, satfio, sancla, shekik, ponchi, sopnau
	;common share_c2cqe
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


if 0 then begin	
	relscores = DBLARR(10)
	
	relScores[0] = abs(scores[0]) + abs(scores[8])
	relScores[1] = abs(scores[0]) + abs(scores[1]) + abs(scores[9])
	relScores[2] = abs(scores[1]) + abs(scores[2]) + abs(scores[10])
	relScores[3] = abs(scores[2]) + abs(scores[3]) + abs(scores[11])
	relScores[4] = abs(scores[3]) + abs(scores[12])
	relScores[5] = abs(scores[4]) + abs(scores[8])
	relScores[6] = abs(scores[4]) + abs(scores[5]) + abs(scores[9])
	relScores[7] = abs(scores[5]) + abs(scores[6]) + abs(scores[10])
	relScores[8] = abs(scores[6]) + abs(scores[7]) + abs(scores[11])
	relScores[9] = abs(scores[7]) + abs(scores[12])
endif	
	
	

;	med = median(scores)
;	sig = robust_sigma(scores)
;	tb = [9,10,6,11,12]	
;	out = where( (abs(scores[tb] - med) gt 3.*sig) , nout)		;find an up down outlier
;	if (nout gt 0) then begin
;	endif

	dev = abs(scores - median(scores))
	out = where( dev gt 3.*robust_sigma(scores) , nout)			;find outliers, for instance, interfaces where a large galaxy sits
	if (nout gt 0) then scores[out[where(dev[out] eq max(dev[out]))]] = 0.		;only get rid of the biggest outlier
		
	;plot, scores
	;print, scores
	;pause
	
	scores[6] = 0.
	return, total(abs(scores))
	

;RETURN, relScores	
END





FUNCTION newtIFy, X, plot=plot
	plot = 1
	common share_c2cqe
	modx = X	
	ret = scoreit(plot=plot)	
	;print, ret
	;print, xxx
	boundDist = abs(median(modx) - 1.)
	IF (boundDist gt 0.2) THEN ret = ret + boundDist 
	++xxx
RETURN, ret
END




FUNCTION fluxScale, img0, img1, scale
	gPix0 = WHERE(img0 ne -32768, countGPix0)
	gPix1 = WHERE(img1 ne -32768, countGPix1)
	;zero = biweight_mean(img0[gPix0])
	;one = biweight_mean(img1[gPix1])
	img0 = TEMPORARY(img0) * scale[0]	
	img1 = TEMPORARY(img1) * scale[1]		
	zero = biweight_mean(img0[gPix0])
	one = biweight_mean(img1[gPix1])
	
RETURN, 1d - (one / zero)
END




FUNCTION subimglr, img0, img1, wid, backwards=backwards
	
	if (n_elements(backwards) ne 0) then begin
		dum = img0
		img0 = img1
		img1 = dum
	endif

	
	for z=0,(n_elements(img1[0,*])-1) do begin
		ind = where(img1[*,z] ne -32768, count)
		if (count gt (wid + 1)) then d_dum = img1[(reverse(ind))[0:wid],z] else d_dum = dblarr(wid+1)-32768
		if (z eq 0) then dum1 = d_dum else dum1 = [[dum1],[d_dum]]
	endfor
	
	for z=0,(n_elements(img0[0,*])-1) do begin
		ind = where(img0[*,z] ne -32768, count)
		if (count gt (wid + 1)) then d_dum = img0[ind[0:wid],z] else d_dum = dblarr(wid+1)-32768
		if (z eq 0) then dum0 = d_dum else dum0 = [[dum0],[d_dum]]
	endfor

	len = min([n_elements(dum0[0,*]),n_elements(dum1[0,*])])-1

	if (n_elements(backwards) eq 0) then return, {reg0:dum0[*,0:len], reg1:dum1[*,0:len]} $
									else return, {reg0:dum1[*,0:len], reg1:dum0[*,0:len]}
END







FUNCTION fluxScalelr, data0, scale, plot=plot

	dum0 = data0.(0)
	gpix = where(dum0 ne -32768)
	dum0[gpix] = dum0[gpix]*scale[0]

	dum1 = data0.(1)
	gpix = where(dum1 ne -32768)
	dum1[gpix] = dum1[gpix]*scale[1]

	dum = median(dum1,dimension=1)/median(dum0,dimension=1)
	;dum = dum[0:min([n_elements(dum0[0,*]),n_elements(dum1[0,*])])-1]
	;dum = dum[where( (dum gt 0) and (dum lt 2.) )]
	dum = dum[where( (dum gt -32768)  )]


	;plot, dum, yrange=[-0.2,0.2]+median(dum)
	for xyz=0,10 do dum = dum[where( abs(dum - median(dum)) le 2.5*stddev(dum - median(dum)) )]
	;oplot, dum, color=fsc_color('red')
	
	mmm, dum, res

	if (res gt 0) then begin
		dum = median(dum0,dimension=1)/median(dum1,dimension=1)
		;dum = dum[0:min([n_elements(dum0[0,*]),n_elements(dum1[0,*])])-1]
		;dum = dum[where( (dum gt 0) and (dum lt 2.) )]
			dum = dum[where( (dum gt -32768)  )]

		;plot, dum, yrange=[-0.2,0.2]+median(dum)
		for xyz=0,10 do dum = dum[where( abs(dum - median(dum)) le 2.5*stddev(dum - median(dum)) )]
		;oplot, dum, color=fsc_color('red')	
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








FUNCTION subimgtb, img0, img1, wid
		
	for z=0,(n_elements(img1[*,0])-1) do begin
		ind = where(img1[z,*] ne -32768, count)
		if (count gt (wid + 1)) then d_dum = transpose(img1[z,(reverse(ind))[0:wid]]) else d_dum = dblarr(wid+1)-32768
		if (z eq 0) then dum1 = d_dum else dum1 = [[dum1],[d_dum]]
	endfor
	
	for z=0,(n_elements(img0[*,0])-1) do begin
		ind = where(img0[z,*] ne -32768, count)
		if (count gt (wid + 1)) then d_dum = transpose(img0[z,ind[0:wid]]) else d_dum = dblarr(wid+1)-32768
		if (z eq 0) then dum0 = d_dum else dum0 = [[dum0],[d_dum]]
	endfor

	len = min([n_elements(dum0[0,*]),n_elements(dum1[0,*])])-1

	return, {reg0:dum0[*,0:len], reg1:dum1[*,0:len]}
END







FUNCTION fluxScaletb, data0, scale, plot=plot

	dum0 = data0.(0)
	gpix = where(dum0 ne -32768)
	dum0[gpix] = dum0[gpix]*scale[0]

	dum1 = data0.(1)
	gpix = where(dum1 ne -32768)
	dum1[gpix] = dum1[gpix]*scale[1]

	dum = median(dum1,dimension=1)/median(dum0,dimension=1)
	;dum = dum[0:min([n_elements(dum0[0,*]),n_elements(dum1[0,*])])-1]
	;dum = dum[where( (dum gt 0) and (dum lt 2.) )]
		dum = dum[where( (dum gt -32768)  )]

	for xyz=0,10 do dum = dum[where( abs(dum - median(dum)) le 2.5*stddev(dum - median(dum)) )]

	mmm, dum, res

	if (res gt 0) then begin
		dum = median(dum0,dimension=1)/median(dum1,dimension=1)
		;dum = dum[0:min([n_elements(dum0[0,*]),n_elements(dum1[0,*])])-1]
		;dum = dum[where( (dum gt 0) and (dum lt 2.) )]
			dum = dum[where( (dum gt -32768)  )]

		
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





PRO chiCla, plot=plot
	common share_c2cqe
	k = 1
	l = k - 1
	m = 0
	scale = [modx[k], modx[l]]
	
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		chicla = subimglr(img0, img1, wid)
	endif 

	scores[m] = fluxscalelr(chicla, scale, plot=plot)

END


PRO claFio, plot=plot
	common share_c2cqe
	k = 2
	l = k - 1
	m = 1
	scale = [modx[k], modx[l]]
	
	if (xxx eq 0) then begin
		img0 = imgs.(k)
		img1 = imgs.(l)
		claFio = subimglr(img0, img1, wid)
	endif 

	scores[m] = fluxscalelr(claFio, scale, plot=plot)

END


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



























PRO chip2chipqe, files, pre



;uncomment these to go back to the .mos file method
;PRO chip2chipqe, inFile, pre
;readcol, inFile, files, xshift, yshift, a, fluxin, format='a,d,d,d,d', /silent					;read in the mos file
;files = pre + files



common share_c2cqe																					;common block


starttime0 = SYSTIME(1)

numFiles 	= N_ELEMENTS(files)																	;number of input image
numChips 	= 10																				;number of chips in the mosaic
numExps 	= numFiles/numChips																	;number of separate exposures


file = strarr(numChips, numExps)
chip = ['chihiro', 'clarisse', 'fio', 'kiki', 'nausicaa', 'ponyo', 'san', 'satsuki', 'sheeta', 'sophie']	;names of each detector

xsize = INTARR(numChips, numExps)
ysize = INTARR(numChips, numExps)


FOR k=0,(numChips - 1) DO BEGIN			;cycle through each chip and ...
	j = 0
	FOR i=0,(numFiles - 1) DO BEGIN
		IF (N_ELEMENTS(STRSPLIT(files[i], chip[k], /regex)) gt 1) THEN begin
			file[k,j] = files[i]
			hdr = headfits(file[k,j])
			xsize[k,j] = sxpar(hdr, 'NAXIS1')
			ysize[k,j] = sxpar(hdr, 'NAXIS2')		
			j++
		ENDIF
	ENDFOR
ENDFOR


;goto, jump1

scores = DBLARR(13)
modx = DBLARR(10) + 1

IF (N_ELEMENTS(wid) EQ 0) THEN wid = 50


addit = FLTARR(numChips, numExps) + 1.
date = fix(strtrim(strsplit(sxpar(headfits(file[0]), 'DATE-OBS'),'-',/extract)))
date = date[0] + date[1]/12. + date[2]/365.

	
FOR i=0,(numExps - 1) DO BEGIN
	
	FOR k=0,(numChips - 1) DO BEGIN
		IF (k EQ 0) THEN imgs = CREATE_STRUCT('num' + STRTRIM(k,1), float(mrdfits(file[k,i], 0, hdr, /silent)) ) $
					ELSE imgs = struct_addtags(imgs, CREATE_STRUCT('num' + STRTRIM(k,1), float(mrdfits(file[k,i], 0, hdr, /silent))))

		IF (k EQ 0) THEN hdrs = CREATE_STRUCT('num' + STRTRIM(k,1), hdr) $
					ELSE hdrs = struct_addtags(hdrs, CREATE_STRUCT('num' + STRTRIM(k,1), hdr))

		xsize[k,i] = N_ELEMENTS(imgs.(k)[*,0])
		ysize[k,i] = N_ELEMENTS(imgs.(k)[0,*])
		
	ENDFOR
	
	IF (i gt 0) THEN BEGIN
					x = addit[*,i-1]
				ENDIF ELSE BEGIN
					if (date le 2009.) then begin									;for older data with lingering non-linearity issues, allow a different scaling for different exposure times
						x = float([1.00000, 0.975279, 0.940946, 0.975610, 0.954359, 1.06361, 0.936070, 0.940543, 0.963282, 1.00205])
					endif else x = float([1.00000, 0.962169, 0.933242, 0.970565, 1.01779, 1.09830, 0.960591, 0.939223, 0.979983, 1.03439])	;initial guesses at the scale, use last images values IF possible
				ENDELSE

	;1.00000     0.975279     0.940946     0.975610     0.954359      1.06361     0.936070     0.940543     0.963282      1.00205		;short exposure 2008
	; 1.00000     0.988884     0.945970     0.989706     0.963803      1.08250     0.940372     0.951366     0.975463     0.985745		;long exposure 2008

;x = DBLARR(numChips) + 1
	xxx = 0
	;xx = tnmin('newtify', x, bestmin=bestmin, /autoderivative)
	;plot, scores
	;stop
	
	;addit[*,i] = NEWTON(x, 'newtify', /double)		;NEWTON-raphson solver
	time, itime
	addit[*,i] = tnmin('newtify', x, bestmin=bestmin, /autoderivative)	
	time, itime
	stop
	addit[*,i] = addit[*,i]/addit[0,i]
	
	
	
	;addit[*,i] = addit[0,i] / addit[*,i]
	;PRINT, addit[*,i]
	if 0 then begin
	FOR k=0,(numChips - 1) DO BEGIN
		;++k
		outFile = 'c' + file[k,i]
		SPAWN, 'rm ' + outFile
		badPix = WHERE(imgs.(k) EQ -32768, countBPix)
		imgs.(k) = TEMPORARY(imgs.(k)) * addit[k,i] / ((exptime = float(sxpar(hdrs.(k), 'EXPTIME'))))			;multiply by the scaling factor and divide by the exposure time
		IF (countBPix GT 0) THEN imgs.(k)[badPix] = -32768
		sxaddpar, hdrs.(k), 'EXPTIME', 1., /saveComment															;reset the exposure time to 1 second
		sxaddhist, 'chipToChipQE.pro: Chip to chip QE differences corrected for.', hdrs.(k)
		sxaddhist, 'chipToChipQE.pro: Image divided by exposure time of ' + strtrim(exptime,2) + ' s', hdrs.(k)
		mwrfits, FLOAT(imgs.(k)), outFile, hdrs.(k)
		;print, k
	ENDFOR
	endif

	print, addit

	print, 'i = ', i
ENDFOR

;result = TOTAL(addit, 2)/float(numExps)

;if (numexps gt 1) then finaddit = median(addit, dimension=2, /even) else finaddit = addit														;final multiplicative scaling between chips




if (date le 2009.) then begin									;for older data with lingering non-linearity issues, allow a different scaling for different exposure times
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
		endif									;otherwise just leave finaddit column as is
	endfor	
endif else finaddit = rebin(median(addit, dimension=2, /even), [numchips, numexps])							;for newer data just take the median


;erase
;cleanplot
;psopen, 'diag_chipscaling', /encapsulated, /color, /heiles, xs=8, ys=6, /inches
;!p.multi=[0,1,1]
;!p.charsize = 1.5
;!p.charthick = 3
;ploterror, indgen(10), finaddit[*,0], stddev(addit - finaddit,dimension=2)/sqrt(n_elements(addit[0,*])), $
;	psym=10, xrange=[-1,10], xstyle=1, xtitle='Chip', ytitle='Scale Factor', title='Chip Scaling relative to Chihiro', thick=5, errthick=3, $
;	xtickname=[' ','chi','cla','fio','kik','nau','pon','san','sat','she','sop',' '], xticks=11
;sharpcorners, thick=4
;psclose


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


mwrfits, arr_struct({name:chip, scale:finaddit[*,0], err:stddev(addit - finaddit,dimension=2)/sqrt(n_elements(addit[0,*]))}), 'diag_supcamchipscale.fits', /create

;jump1:
;finaddit = rebin((mrdfits('supcamchipscale.fits',1,/silent)).scale, [numchips,numexps])

;if newer data then use finaddit, but if older data, then use finaddit for same exposure times and addit otherwise?

FOR i=0,(numExps - 1) DO BEGIN

	FOR k=0,(numChips - 1) DO BEGIN
		
		img = float(mrdfits(file[k,i], 0, hdr, /silent))				
		outFile = 'c' + file[k,i]
		SPAWN, 'rm ' + outFile
		badPix = WHERE(img EQ -32768, countBPix)
		img = TEMPORARY(img) * finaddit[k,i] / ((exptime = float(sxpar(hdr, 'EXPTIME'))))							;multiply by the scaling factor and divide by the exposure time
		IF (countBPix GT 0) THEN img[badPix] = -32768
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

pre = 'c'

END





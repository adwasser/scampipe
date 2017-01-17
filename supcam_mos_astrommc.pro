

FUNCTION propDens, x, seed
	
	nChips = 10
	mult = [ 	1d-6, $						;CDELT1
				1d-6, $						;CDELT2
				0d, $						;CRVAL1
				0d, $						;CRVAL2
				replicate(0.1d,nChips), $	;CRPIX1
				replicate(0.1d,nChips), $	;CRPIX2
				replicate(1d-4,nChips) ]	;CROTA2

	nChips = 10
	mult = [ 	1d-8, $						;CDELT1
				1d-8, $						;CDELT2
				1d-15, $						;CRVAL1
				1d-15, $						;CRVAL2
				replicate(0.5d,nChips), $	;CRPIX1
				replicate(0.5d,nChips), $	;CRPIX2
				replicate(1d-4,nChips) ]	;CROTA2

	
	nPar = n_elements(x)
	par = dblarr(nPar)
	for i=0,(nPar - 1) do par[i] = randomn(seed)*mult[i] + x[i]
	return, par

END

;-----------------------------------------------------------------------------------------------

FUNCTION jitter, x, mult=mult
	
	nChips = 10
	if (n_elements(mult) eq 0) then begin
		mult = [ 	1d-7, $						;CDELT1
					1d-7, $						;CDELT2
					0d-4, $						;CRVAL1
					0d-4, $						;CRVAL2
					replicate(100d,nChips), $	;CRPIX1
					replicate(100d,nChips), $	;CRPIX2
					replicate(0d-4,nChips), $	;CROTA2
					1d-1]
	endif
	
	nPar = n_elements(x)
	par = dblarr(nPar)
	for i=0,(nPar - 1) do par[i] = (randomu(seed) - 0.5)*mult[i]
	return, par

END

;-----------------------------------------------------------------------------------------------

FUNCTION getParStruc, param0
	param = param0
	cdelt1 = param[0]
	cdelt2 = param[1]
	crval1 = param[2]
	crval2 = param[3]
	crpix1 = param[4:13]
	crpix2 = param[14:23]
	crota2 = param[24:33]
	angle = param[34]
	param = {CDELT1:CDELT1, CDELT2:CDELT2, CRVAL1:CRVAL1, CRVAL2:CRVAL2, CRPIX1:CRPIX1, CRPIX2:CRPIX2, CROTA2:CROTA2, ANGLE:ANGLE}
	return, param	
END

;-------------------------------------------------------------------------------------------

PRO supcam_mos_astrom_writeastrom, data, filenames, nfiles, astr, post=post
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
		sxaddpar, hdr, 'CROTA2', (astr.angle + astr.crota2[i])*!radeg, /saveComment, BEFORE='HISTORY'
	
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

PRO supcam_mos_astrom_plot_resid, res, ps=ps, outfile=outfile

	if (n_elements(ps) gt 0) then begin
		psopen, outfile, /encapsulated, /color, xs=8, ys=6, /inches
		!p.charsize = 1.2
		!p.charthick = 3
	endif

	xi = res[*,0]
	eta = res[*,1]
	xiDif = res[*,2]
	etaDif = res[*,3]


	plotsym, 0, 0.4, /fill
	cleanplot, /silent
	erase
	!p.multi=[0,2,2]
	yrange = [-1,1]*0.59
	multiplot
	plot, xi, xiDif, xrange=reverse(minmax(xi)), yrange=yrange, xstyle=1, ystyle=1, psym=8, ytitle=textoidl('\Delta \chi'), title=textoidl('\chi_{RMS} = ')+roundx(robust_sigma(xiDif),3)+'"'
	oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
	sharpcorners, thick=3
	multiplot
	plot, eta, xiDif, xrange=minmax(eta), yrange=yrange, xstyle=1, ystyle=1, psym=8, title=textoidl('\eta_{RMS} = ')+roundx(robust_sigma(etaDif),3)+'"'
	oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
	sharpcorners, thick=3
	multiplot
	plot, xi, etaDif, xrange=reverse(minmax(xi)), yrange=yrange, xstyle=1, ystyle=1, psym=8, ytitle=textoidl('\Delta \eta'), xtitle=textoidl('\chi [degrees]')
	oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
	sharpcorners, thick=3
	multiplot
	plot, eta, etaDif, xrange=minmax(eta), yrange=yrange, xstyle=1, ystyle=1, psym=8, xtitle=textoidl('\eta [degrees]')
	oplot, (findgen(10)-5)*1000, intarr(10), linestyle=2, color=fsc_Color('red')
	sharpcorners, thick=3

	if (n_elements(ps) gt 0) then psclose 
	cleanplot, /silent
END

;-------------------------------------------------------------------------------------------

FUNCTION getGZ, seed, a, nSamp
	if (n_elements(nSamp) eq 0) then nSamp = 1l	
	;a = 2d
	x0 = 1d/a
	x1 = a
	alpha = -0.5d
	beta = 1d + alpha
	gz = (randomu(seed,nSamp) * (x1^beta - x0^beta) + x0^beta )^(1d/beta)
	if (nSamp eq 1) then gz = gz[0]
	return, gz
END		 

;-------------------------------------------------------------------------------------------

FUNCTION model, data, param, res=res
	p = getParStruc(param)
	nChips = 10
	res = dblarr(1d4,4)
	k = 0l
	for i=0,(nChips - 1) do begin
		data_chip = data.(i)
		wcssph2xy, data_chip.ra, data_chip.dec, xiw, etaw, 2, latpole = 90.0, crval=[p.crval1,p.crval2]		;convert RA/DEC to xiw and etaw
		xyrot, data_chip.x - p.crpix1[i], data_chip.y - p.crpix2[i], p.angle + p.crota2[i], xi, eta, mult=[p.cdelt1,p.cdelt2]	;convert X/Y to xi and eta

		;WCSXY2SPH, xi, eta, ra, dec, CRVAL = [p.crval1,p.crval2], latpole=90.								;
		;dif = [xi - xiw, eta - etaw]*3600d											;calculate the deviation between each set of xi and eta
		res[k,0] = [[xiw], [etaw], [(xi - xiw)*3600d], [(eta - etaw)*3600d]]
		k += n_elements(data_chip)
	endfor
	res = res[0:k,*]
	difs = [res[*,2],res[*,3]]
	;plot, xiDifs, etaDifs, psym=1, /iso, xrange=[-1,1]*50, yrange=[-1,1]*50
	;pause
	;difs = [xiDifs, etaDifs]
	return, difs
END

;-------------------------------------------------------------------------------------------

FUNCTION likelihood, data, param, noFlag
	sigmaSquared = 0.5^2d
	noFlag = 1
	nWalkers = n_elements(param[0,*])
	prob = dblarr(nWalkers)
	for i=0l,(nWalkers - 1l) do begin
		model = model(data, param[*,i])
	
		;pr = losvd_mod(model, [0d,sigma,0d,0d])
		;stop	
		pr = model^2d
		;pr = (pr > 1d-323)			;gaussian likelihood		
		
		;prob[i] = total(alog(pr), /nan)
		prob[i] = -0.5d*total(pr)/sigmaSquared
		if ~finite(prob[i]) then stop
	endfor
	;plot, prob
	;stop
	if (nWalkers eq 1) then prob = prob[0]
	;if (prob lt -4.075*10^5d) then stop
	return, prob
END

;-----------------------------------------------------------------------------------------------

FUNCTION getParams, data, date, filenames=filenames, start=start
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
			getrot, headfits(filenames[0]), imrot, cdelt
			ANGLE = imRot*!dtor
			CROTA2 = CROTA2 * 0d
		endif else begin
			if (strmid(filenames[0],0,1) eq 'n') then begin
			endif else begin
			endelse
			CRPIX1 = [5425.8691,3351.0452,1241.0316,-871.95972,-2964.5720,5423.8774,3351.3420,1242.9849,-877.50110,-2966.1489]
			CRPIX2 = [111.32944,111.39816,111.15529,113.26266,113.70063,4315.5620,4330.3096,4338.3032,4330.0840,4318.5273]
			CROTA2 = [-0.0055500342,-0.0053783910,-0.0054863331,-0.0053572311,-0.0051588219,-0.0050752364,-0.0053784340,-0.0049788447,-0.0052406706,-0.0055207304]
			getrot, headfits(filenames[0]), imrot, cdelt
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

	param = {CDELT1:CDELT1, CDELT2:CDELT2, CRVAL1:CRVAL1, CRVAL2:CRVAL2, CRPIX1:CRPIX1, CRPIX2:CRPIX2, CROTA2:CROTA2, ANGLE:ANGLE}
	
	
	angleRange = !pi/(sqrt(findgen(10)+1))^3d
	angleRange = dblarr(3) + 1d-1	;angleRange[0]
	; testing this out
	;angleRange = dblarr(3) + 6d   ;angleRange[0]
	;angleRange = dblarr(3) + 6d

a = 0	
while (a lt n_elements(angleRange)) do begin
	; testing with this off
	for i=0,10 do begin

		param = [CDELT1, CDELT2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CROTA2, ANGLE]
	
		z = model(data, param, res=res)	
		if (i eq 0) then res0 = res
		;window, 2
		;plot, res[*,0], res[*,2], psym=1, title=roundx(likelihood(data, param),2)
		;window, 0
		;plot, res[*,0], res[*,3], psym=1, title=roundx(likelihood(data, param),2)
		;oplot, res0[*,0], res0[*,3], psym=1, color=fsc_color('red')
		
		;stop
		;pause
		
		;print, likelihood(data, param)
		
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
	;plot, angles*!radeg, like
	;wset, 0
	



	dummy = max(like,imax)
	if ((imax lt 200) or (imax gt 300)) then stop
	;print, dummy
	ANGLE = (angles[imax])[0]
	param = [CDELT1, CDELT2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CROTA2, ANGLE]

	;pause

++a
endwhile
	
	;stop	

	return, param
END

;-----------------------------------------------------------------------------------------------

FUNCTION prior, y, noFlag
	badVal = -1d99
	nWalkers = n_elements(y[0,*])
	
	noFlag = finite(nWalkers)*0 + 1
	prob = dblarr(nWalkers)

	val = y[0,*]	;CDELT1
	iBad = where((abs(val) lt 5.5d-5) or (abs(val) gt 5.7d-5),nBad)
	if (nBad gt 0) then noFlag[iBad] *= 0 
	
	val = y[1,*]	;CDELT2
	iBad = where((abs(val) lt 5.5d-5) or (abs(val) gt 5.7d-5),nBad)
	if (nBad gt 0) then noFlag[iBad] *= 0
		
	val = y[2,*]	;CRVAL1
	iBad = where((val lt 0d) or (val gt 360d),nBad)	
	if (nBad gt 0) then noFlag[iBad] *= 0 

	val = y[3,*]	;CRVAL2
	iBad = where((val lt -90d) or (val gt 90d),nBad)	
	if (nBad gt 0) then noFlag[iBad] *= 0 

	if (nWalkers eq 1) then prob = prob[0]
		
	return, prob
END

;-----------------------------------------------------------------------------------------------

FUNCTION wander, nIters, data, start, nParams, seed, nWalkers=nWalkers, stretch=stretch, accep=accep, $
		prob=prob_data, burn=burn, mult=mult
	
	;nParams = n_elements(start)
	;nWalkers = 10

	if (nWalkers gt 1) then begin
		device, decomposed=0
		loadct, 39, /silent
		maxcolor = 200
		colors = round(findgen(nWalkers)*maxColor/(float(nWalkers) - 1) + (255 - maxColor)/2d)
	endif else colors = 255



	;a = 1.1d																	;parameter controlling the step
	x = dblarr(nParams,nWalkers,nIters)

	;x[0,0,0] = cmreplicate(start,nWalkers)										;use the initial guesses as the starting point
	if (n_elements(burn) gt 0) then begin
		x[*,*,0] = cmreplicate(start,nWalkers)									;use the initial guesses as the starting point
		for i=1,(nWalkers - 1) do $
			x[*,i,0] += jitter(start, mult=mult)
		;for i=1,(nWalkers - 1) do $
		;	x[*,i,0] += randomn(seed,4)*1d-5
	endif else begin
		x[*,*,0] = cmreplicate(start,nWalkers)									;use the initial guesses as the starting point
		for i=1,(nWalkers - 1) do $
			x[*,i,0] += jitter(start)*1d-2
		
		
		;x[*,*,0] = start
	endelse






	prior = dblarr(nWalkers,nIters)
	prior[0,0] = prior(x[*,*,0])												;set the initial prior values

	prob_data = dblarr(nWalkers,nIters)
	prob_data[0,0] = likelihood(data, x[*,*,0])									;probability of data given parameters x[*,t]



	accep = intarr(nWalkers,nIters)
	rArr = fltarr(nWalkers,nIters)
	indWalk = lindgen(nWalkers)
	k = 0l
	for t=0l,(nIters - 2l) do begin
	
		if (nWalkers eq 1) then begin

			for k=0,(nWalkers - 1) do begin
	
	
				y = propDens(x[*,k,t], seed)											;get the new proposed mu
				prior_y = prior(y, noFlag)												;prior on y
				
				
			
				if noFlag then begin
					if (y[1] lt 0) then stop

					
					prob_data_y = likelihood(data, y, noFlag)							;probability of data given proposed parameters y
					r = min([1d,exp(prob_data_y + prior_y - prob_data[k,t] - prior[k,t])])	;acceptance ratio
					;print, r
					
					rArr[k,t] = r
					u = randomu(seed)													;random number between 0 and 1
					accep[k,t] = (u lt r) and noFlag
				endif else accep[k,t] = 0
				if accep[k,t] then begin
					if (y[1] gt 500d) then stop
					x[*,k,t+1] = y 														;if u lt r, then accept the proposed parameters y
					prior[k,t+1] = prior_y												;and the prior on y
					prob_data[k,t+1] = prob_data_y										;and the probability of data given proposed parameter y
				endif else begin
					x[*,k,t+1] = x[*,k,t]												;otherwise, keep the same parameters
					prior[k,t+1] = prior[k,t]											;and the same prior still applies
					prob_data[k,t+1] = prob_data[k,t]		
				endelse
			endfor
		
			;if ((t mod 10) eq 0) then plot, x[0,*,t], x[1,*,t], psym=1, xrange=[-1,1]*200, yrange=[0,1]*200

			if 0*((t mod 500) eq 0) then begin
				plot, x[0,*,0:t], x[1,*,0:t], /nodata
				oplot, x[0,0,0:t], x[1,0,0:t]
				if (nWalkers gt 1) then oplot, x[0,1,0:t], x[1,1,0:t], color=fsc_color('red')			
			endif


		endif else begin
			
			for k=0,(nWalkers - 1) do begin
				ind = indWalk[where(indWalk ne k)]										;get the indices of the other walkers
				ind = ind[sort(randomu(seed,nWalkers-1l))]								;randomly sort them
				ind = ind[0]															;select one (effectively at random)

				xj = x[*,ind,t]															;select the position of a random walker
				
				z = getGZ(seed, stretch, 100)														;draw a random number from a power law
				
				y = xj + z*(x[*,k,t] - xj)												;take a step towards the other walker
				
				prior_y = prior(y, noFlag)												;prior on y
				

				
				if noFlag then begin
					
					prob_data_y = likelihood(data, y, noFlag)							;probability of data given proposed parameters y
		
					;r = min([1d,exp(prob_data_y + prior_y - prob_data[t] - prior[t])])	;acceptance ratio
					r = min([1d,(z^(nParams - 1d))*exp(prob_data_y + prior_y - prob_data[k,t] - prior[k,t])])	;acceptance ratio
					rArr[k,t] = r
					;r = min([1, (prob_data_y * prior_y) / (prob_data_xt * prior_xt)])	;acceptance ratio
					u = randomu(seed)													;random number between 0 and 1
					accep[k,t] = (u le r) and noFlag
				endif else accep[k,t] = 0
		
				if accep[k,t] then begin
					if (y[1] gt 500d) then stop
					x[*,k,t+1] = y 														;if u lt r, then accept the proposed parameters y
					prior[k,t+1] = prior_y												;and the prior on y
					prob_data[k,t+1] = prob_data_y										;and the probability of data given proposed parameter y
				endif else begin
					x[*,k,t+1] = x[*,k,t]												;otherwise, keep the same parameters
					prior[k,t+1] = prior[k,t]											;and the same prior still applies
					prob_data[k,t+1] = prob_data[k,t]		
				endelse
			endfor
			
			;pause
		endelse
		statusbar, t, nIters
		if 0*(((t mod 50) eq 0) and (t gt 0)) then begin
			if (t eq 50) then window, 2 else wset, 2
			plot, x[4,0,0:t-1], yrange=minmax(x[4,*,0:t-1])
			for w=1,(nWalkers - 1) do oplot, x[4,w,0:t-1], color=colors[w]

			if (t eq 50) then window, 3 else wset, 3
			plot, x[-1,0,0:t-1], yrange=minmax(x[-1,*,0:t-1])
			for w=1,(nWalkers - 1) do oplot, x[-1,w,0:t-1], color=colors[w]

			if (t eq 50) then window, 1 else wset, 1
			plot, x[-2,0,0:t-1], yrange=minmax(x[-2,*,0:t-1])
			for w=1,(nWalkers - 1) do oplot, x[-2,w,0:t-1], color=colors[w]


			
			if (t eq 50) then window, 0 else wset, 0
			
			bestProb = max(prob_data[*,0:t],imax)
			bestPar = x[*, imax mod n_elements(prob_data[*,0]), imax / n_elements(prob_data[*,0])]
			
			if (likelihood(data, bestPar) ne bestProb) then stop
			
			z = model(data, bestPar, res=res)
			
	 		supcam_mos_astrom_plot_resid, res, ps=ps, outfile=outfile
					
						
			;plot, prob_data[0,0:t-1], yrange=minmax(prob_data[*,0:t-1])
			;for w=1,(nWalkers - 1) do oplot, prob_data[w,0:t-1], color=colors[w]			
			
		endif
	
	endfor
	statusbar, t, nIters
	cleanplot, /silent
	;erase
	return, x
END

;-----------------------------------------------------------------------------------------------

PRO makePlot, x, accep, data, z=z, prob=prob
	nWalkers = n_elements(x[0,*,0])
	if (nWalkers gt 1) then begin
		device, decomposed=0
		loadct, 39, /silent
		maxcolor = 200
		colors = round(findgen(nWalkers)*maxColor/(float(nWalkers) - 1) + (255 - maxColor)/2d)
	endif else colors = 255

	window, 1
	plot, x[4,*,*], x[14,*,*], /iso, psym=-3, /nodata, xrange=[-1,1]*6d3, yrange=[-1d3,5d3]
	for i=0,(nWalkers - 1) do for j=4,13 do oplot, x[j,i,*], x[j+10,i,*], color=colors[i]

	
	window, 2
	bin = 0.6
	hist = histogram(float(accep[0,*]),bin=bin)
	plothist, accep[0,*], bin=bin, title=roundx(hist[1]/double(max(hist)),2)
	for i=0,(nWalkers - 1) do plothist, accep[i,*], bin=bin, color=colors[i], /overplot


	window, 3
	plot, x[14,*,*], x[24,*,*], psym=3                                   
	oplot, x[14,*,500:*], x[24,*,500:*], psym=3, color=fsc_color('red')	
	
	stop
	
	
	param = x[*,0,-1]
	p = getParStruc(param)
	nChips = 10
	for i=0,(nChips - 1) do begin
		wcssph2xy, data.(i).ra, data.(i).dec, xiw, etaw, 2, latpole = 90.0, crval=[p.crval1,p.crval2]		;convert RA/DEC to xiw and etaw
		xyrot, data.(i).x - p.crpix1[i], data.(i).y - p.crpix2[i], p.crota2[i], xi, eta, mult=[p.cdelt1,p.cdelt2]	;convert X/Y to xi and eta
		;WCSXY2SPH, xi, eta, ra, dec, CRVAL = [p.crval1,p.crval2], latpole=90.								;
		dif = [xi - xiw, eta - etaw]*3600d											;calculate the deviation between each set of xi and eta
		if (i eq 0) then difs = dif else difs = [difs, dif]
	endfor

	
	return

END

;-----------------------------------------------------------------------------------------------

FUNCTION astrom_sigClip, data, res, nSig
	sigChi = robust_sigma(res[*,2])
	sigEta = robust_sigma(res[*,3])

	k = 0l
	for t=0,9 do begin
		nData = n_elements(data.(t))
		keep = intarr(nData)
		for i=0,(nData - 1) do begin
			if ((abs(res[k,2]) lt nSig*sigChi) and (abs(res[k,3]) lt nSig*sigEta)) then keep[i] = 1
			++k
		endfor
		hold = data.(t)
		hold = hold[where(keep,nKeep)]
		if (nKeep lt 10) then stop
		dataNew = struct_addtags(dataNew, create_struct('NUM'+roundx(t), hold))
	endfor
	return, dataNew
END

;-----------------------------------------------------------------------------------------------

FUNCTION supcam_mos_astrommc, data, date, filenames=filenames, maxiter=maxiter, start=start, distcorr=distcorr, post=post
	COMPILE_OPT idl2, HIDDEN 
	if (n_elements(maxiter) eq 0) then maxiter = 2
	;common sharemos
		
	data0 = data
	if (n_elements(distcorr) gt 0) then begin
		for i=0,(n_tags(data)-1) do begin
			hdr = headfits(filenames[i])
			chip = strtrim(sxpar(hdr, 'DETECTOR'),2)
			supcam_distcorr, data.(i).x, data.(i).y, chip, date, x1, y1
			data0.(i).x = x1
			data0.(i).y = y1
			stop
		endfor
		;stop
	endif
	
	;initial guesses at the offsets and rotations
	start = getParams(data, date, filenames=filenames, start=start)
	z = model(data, start, res=res)	
	;supcam_mos_astrom_plot_resid, res
	;stop
	
	nParams = n_elements(start)

	nBurn = 3d3
	nWalkers = 10
	stretch = 1.05
	;run the affine invariant MCMC
	x = wander(nBurn, data, start, nParams, seed, nWalkers=nWalkers, stretch=stretch, accep=accep, prob=prob, mult=mult, /burn)


	;find the maximum likelihood
	maxLikelihood = max(prob,imax)
	;the best fitting parameters
	bestPar = x[*, imax mod n_elements(prob[*,0]), imax / n_elements(prob[*,0])]
	z = model(data, bestPar, res=res)
	
	;testing out turning off sigclip
	dataNew = astrom_sigClip(data, res, 5.0)
	
	start = bestPar
	nChips = 10
	mult = [ 	1d-7, $						;CDELT1
				1d-7, $						;CDELT2
				0d-4, $						;CRVAL1
				0d-4, $						;CRVAL2
				replicate(10d,nChips), $	;CRPIX1
				replicate(10d,nChips), $	;CRPIX2
				replicate(0d-4,nChips), $	;CROTA2
				1d-1]
	
	nBurn = 5d3
	x = wander(nBurn, dataNew, start, nParams, seed, nWalkers=nWalkers, stretch=stretch, accep=accep, prob=prob, mult=mult, /burn)

;	makePlot, x, accep, data, prob=prob

	;find the maximum likelihood
	maxLikelihood = max(prob,imax)

	;the best fitting parameters
	bestPar = x[*, imax mod n_elements(prob[*,0]), imax / n_elements(prob[*,0])]
	z = model(dataNew, bestPar, res=res)	
	bestPar = getParStruc(bestPar)										;store the best parameters to a structure
	
	if ((nfiles = n_elements(filenames)) gt 0) then supcam_mos_astrom_writeastrom, data, filenames, nfiles, bestPar, post=post

	supcam_mos_astrom_plot_resid, res, /ps, outfile='astrometry_'+strmid(filenames[0],0,(strsplit(filenames[0],'_'))[-1]-1)
	
	dummy = file_search('dir_astrometry', /test_directory, count=nMatch)
	if ~nMatch then spawn, 'mkdir dir_astrometry'
	spawn, 'mv astrometry*eps dir_astrometry/'
	
	start = bestPar
	return, bestPar

END



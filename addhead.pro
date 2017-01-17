;+
;PURPOSE
;	Add header parameters to suprime-cam images processed with supcampipe.pro
;SYNTAX
;	addhead, mosaic_image, component_images
;Written by Jacob A. Arnold, 12-21-09, UCSC
;
;-


pro addhead, mosaic, files, hdr=hdr, xoff=xoff, yoff=yoff


numFiles = n_elements(files)														;number of files used in making mosaic


ind = indgen(numFiles)																;generate an index for each file
ind = ind[where( (ind mod 10) eq 0, numExp)]										;isolate just the first file in each mosaic (chihiro)


if (n_elements(xoff) eq 0) then xoff = 0.
if (n_elements(yoff) eq 0) then yoff = 0.



exptimes  = strarr(numExp)															;exposure times
airmasses = strarr(numExp)															;airmasses
obsDates  = strarr(numExp)															;observation dates
UTs		  = strarr(numExp)															;universal times at middle of exposure
HSTs	  = strarr(numExp)															;hawaii standard times at middle of exposure
LSTs	  = strarr(numExp)															;local (sidereal time?) at middle of exposure
MJDs	  = strarr(numExp)															;modified julian day at middle of exposure

header = headfits(mosaic)															;header of the mosaic image

for i=0,(numExp - 1) do begin														;cycle through each mosaic, just looking at first ccd

	hdr = headfits(files[ind[i]])													;read in header

	if (i eq 0) then begin															;for just the very first ccd record the following
		sxaddhist, 'Observation params from ' + files[ind[i]], header				;add history to the header


		SXADDPAR, header, 'EXPTIME', strtrim(sxpar(hdr, 'EXPTIME', comment=cexptime)), cexptime, /SaveComment 				
		oexptime = sxpar(hdr, 'OEXPTIME', comment=coexptime)
		if (oExpTime gt 0) then sxaddpar, header, 'OEXPTIME', oexptime, coexptime, /saveComment		
		SXADDPAR, header, 'FILTER01', strtrim(sxpar(hdr, 'FILTER01', comment=cfilter)), cfilter, /SaveComment                            
		SXADDPAR, header, 'OBJECT', strtrim(sxpar(hdr, 'OBJECT', comment=cobject)), cobject, /SaveComment                            
		SXADDPAR, header, 'OBSERVER', strtrim(sxpar(hdr, 'OBSERVER', comment=cobserver)), cobserver, /SaveComment                            
		SXADDPAR, header, 'PROP-ID', strtrim(sxpar(hdr, 'PROP-ID', comment=cpropid)), cpropid, /SaveComment                            
		SXADDPAR, header, 'RA2000', strtrim(sxpar(hdr, 'RA2000', comment=cra2000)), cra2000, /SaveComment                            
		SXADDPAR, header, 'DEC2000', strtrim(sxpar(hdr, 'DEC2000', comment=cdec2000)), cdec2000, /SaveComment 
		SXADDPAR, header, 'CRVAL1', float(strtrim(sxpar(hdr, 'CRVAL1', comment=ccrval1))), ccrval1, /SaveComment 
		SXADDPAR, header, 'CRVAL2', float(strtrim(sxpar(hdr, 'CRVAL2', comment=ccrval2))), ccrval2, /SaveComment 		                           
		SXADDPAR, header, 'CRPIX1', float(strtrim(sxpar(hdr, 'CRPIX1', comment=ccrpix1))) + xoff, ccrpix1, /SaveComment 
		SXADDPAR, header, 'CRPIX2', float(strtrim(sxpar(hdr, 'CRPIX2', comment=ccrpix2))) + yoff, ccrpix2, /SaveComment 
		SXADDPAR, header, 'OBSERVAT', strtrim(sxpar(hdr, 'OBSERVAT', comment=cobservat)), cobservat, /SaveComment                            
		SXADDPAR, header, 'TELESCOP', strtrim(sxpar(hdr, 'TELESCOP', comment=ctelescop)), ctelescop, /SaveComment 
		SXADDPAR, header, 'NUMEXPS', numEXP, 'Number of exposures used to make final mosaic.', /SaveComment 
		saturate = sxpar(hdr, 'SATURATE', comment=csaturate)
		if (saturate gt 0) then sxaddpar, header, 'SATURATE', saturate, csaturate, /saveComment
		gain = sxpar(hdr, 'GAIN', comment=comment_GAIN)
		if (gain gt 0) then sxaddpar, header, 'GAIN', gain, comment_GAIN, /saveComment





		SXADDPAR, header, 'CDELT1'	,	float(strtrim(sxpar(hdr, 'CDELT1'	, comment=com))), com
		SXADDPAR, header, 'CDELT2'	,	float(strtrim(sxpar(hdr, 'CDELT2'	, comment=com))), com
		SXADDPAR, header, 'LONGPOLE',	float(strtrim(sxpar(hdr, 'LONGPOLE'	, comment=com))), com
		SXADDPAR, header, 'CTYPE1'	,	strtrim(sxpar(hdr, 'CTYPE1'			, comment=com))	, com
		SXADDPAR, header, 'CTYPE2'	,	strtrim(sxpar(hdr, 'CTYPE2'			, comment=com))	, com
		SXADDPAR, header, 'CUNIT1'	,	strtrim(sxpar(hdr, 'CUNIT1'			, comment=com))	, com
		SXADDPAR, header, 'CUNIT2'	,	strtrim(sxpar(hdr, 'CUNIT2'			, comment=com))	, com
		SXADDPAR, header, 'WCS-ORIG',	strtrim(sxpar(hdr, 'WCS-ORIG'		, comment=com))	, com
		SXADDPAR, header, 'RADECSYS',	strtrim(sxpar(hdr, 'RADECSYS'		, comment=com))	, com
		SXADDPAR, header, 'CD1_1'	,	float(strtrim(sxpar(hdr, 'CD1_1'	, comment=com))), com
		SXADDPAR, header, 'CD1_2'	,	float(strtrim(sxpar(hdr, 'CD1_2'	, comment=com))), com
		SXADDPAR, header, 'CD2_1'	,	float(strtrim(sxpar(hdr, 'CD2_1'	, comment=com))), com
		SXADDPAR, header, 'CD2_2'	,	float(strtrim(sxpar(hdr, 'CD2_2'	, comment=com))), com







		if 0 then begin
	
	
			filter 	 = strtrim(sxpar(hdr, 'FILTER01', comment=cfilter))					;	filter
			object	 = strtrim(sxpar(hdr, 'OBJECT', comment=cobject))					;	object
			observer = strtrim(sxpar(hdr, 'OBSERVER', comment=cobserver))				;	observers
			propID	 = strtrim(sxpar(hdr, 'PROP-ID', comment=cpropid))					;	proposal id
			RA2000	 = strtrim(sxpar(hdr, 'RA2000', comment=cra2000))					;	right ascension 2000
			DEC2000	 = strtrim(sxpar(hdr, 'DEC2000', comment=cdec2000))					;	declination 2000
			CRVAL1	 = strtrim(sxpar(hdr, 'CRVAL1', comment=ccrval1))
			CRVAL2	 = strtrim(sxpar(hdr, 'CRVAL2', comment=ccrval2))
			CRPIX1	 = strtrim(sxpar(hdr, 'CRPIX1', comment=ccrpix1))
			CRPIX2	 = strtrim(sxpar(hdr, 'CRPIX2', comment=ccrpix2))
			observat = strtrim(sxpar(hdr, 'OBSERVAT', comment=cobservat))				;	observatory
			telescop = strtrim(sxpar(hdr, 'TELESCOP', comment=ctelescop))				;	telescope
			;effitime = strtrim(sxpar(hdr, 'EFFITIME', comment=ceffitime)	)
			exptime	 = strtrim(sxpar(hdr, 'EXPTIME', comment=cexptime)	)				;	exposure time
			
			;SXADDPAR, header, 'EFFITIME', effitime, ceffitime, /SaveComment
			;ceffitime = 'Effective integration time. [e-/sec]'
			;effitime = 1  
			;SXADDPAR, header, 'EFFITIME', effitime, ceffitime, /SaveComment  
			SXADDPAR, header, 'EXPTIME', exptime, cexptime, /SaveComment 				
			SXADDPAR, header, 'FILTER01', filter, cfilter, /SaveComment                            
			SXADDPAR, header, 'OBJECT', object, cobject, /SaveComment                            
			SXADDPAR, header, 'OBSERVER', observer, cobserver, /SaveComment                            
			SXADDPAR, header, 'PROP-ID', propID, cpropid, /SaveComment                            
			SXADDPAR, header, 'RA2000', RA2000, cra2000, /SaveComment                            
			SXADDPAR, header, 'DEC2000', DEC2000, cdec2000, /SaveComment 
			SXADDPAR, header, 'CRVAL1', float(CRVAL1), ccrval1, /SaveComment 
			SXADDPAR, header, 'CRVAL2', float(CRVAL2), ccrval2, /SaveComment 		                           
			SXADDPAR, header, 'CRPIX1', float(CRPIX1), ccrpix1, /SaveComment 
			SXADDPAR, header, 'CRPIX2', float(CRPIX2), ccrpix2, /SaveComment 
			SXADDPAR, header, 'OBSERVAT', observat, cobservat, /SaveComment                            
			SXADDPAR, header, 'TELESCOP', telescop, ctelescop, /SaveComment 
			SXADDPAR, header, 'NUMEXPS', numEXP, 'Number of exposures used to make final mosaic.', /SaveComment 

		endif

	endif

	;for each ccd do the following
	sxaddhist, '', header															;add history to the header
	sxaddhist, 'Exposure parameters from exposure #' + strtrim(i,1) + '.', header	;add history to the header
	sxaddhist, 'Image = ' + files[ind[i]], header									;add history to the header

	
	exptimes[i]	 = strtrim(sxpar(hdr, 'EXPTIME', comment=cexptime)	)				;exposure time for each separate exposure
	airmasses[i] = strtrim(sxpar(hdr, 'AIRMASS', comment=cairmass)	)				;airmass for each exposure
	obsDates[i]	 = strtrim(sxpar(hdr, 'DATE-OBS', comment=cdateObs)	)				;observation date of each exposure
	UTs[i] 		 = strtrim(sxpar(hdr, 'UT', comment=cut))							;UT of each exposure
	HSTs[i]		 = strtrim(sxpar(hdr, 'HST', comment=chst))							;HST of each exposure
	LSTs[i]		 = strtrim(sxpar(hdr, 'LST', comment=clst))							;LST of each exposure
	MJDs[i]		 = strtrim(sxpar(hdr, 'MJD', comment=cmjd))							;MJD of each exposure


	SXADDPAR, header, 'EXPTIM' + strtrim(i,1), exptimes[i], cexptime, /SaveComment	;EXPTIME for a single exposure                  
	SXADDPAR, header, 'AIRMAS' + strtrim(i,1), airmasses[i], cairmass, /SaveComment	;AIRMASS for a single exposure
	SXADDPAR, header, 'DATEOB' + strtrim(i,1), obsDates[i], cdateobs, /SaveComment	;DATEOBS for a single exposure
	SXADDPAR, header, 'UT' + strtrim(i,1), UTs[i], cut, /SaveComment				;UT 	 for a single exposure
	SXADDPAR, header, 'HST' + strtrim(i,1), HSTs[i], chst, /SaveComment				;HST 	 for a single exposure
	SXADDPAR, header, 'LST' + strtrim(i,1), LSTs[i], clst, /SaveComment				;LST 	 for a single exposure
	SXADDPAR, header, 'MJD' + strtrim(i,1), MJDs[i], cmjd, /SaveComment				;MJD	 for a single exposure



endfor


sxaddhist, '', header										;spacing

if (i gt 1) then begin
	SXADDPAR, header, 'AIRMASS', median(airmasses), 'Median airmass of all exposures.', /SaveComment                 					;Median airmass for all exposures
	SXADDPAR, header, 'UT', jaa_dec2hms(median(hms2dec(uts)), /double, /colon), 'Median HH:MM:SS.S UTC of all exposures.', /SaveComment		;Median UT for all exposures                   
	SXADDPAR, header, 'HST', jaa_dec2hms(median(hms2dec(HSTs)), /double, /colon), 'Median HH:MM:SS.S HST of all exposures.', /SaveComment	;Median HST for all exposures
	SXADDPAR, header, 'LST', jaa_dec2hms(median(hms2dec(lsts)), /double, /colon), 'Median HH:MM:SS.SSS LST of all exposures.', /SaveComment	;Median LST for all exposures
	SXADDPAR, header, 'MJD', median(mjds), '[d] Median Mod. Julian Day of all exposures.', /SaveComment									;Median MJD for all exposures
endif 
if (i eq 1) then begin
	SXADDPAR, header, 'AIRMASS', float(airmasses[0]), 'Median airmass of all exposures.', /SaveComment                 					;Median airmass for all exposures
	SXADDPAR, header, 'UT', jaa_dec2hms((hms2dec(uts[0])), /double, /colon), 'Median HH:MM:SS.S UTC of all exposures.', /SaveComment		;Median UT for all exposures                   
	SXADDPAR, header, 'HST', jaa_dec2hms((hms2dec(HSTs[0])), /double, /colon), 'Median HH:MM:SS.S HST of all exposures.', /SaveComment	;Median HST for all exposures
	SXADDPAR, header, 'LST', jaa_dec2hms((hms2dec(lsts[0])), /double, /colon), 'Median HH:MM:SS.SSS LST of all exposures.', /SaveComment	;Median LST for all exposures
	SXADDPAR, header, 'MJD', float(mjds[0]), '[d] Median Mod. Julian Day of all exposures.', /SaveComment									;Median MJD for all exposures
endif



modfits, mosaic, 0, header															;put modified head back into image







end


















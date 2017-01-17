;+
;
; NAME
;	sCamPipe
;
; PURPOSE
;	Reduce raw Suprime-Cam imaging data into a final co-added mosaic.
;
; CALLING SEQUENCE
;	sCamPipe, filter=filter, flatType=flatType, object=object, pre=pre, sm=sm, $
;		outfile=outfile, batch=batch, redoastrom=redoastrom, noastrom=noastrom, $
;		catalog=catalog, nthreads=nthreads, nobatch=nobatch, noinherit=noinherit, $
;		distcorr=distcorr
;
; KEYWORDS:
;	FILTER: name of the filter used to take imaging (e.g., W-S-G+, W-S-R+, W-S-I+)
;	FLATTYPE: type of flatfielding (e.g., DOMEFLAT, SKYFLAT, SUPERSKYFLAT) 
;	OBJECT: name of the targeted object (e.g., NGC 821)
;	PRE: prefix of the files to reduce. Use this keyword when starting the pipeline
;		at an intermediate step and working on a non-standard file prefix. 
;		For instance, if you wanted to mask the AG probe of the non-distortion-corrected
;		files (prefix 'c'), then set pre='c'. Otherwise, the code will operage on files
;		with the prefix 'g' (the distortion corrrected files).
;	OUTFILE: name of the final mosaiced image.
;	BATCH: set this keyword to ...		
;	SM:
;	REDOASTROM:
;	NOASTROM: 
;	CATALOG:
;	NTHREADS:
;	NOBATCH
;	NOINHERIT:
;	DISTCORR:
;
; PROCESS:
;	(1) RENAME FILES
;	 		(SUPA* -> H*)
;			Change file names to reflect which chip (of 10) they belong to.
;	(2) SUBTRACT OVERSCAN
;			(-> To_RH*)
;			Subtract off bias as measured from the overscan regions on each chip.
;			Create and subtract a superbias created from overscan subtracted bias frames
;			Trim off unused portions of each image. 
;	(3) MAKE FLAT-FIELD FRAMES	
;			(obj_mflat*)
;			Calculate the gain and readnoise from the bias and domeflat frames.
;			Create sigma-clipped median combined flat-field frames.
;			Mask the auto-guider (AG) probe.
;	(4) DIVIDE BY FLATS
;			(-> f*)
;			Divide by the flat-field frames.
;	(5) REMOVE COSMIC RAYS
;			(-> z*)
;			Remove cosmic rays using L.A. Cosmic (http://www.astro.yale.edu/dokkum/lacosmic/)
;	(6) MULTIPLICATIVELY SCALE CCDS TO COMMON SENSITITIVY
;			(-> c*)
;			Correct for chip-to-chip QE variations, and any time-variable gain by scaling
;			invidividual chips. Flux levels across chip boundaries are minimized using
;			a Gauss-Newton minimizer.
;	(7) DISTORTION CORRECTION
;			(-> g*)
;			Correct for the predominantly radial distortion of S-Cam.
;	(8) MASK AUTO-GUIDER PROBE
;			(-> A*)
;			Mask out the auto-guider probe on each frame.
;	(9) SUBTRACT SKY
;			(-> s*)
;			Currently a simple median of medians.
;	(10) CALIBRATE ASTROMETRY
;			(-> as*)
;			Calibrate astrometry using SDSS or USNO.
;	(11) REPROJECT IMAGES
;			(-> projdir/hdu*)
;			Project images to the tangent-plane.
;	(12) WEIGHT EACH PIXEL
;			Determine a weight for each pixel using the exposure time, the drizzled area, and a bad pixel mask.
;	(13) CO-ADD INTO MOSAIC
;			Create a final median co-added mosaic image using Montage (not currently using weighting).
;
;
; RESULTS:
;	Reduced and mosaiced image.
;
;-----------------------------------------
;
; REQUIRED ROUTINES:
;
;
; MODIFICATION HISTORY:
;	Written by Jacob Arnold, UC-Santa Cruz, October 22, 2012
;
;
;-
;-------------------------------------------------------------------------------------------

pro sCamPipe, filter=filter, flatType=flatType, object=object, pre=pre, outfile=outfile, $
		batch=batch, sm=sm, redoastrom=redoastrom, noastrom=noastrom, catalog=catalog, $
		nthreads=nthreads, nobatch=nobatch, noinherit=noinherit, distcorr=distcorr, $
		coaddType=coaddType
	COMPILE_OPT idl2, HIDDEN

;	resolve_routine, 'supCamData__define', /compile
;	resolve_all;, /quiet
	;save, /routines, filename='supCamPipe_2012Oct25.sav'
		
	if (n_elements(nthreads) eq 0) then nthreads = 0 ;default to 6 threads when running in parallel
	astrom = 1 - n_elements(noastrom)

	data = obj_new('supCamData', filter=filter, object=object, $ ;create a new Suprime-Cam data object
		flatType=flatType, batch=batch, coAddType=coAddType, $
		outFile=outFile)

	done = 'no'
	k = 1
	WHILE (done ne 'yes') DO BEGIN
		print, '----------------------------------------------'
		print, 'Choose an option.'
		print, '----------------------------------------------'
		print
		options = [	$
			'Rename Files', $
			'Subtract Overscan', $
			'Make Flats', $
			'Flatfield Data', $
			'Remove cosmic rays using L.A. Cosmic.', $
			'Multiplicatively Scale Chips (correct for a time dependent gain)', $
			'Correct for Distortion and Atmospheric Dispersion', $
			'Mask the Auto-Guiding Probe', $
			'Subtract Sky']
		options = [options, 'Calibrate Astrometry', $
							'Scale exposures to common flux level.', $
							'Reproject Images', $
							'Determine Image weighting']
		options = [options, 'Combine Frames']
		options = [options, 'Perform All', $
							'Perform a Subset', $
							'Exit'	]
	
		for i=1,n_elements(options) do print, i, options[i-1], format='(I-10,A)'
		print
		print, '----------------------------------------------'
		
		choices = fix(choose('Choice(s):',[k]))
		k = choices[0]
	
	
		i = 0
		WHILE 1 DO BEGIN	
	
			choice = choices[i]
	
			case choice of
	
				1:	data->rename_files	;Rename Files SUPA* 			-> H*
	
				2:	data->overscan_subtract	;Subtract Overscan H* 			-> To_RH*
	
				3:	data->make_flats, sm=sm	;Create Flats				-> obj_mflat*
	
				4:	data->flatfield_data	;Divide by Flats			-> f*
				
				5:	data->clean_cosmics	;Remove Cosmic Rays			-> z*
					
				6:	data->scale_chips	;Scale Chips to common level		-> c*
	
				7:	data->distortion_correction	;Correct Distortion		-> g*

				8:	data->mask_AGprobe	;Mask AG Probe				-> A*
				
				9:	data->subtract_sky	;Subtract Sky				-> s*

				10:	data->astrom, redoastrom=redoastrom, catalog=catalog, $	;Calibrate Astrometry 	-> as*
						nthreads=nthreads, nobatch=nobatch, noinherit=noinherit, $
						distcorr=distcorr			

				11:	data->scale_exps	;Scale exposures to same level		-> x* 
	
				12:	data->reproject		;Reproject Images			-> projdir/hdu*
	
				13:	data->weight_area	;Weight Each exposure			-> projdir/hdu/area/hdu*area.fits	
					
				14:	BEGIN
						data->coAdd	;Co-add exposures
						done = 'yes'
					END
									
				15:	BEGIN
						choices = [0, indgen(12)+1]
						done = 'yes'
					END
		
				16:	BEGIN
						choices = [0, choose('Which options?', [1,2,3,4,5], /preserve)]
						done = 'yes'
					END
			
				17:	done = 'yes'
				
				else:	done = 'yes'
	
			ENDCASE
		
			if (i eq (n_elements(choices)-1)) then break
			i++
		ENDWHILE
	
	++k
	ENDWHILE

	obj_destroy, data
	heap_gc


end

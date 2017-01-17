


	objects = runse('tmos.fits', zp=29., exptime=1d, numexps=1d, /batch, saturate=2944d, oexptime=74., brightStars=bStars)
	d0 = seSelection(objects, {magcut:[18d, 22d], CLASS_STAR:0.0d, fwhmcut:[90d, 70d], apercut:[90d, 90d], sncut:[40d, 80d]}, /simple, /batch)


	objects = runse('tmos1.fits', zp=29., exptime=1d, numexps=1d, /batch, saturate=2944d, oexptime=74., brightStars=bStars)
	d1 = seSelection(objects, {magcut:[18d, 22d], CLASS_STAR:0.0d, fwhmcut:[90d, 70d], apercut:[90d, 90d], sncut:[40d, 80d]}, /simple, /batch)


	match = match_2d(d0.x_world, d0.y_world, d1.x_world, d1.y_world, 5d/3600d, match_distance=mindist)
	ind = where( (match ne -1) and (mindist*3600d lt 2d), ngood)
	if (nGood eq 0) then stop
	
	
	c = d0[ind]
	d = d1[match[ind]]
	
	
	ratio = median(d.flux_aper2/c.flux_aper2)
	
	
	


stop



end
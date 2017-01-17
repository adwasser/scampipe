

					
function maskgrow, subImg0, ind0, gr, value, newind=newind, mask=mask, xpoz=xpoz, ypoz=ypoz

	if (n_elements(xpoz) eq 0) then xpoz = 1 else xpoz = 0
	if (n_elements(ypoz) eq 0) then ypoz = 1 else ypoz = 0
	if (n_elements(gr) eq 1) then gr = [gr,gr]
	grx = gr[0]
	gry = gr[1]
	
	mask = finite(subImg0)*0																	;make an empty integer array the same size as the input image
	mask[ind0] = 1																				;set masked pixels = 1
	
	supermask = intarr(n_elements(subImg0[*,0]) + 2*grx, n_elements(subImg0[0,*]) + 2*gry)	;create a larger image to work within
	supermask[grx,gry] = mask																		;place the original mask inside the supermask
	maskblank = supermask																		;create a dummy mask to use
	;for gry=-gr,gr do for grx=-gr,gr do print, minmax((mask1 = mask1 + shift(maskblank,grx,gry)))
	for y=-gry*ypoz,gry do begin
		for x=-grx*xpoz,grx do begin
			supermask = supermask + shift(maskblank,x,y)				;step left-right and up-down and shift the dummy mask around, effectively growing the masked regions
		endfor
	endfor
	;display, mask1, min=0, max=1, top=254 
	;mask1 = mask1 + shift(mask1,-1*gr,0) + shift(mask1,-1*gr,1*gr) + shift(mask1,0,1*gr) + shift(mask1,1*gr,1*gr) + shift(mask1,1*gr,0) + shift(mask1,1*gr,-1*gr) + shift(mask1,0,-1*gr) + shift(mask1,-1*gr,-1*gr)
	mask = supermask[grx:-grx-1,gry:-gry-1]															;trim off the unneeded parts
	mask = mask + 1																				;add 1 to the whole mask
	mask[(newind = where(mask gt 1))] = 0														;set any masked pixels to 0, unmasked pixels are now 1
	
	result = subImg0*mask																		;multiply the input image by the mask
	result[where(mask eq 0)] = value															;give the masked pixels the specified value
	;display, subImg0*mask, min=mode, max=5000, top=254
	
	return, result																				;return the masked image

end

function getsubimgs, file, img=img, hdr=hdr, xy=xy
COMPILE_OPT idl2, HIDDEN 

img = mrdfits(file, 0, hdr, /silent)																;primary chip

BZERO = sxpar(hdr, 'BZERO')
xy = {lowx:fltarr(4), highx:fltarr(4), lowy:fltarr(4), highy:fltarr(4)}

for region=1,4 do begin																		;isolate each subchip
	lowx = 	sxpar(hdr, 'S_EFMN' + strtrim(region,2) + '1') - 1
	highx = sxpar(hdr, 'S_EFMX' + strtrim(region,2) + '1') - 1
	lowy = 	sxpar(hdr, 'S_EFMN' + strtrim(region,2) + '2') - 1
	highy = sxpar(hdr, 'S_EFMX' + strtrim(region,2) + '2') - 1

	xy.lowx[region-1] = lowx
	xy.highx[region-1] = highx
	xy.lowy[region-1] = lowy
	xy.highy[region-1] = highy

	dsubImg = create_struct('num'+strtrim(region,2), img[lowx:highx,lowy:highy] + BZERO)
	if (region eq 1) then subImg = dsubImg else subImg = struct_addtags(subImg, dsubImg)
endfor

return, subImg

end

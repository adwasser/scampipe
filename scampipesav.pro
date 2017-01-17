
pro scampipesav
	procList = ['scampipe','supCamData__define']
	foreach progName, procList do resolve_routine, progName, /compile
	funcList = ['arr_struct']
	foreach progName, funcList do resolve_routine, progName, /compile, /is_function
	;resolve_routine, 'scampipe', /compile
	;resolve_routine, 'supCamData__define', /compile
	resolve_all
	time = (strsplit(systime(),' ',/extract))[1:*]
	time = '_'+time[3]+time[0]+time[1]+'_'+time[2]
	time = strsplit(time,':',/extract)
	time = time[0]+'.'+time[1]+'.'+time[2]
	dir = file_search('savFiles',count=nDir)
	if (nDir eq 0) then spawn, 'mkdir savFiles'
	save, /routines, filename='savFiles/sCamPipe'+time+'.sav'
end
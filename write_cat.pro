;Quick program to take in the catalog structure file from scampipe
;and output to .csv files
;
;
;
;

pro write_cat,cat,file_suffix=file_suffix

tags = strsplit(tag_names(cat),' ',/extract)

for i=0,n_elements(tags)-1 do begin
   x = cat.(i).x
   y = cat.(i).y
   ra = cat.(i).ra
   dec = cat.(i).dec

   file = cat.(i).file
   this_file = file[0]
   split = strsplit(this_file,'_',/extract)

   chip_file = strsplit(split[2],'.',/extract)
   chip_name = chip_file[0]

   filename=split[1]+'_'+chip_name+'.csv'

   openw,unit,filename,/get_lun
   for j=0,n_elements(x)-1 do begin
      printf,unit,strtrim(x[j],2),',',strtrim(y[j],2),',',$
             strtrim(ra[j],2),',',strtrim(dec[j],2)
   endfor
   close,unit
   free_lun,unit

endfor

return
end

begin

f1 = addfile("/mnt/hgfs/iap_down0711/t2m.nc","r")
print(f1)
f2 = addfile("/mnt/hgfs/iap_down0711/cru_ts4.07.1901.2022.tmp.dat.nc","r")

t2m = f1->t2m
low_lat = f1->lat
low_lon1 = f1->lon
low_lon = low_lon1 - 179.75
printVarSummary(low_lon)

tmp = f2->tmp
printVarSummary(tmp)

high_lat = f2->lat
high_lon = f2->lon
printVarSummary(high_lon)

t2m_new = linint2_Wrap(low_lon, low_lat, t2m, True, high_lon, high_lat, 0)
printVarSummary(t2m_new)


t2m_2 = t2m_new - 273.15
copy_VarCoords(t2m_new, t2m_2)

out = addfile("/mnt/hgfs/iap_down0711/t2m_new.nc","c")
out->t2m = t2m_2

; out->LON = high_lon


; wks_type = "png"
; wks_type@wkWidth = 2500
; wks_type@wkHeight = 2500

; wks = gsn_open_wks(wks_type,"DJF_low")

; res = True


; plot = gsn_csm_contour_map_ce(wks, t2m(0,:,:), res)

; wks = gsn_open_wks(wks_type,"DJF_high")

; res = True


; plot = gsn_csm_contour_map_ce(wks, t2m_new(0,:,:), res)
end
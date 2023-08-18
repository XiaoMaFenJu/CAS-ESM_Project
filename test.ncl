begin

f = addfile("/mnt/hgfs/iap_down0711/B20TRC5X_IAPDGVM-colm-2001-01.nc","r")
print(f)

t2m = f->prl
printVarSummary(t2m)

lat = f->lat
lon = f->lon

lat@units = "degrees_north"
lon@units = "degrees_east"

t2m!1 = "lat"
t2m!2 = "lon"

t2m&lat = lat
t2m&lon = lon

; t2m@_FillValue = -9999.0

printVarSummary(t2m)
; printVarSummary(lat)

; f = addfile("/home/xiaoma/nc/data/ncep.air.mon.mean.nc","r")
; print(f)

wks_type = "png"
wks_type@wkWidth = 2500
wks_type@wkHeight = 2500

wks = gsn_open_wks(wks_type,"DJFtt")

res = True
res@cnFillOn = True
res@cnLineColor = "black"
; res@cnLineThicknessF = 10
res@cnLevelSelectionMode = "Explicitlevels"
res@cnLevels = (/0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10/)


plot = gsn_csm_contour_map_ce(wks, t2m(0,:,:), res)

end
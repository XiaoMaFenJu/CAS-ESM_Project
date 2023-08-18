begin

f1 = addfile("/mnt/hgfs/iap_down0711/t2m_CRU.nc","r")
f2 = addfile("/mnt/hgfs/iap_down0711/t2m_new.nc","r")
; print(f1)
; print(f2)
lon1 = f1->lon
lat = f1->lat
lon = lon1
lon(0:359) = lon1(360:719)
lon(360:719) = lon1(0:359)

t2m_C1 = f1->t2m(420:779,{50:90},:)
t2m_M2 = f2->t2m(420:779,{50:90},:)
; printVarSummary(t2(180:779,:,:)m_C(419:780,:,:))
t2m_M1 = t2m_M2
t2m_M1(:,:,0:359) = t2m_M2(:,:,360:719)
t2m_M1(:,:,360:719) = t2m_M2(:,:,0:359)

t2m_M1(:,{85:90},:) = -9999.0
t2m_M1!0 = "time"
t2m_M1!1 = "lat"
t2m_M1!2 = "lon"
t2m_C1!0 = "time"
t2m_C1!1 = "lat"
t2m_C1!2 = "lon"

t2m_M = t2m_M1(lat|:,lon|:,time|:)
t2m_C = t2m_C1(lat|:,lon|:,time|:)

corc = escorc(t2m_C, t2m_M)
copy_VarCoords(t2m_C,corc)
printVarSummary(corc)
prob = rtest(corc, 30*12, 0)
; print(prob(40:60,:))
copy_VarCoords(corc,prob)

print(max(prob))
print(min(prob))
; siglvl = 0.05
; yes_no = where(prob.lt.siglvl, True, False)
; p = where(prob.gt.0.05, prob@_FillValue, prob)
; copy_VarCoords(corc,yes_no)
; printVarSummary(p)

wks_type = "png"
wks_type@wkWidth = 5000
wks_type@wkHeight = 5000

wks = gsn_open_wks(wks_type,"r")

res = True
; res@vpWidthF=0.45
; res@vpHeightF=0.45
res@gsnFrame=False
res@gsnDraw=False
res@cnFillOn = True
res@cnLinesOn = False
; res@lbLabelBarOn = False
res@tmBorderThicknessF = 15
res@cnLineThicknessF = 2
res@cnFillPalette = "MPL_GnBu"
res@gsnPolar = "NH"
res@gsnLeftStringFontHeightF = 0.013
res@gsnRightStringFontHeightF = 0.013
res@mpMinLatF = 50
res@mpShapeMode             = "FreeAspect"
res@gsnMaximize  = False

res@gsnLeftString = "Corelation Coefficient"
res@gsnRightString = "2m Temperature"
res@gsnLeftStringFontHeightF = 0.018
res@gsnRightStringFontHeightF = 0.018
; res@gsnRightString = "~S~o~N~C"
res@gsnRightStringOrthogonalPosF = 0.04
res@gsnLeftStringOrthogonalPosF = 0.05

res@pmLabelBarWidthF = 0.08;colorbar width
res@pmLabelBarHeightF = 0.4;colorbar height
res@lbOrientation = "Vertical"
res@lbBoxEndCapStyle = "TriangleBothEnds"
res@pmLabelBarOrthogonalPosF = 0.08
res@pmLabelBarParallelPosF = 0.5

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.8
res@cnMaxLevelValF = 0.98
res@cnLevelSpacingF = 0.03
res@cnFillDrawOrder = "PreDraw"

plot = gsn_csm_contour_map_polar(wks, corc, res)

resp = True
resp@gsLineThicknessF = 10
y = new(720,"float")
y = 50.1
; print(y)
pl1 = gsn_add_polyline(wks,plot,lon,y,resp)

rest = True
rest@gsnFrame=False
rest@gsnDraw=False
rest@cnLevelSelectionMode="ExplicitLevels"
rest@cnLevels=(/0,0.01/)
rest@cnFillPalette="GMT_gray"
rest@cnFillColors         = (/1,1,1/)
rest@cnFillPattern        =17
rest@cnFillScaleF=0.5
rest@cnInfoLabelOn=False
rest@lbLabelBarOn=False
plot1 = gsn_csm_contour(wks,prob,rest)

; printVarSummary(prob)
; draw(plot)
draw(plot1)
; overlay(plot, plot1)
frame(wks)
; 
end


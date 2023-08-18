begin

f1 = addfile("/mnt/hgfs/iap_down0711/t2m_CRU.nc","r")
f2 = addfile("/mnt/hgfs/iap_down0711/t2m_new.nc","r")
; print(f1)
; print(f2)
t2m_C = f1->t2m(420:779,{50:90},:) ;420:779 后三十年
t2m_M1 = f2->t2m(420:779,{50:90},:)
lon1 = f1->lon
lon = lon1
lon(0:359) = lon1(360:719)
lon(360:719) = lon1(0:359)

t2m_M1(:,{85:90},:) = -9999.0
t2m_M = t2m_M1
t2m_M(:,:,0:359) = t2m_M1(:,:,360:719)
t2m_M(:,:,360:719) = t2m_M1(:,:,0:359)
copy_VarMeta(t2m_M1,t2m_M)

lat = f1->lat({50:90})

t2m_C_DJF = month_to_season(t2m_C, "DJF")
t2m_C_JJA = month_to_season(t2m_C, "JJA")
t2m_M_DJF = month_to_season(t2m_M, "DJF")
t2m_M_JJA = month_to_season(t2m_M, "JJA")

t2m_C_MAM = month_to_season(t2m_C, "MAM")
t2m_C_SON = month_to_season(t2m_C, "SON")
t2m_M_MAM = month_to_season(t2m_M, "MAM")
t2m_M_SON = month_to_season(t2m_M, "SON")

printVarSummary(t2m_C_DJF)

corc = escorc_n(t2m_C_SON, t2m_M_SON,0,0)
printVarSummary(corc)
copy_VarCoords(t2m_C_DJF(0,:,:),corc)
prob = rtest(corc, 30, 0)
printVarSummary(prob)
copy_VarCoords(corc,prob)

; p = prob
; p = where(prob.gt.0.05, prob@_FillValue, prob)
; copy_VarCoords(t2m_C,corc)
; copy_VarCoords(t2m_C,prob)

wks_type = "png"
wks_type@wkWidth = 5000
wks_type@wkHeight = 5000

wks = gsn_open_wks(wks_type,"r_SON")


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
res@cnFillPalette = "MPL_RdBu"
res@gsnPolar = "NH"
res@gsnLeftStringFontHeightF = 0.013
res@gsnRightStringFontHeightF = 0.013
res@mpMinLatF = 50
res@mpShapeMode = "FreeAspect"
res@gsnMaximize  = False

res@gsnLeftString = "Corelation Coefficient"
res@gsnRightString = "2m Temperature (SON)"
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
res@cnMinLevelValF = -0.8
res@cnMaxLevelValF = 0.8
res@cnLevelSpacingF = 0.2

res@trGridType = "TriangularMesh"
res@mpGridLineThicknessF = 10

res@lbBoxLineThicknessF = 10

res@mpOceanFillColor = "white"
res@mpLandFillColor = -1
res@cnFillDrawOrder = "PreDraw"
res@mpFillDrawOrder = "Draw"

; res@mpAreaMaskingOn=True
; res@mpDataBaseVersion = "Ncarg4_0" 
; res@mpDataSetName = "Earth..4"
; res@mpSpecifiedFillColors = -1
; res@mpMaskAreaSpecifiers=(/"ocean"/)
; res@mpFillDrawOrder = "PostDraw"

plot = gsn_csm_contour_map_polar(wks, corc, res)

rest = True
rest@gsnFrame=False
rest@gsnDraw=False
rest@cnFillOn = True
rest@cnLinesOn = False
rest@cnLineLabelsOn = False
rest@cnLevelSelectionMode="ExplicitLevels"
rest@cnFillPalette="GMT_gray"
rest@cnLevels=(/0.05/)
rest@cnFillColors         = (/1,-1/)
rest@cnFillPattern        =17
rest@cnFillScaleF = 1.5
rest@cnFillDotSizeF = 0.004
rest@cnInfoLabelOn=False
rest@lbLabelBarOn=False
rest@trGridType = "TriangularMesh"
rest@cnFillDrawOrder = "PreDraw"


; prob({50.25},0) = 0.00001
; 
plot1 = gsn_csm_contour(wks,prob,rest)

resp = True
resp@gsLineThicknessF = 15
resp@tfPolyDrawOrder = "PostDraw"
y = new(720,"float")
y = 50.1
; print(y)
pl1 = gsn_add_polyline(wks,plot,lon,y,resp)


overlay(plot, plot1)
draw(plot)
frame(wks)

print(max(prob))
print(min(prob))
end


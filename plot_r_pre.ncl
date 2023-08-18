begin

f1 = addfile("/mnt/hgfs/iap_down0711/pre_CRU.nc","r")
f2 = addfile("/mnt/hgfs/iap_down0711/pre_new.nc","r")
; print(f1)
; print(f2)
lon1 = f1->lon
lat = f1->lat
lon = lon1
lon(0:359) = lon1(360:719)
lon(360:719) = lon1(0:359)

t2m_C1 = f1->pre(420:779,{50:90},:)
t2m_M2 = f2->pre(420:779,{50:90},:)
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
prob = rtest(corc, 30*12, 0)
copy_VarMeta(corc,prob)
printVarSummary(prob)
; p = prob
; p = where(prob.gt.0.05, prob@_FillValue, prob)
; copy_VarCoords(t2m_C,corc)
; copy_VarCoords(t2m_C,prob)

wks_type = "png"
wks_type@wkWidth = 5000
wks_type@wkHeight = 5000

wks = gsn_open_wks(wks_type,"r_pre")

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
res@mpShapeMode = "FreeAspect"
res@gsnMaximize  = False

res@gsnLeftString = "Corelation Coefficient"
res@gsnRightString = "Precipitation"
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
res@cnMinLevelValF = -0.4
res@cnMaxLevelValF = 0.8
res@cnLevelSpacingF = 0.2

res@trGridType = "TriangularMesh"
res@mpGridLineThicknessF = 5

res@lbBoxLineThicknessF = 5

; res@mpAreaMaskingOn=True
; res@mpDataBaseVersion = "Ncarg4_0" 
; res@mpDataSetName = "Earth..4"
; res@mpSpecifiedFillColors = -1
; res@mpMaskAreaSpecifiers=(/"ocean"/)
; res@mpFillDrawOrder = "PostDraw"

plot = gsn_csm_contour_map_polar(wks, corc, res)

resp = True
resp@gsLineThicknessF = 5
y = new(720,"float")
y = 50.1
; print(y)
pl1 = gsn_add_polyline(wks,plot,lon,y,resp)

rest = True
rest@gsnFrame=False
rest@gsnDraw=False
rest@cnFillOn = True
rest@cnLinesOn = False
rest@cnLineLabelsOn = False
rest@cnLevelSelectionMode="ExplicitLevels"
rest@cnLevels=(/0.05/)
rest@cnFillPalette="GMT_gray"
rest@cnFillColors         = (/1,-1/)
rest@cnFillPattern        =17
; rest@cnFillScaleF= 
rest@cnFillDotSizeF = 0.004
rest@cnInfoLabelOn=False
rest@lbLabelBarOn=False
rest@trGridType = "TriangularMesh"


; ; printVarSummary(prob)
; ; lsm  = landsea_mask(prob,prob&lat,prob&lon)
; ; prob_land = mask(prob, lsm.eq.1 ,False)
; a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
; lsdata1 = a->LSMASK
; ; printVarSummary(lsdata1)
; lsdata = lsdata1
; ; lsdata(:,0:179) = lsdata1(:,180:359)
; ; lsdata(:,180:359) = lsdata1(:,0:179)
; lsm  = landsea_mask(lsdata,prob&lat,prob&lon)
; prob_land1 = mask(prob, lsm.eq.1 ,False)
; prob_land = prob_land1
; ; prob_land(:,0:359) = prob_land1(:,360:719)
; ; prob_land(:,360:719) = prob_land1(:,0:359)

plot1 = gsn_csm_contour_map_ce(wks,prob,rest)


; draw(plot1)
overlay(plot, plot1)
draw(plot)
frame(wks)

end


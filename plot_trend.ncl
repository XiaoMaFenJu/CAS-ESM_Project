begin
f1 = addfile("/mnt/hgfs/iap_down0711/t2m_CRU.nc","r")
f2 = addfile("/mnt/hgfs/iap_down0711/t2m_new.nc","r")
; print(f1)
; print(f2)

t2m_C = f1->t2m(180:779,{50:90},:)
t2m_M1 = f2->t2m(180:779,{50:90},:)

t2m_M1(:,{85:90},:) = -9999.0
t2m_M = t2m_M1
t2m_M(:,:,0:359) = t2m_M1(:,:,360:719)
t2m_M(:,:,360:719) = t2m_M1(:,:,0:359)
printVarSummary(t2m_M)

lon1 = f1->lon
lat = f1->lat
lon = lon1
lon(0:359) = lon1(360:719)
lon(360:719) = lon1(0:359)
; print(lon)

t2m_C_annual = month_to_annual(t2m_C, 1);0是求和
t2m_M_annual = month_to_annual(t2m_M, 1);0是求和 

t2m_C_DJF = month_to_season(t2m_C, "DJF")
t2m_C_JJA = month_to_season(t2m_C, "JJA")
t2m_M_DJF = month_to_season(t2m_M, "DJF")
t2m_M_JJA = month_to_season(t2m_M, "JJA")

t2m_C_MAM = month_to_season(t2m_C, "MAM")
t2m_C_SON = month_to_season(t2m_C, "SON")
t2m_M_MAM = month_to_season(t2m_M, "MAM")
t2m_M_SON = month_to_season(t2m_M, "SON")

; x = ispan(0, 49, 1)
; rc = regline(x, T_annual(0,:))
; y = rc*x + rc@yintercept
; printVarSummary(y)
p_M1 = trend_manken(t2m_M_annual,False,0)
p_M2 = trend_manken(t2m_M_DJF,False,0)
p_M3 = trend_manken(t2m_M_MAM,False,0)
p_M4 = trend_manken(t2m_M_JJA,False,0)
p_M5 = trend_manken(t2m_M_SON,False,0)

p_C1 = trend_manken(t2m_C_annual,False,0)
p_C2 = trend_manken(t2m_C_DJF,False,0)
p_C3 = trend_manken(t2m_C_MAM,False,0)
p_C4 = trend_manken(t2m_C_JJA,False,0)
p_C5 = trend_manken(t2m_C_SON,False,0)

printVarSummary(p_M1)

alpha1 = 1 -p_M1(0,:,:)
alpha2 = 1 -p_M2(0,:,:)
alpha3 = 1 -p_M3(0,:,:)
alpha4 = 1 -p_M4(0,:,:)
alpha5 = 1 -p_M5(0,:,:)
alpha6 = 1 -p_C1(0,:,:)
alpha7 = 1 -p_C2(0,:,:)
alpha8 = 1 -p_C3(0,:,:)
alpha9 = 1 -p_C4(0,:,:)
alpha10 = 1 -p_C5(0,:,:)

copy_VarCoords(t2m_C_annual(0,:,:),alpha1)
copy_VarCoords(t2m_C_annual(0,:,:),alpha2)
copy_VarCoords(t2m_C_annual(0,:,:),alpha3)
copy_VarCoords(t2m_C_annual(0,:,:),alpha4)
copy_VarCoords(t2m_C_annual(0,:,:),alpha5)
copy_VarCoords(t2m_C_annual(0,:,:),alpha6)
copy_VarCoords(t2m_C_annual(0,:,:),alpha7)
copy_VarCoords(t2m_C_annual(0,:,:),alpha8)
copy_VarCoords(t2m_C_annual(0,:,:),alpha9)
copy_VarCoords(t2m_C_annual(0,:,:),alpha10)

trend1 = p_M1(1,:,:)
trend2 = p_M2(1,:,:)
trend3 = p_M3(1,:,:)
trend4 = p_M4(1,:,:)
trend5 = p_M5(1,:,:)
trend6 = p_C1(1,:,:)
trend7 = p_C2(1,:,:)
trend8 = p_C3(1,:,:)
trend9 = p_C4(1,:,:)
trend10 = p_C5(1,:,:)

trend_delta1 = trend1 - trend6
trend_delta2 = trend2 - trend7
trend_delta3 = trend3 - trend8
trend_delta4 = trend4 - trend9
trend_delta5 = trend5 - trend10

copy_VarCoords(t2m_C_annual(0,:,:),trend1)
copy_VarCoords(t2m_C_annual(0,:,:),trend2)
copy_VarCoords(t2m_C_annual(0,:,:),trend3)
copy_VarCoords(t2m_C_annual(0,:,:),trend4)
copy_VarCoords(t2m_C_annual(0,:,:),trend5)
copy_VarCoords(t2m_C_annual(0,:,:),trend6)
copy_VarCoords(t2m_C_annual(0,:,:),trend7)
copy_VarCoords(t2m_C_annual(0,:,:),trend8)
copy_VarCoords(t2m_C_annual(0,:,:),trend9)
copy_VarCoords(t2m_C_annual(0,:,:),trend10)

copy_VarCoords(t2m_C_annual(0,:,:),trend_delta1)
copy_VarCoords(t2m_C_annual(0,:,:),trend_delta2)
copy_VarCoords(t2m_C_annual(0,:,:),trend_delta3)
copy_VarCoords(t2m_C_annual(0,:,:),trend_delta4)
copy_VarCoords(t2m_C_annual(0,:,:),trend_delta5)

printVarSummary(trend1)
print(dim_avg_n(dim_avg_n(trend1,0), 0))
print(dim_avg_n(dim_avg_n(trend6,0), 0))

wks_type = "png"
wks_type@wkWidth = 5000
wks_type@wkHeight = 5000

wks = gsn_open_wks(wks_type,"t2m_space_trend_season_SON")

plot1 =  new(2,graphic)

res = True
res@vpWidthF=0.31
res@vpHeightF=0.3
res@gsnFrame=False
res@gsnDraw=False
res@cnFillOn = True
res@lbLabelBarOn = False
res@tmBorderThicknessF = 15
res@cnLineThicknessF = 2
res@cnFillPalette = "BlueRed"
res@gsnPolar = "NH"
res@gsnLeftStringFontHeightF = 0.013
res@gsnRightStringFontHeightF = 0.013
res@mpMinLatF = 50
res@mpShapeMode             = "FreeAspect"
res@gsnMaximize  = False
res@mpGridLineThicknessF = 5
; res@cnLevelSelectionMode = "Explicitlevels"
; res@cnLevels = (/-0.05,-0.035,-0.02,-0.01,0,0.01,0.02,0.035,0.05/)
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.12
res@cnMaxLevelValF = 0.12
res@cnLevelSpacingF = 0.015
res@lbBoxLineThicknessF = 5
res@cnLinesOn = False
res@cnLineLabelsOn = False
; res@gsnAddCyclic = False
; res@trGridType = "TriangularMesh"
res@mpOceanFillColor = "white"
res@mpLandFillColor = -1
res@cnFillDrawOrder = "PreDraw"
res@mpFillDrawOrder = "Draw"

res1 = res
res1@gsnLeftString = "OBS(SON)"
res1@gsnRightString = "~S~o~N~C/a"
res1@gsnCenterString = "~S~"
res1@gsnCenterStringOrthogonalPosF = 0.122
; res1@vpXF = 0.01
; res1@vpYF = 0.4

res2 = res
res2@gsnLeftString = "CAS-ESM(SON)"
res2@gsnRightString = "~S~o~N~C/a"
res2@gsnCenterString = "~S~"
res2@gsnCenterStringOrthogonalPosF = 0.122
; res2@vpXF = 0.34
; res2@vpYF = 0.4

res3 = res
res3@gsnLeftString = "Bias(b-a)"
res3@gsnRightString = "~S~o~N~C/a"
; res3@vpWidthF=0.13
; res3@vpHeightF=0.16
res3@lbLabelBarOn = True
; res3@lbOrientation = "Vertical"
res3@pmLabelBarWidthF = 0.3;colorbar width
res3@pmLabelBarHeightF = 0.066;colorbar height
res3@lbLabelFontHeightF = 0.013
res3@gsnLeftStringFontHeightF = 0.013
res3@gsnRightStringFontHeightF = 0.013
res3@gsnMaximize  = True

; res3@cnLevelSelectionMode = "ManualLevels"
; res3@cnMinLevelValF = -1
; res3@cnMaxLevelValF = 1
; res3@cnLevelSpacingF = 0.2
; delete(res3@cnLevels)
; res3@cnLevelSelectionMode = "Explicitlevels"
; res3@cnLevels = (/-0.05,-0.035,-0.02,-0.01,0,0.01,0.02,0.035,0.05/)
res3@cnLevelSelectionMode = "ManualLevels"
res3@cnMinLevelValF = -0.05
res3@cnMaxLevelValF = 0.05
res3@cnLevelSpacingF = 0.01

plot1(0) = gsn_csm_contour_map_polar(wks, trend5, res1)


plot1(1) = gsn_csm_contour_map_polar(wks, trend10, res2)


plot2 = gsn_csm_contour_map_polar(wks, trend_delta5, res3)

resp1 = True
resp1@gsLineThicknessF = 10
y = new(720,"float")
y = 50.1
lon = f1->lon
pl1 = gsn_add_polyline(wks,plot1(0),lon,y,resp1)
pl2 = gsn_add_polyline(wks,plot1(1),lon,y,resp1)
pl3 = gsn_add_polyline(wks,plot2,lon,y,resp1)

resp = True
resp@gsnFrame = False 
resp@txString = ""
resp@gsnPanelLabelBar = True
resp@gsnPanelFigureStrings = (/"a","b"/)
resp@gsnPanelXWhiteSpacePercent = 1 ;左右两边的空白区域占比
resp@gsnPanelYWhiteSpacePercent = 0.1 ;上下两边的空白区域占比
resp@gsnPanelFigureStringsFontHeightF = 0.012
resp@gsnPanelFigureStringsPerimOn = False ;去掉方框的边框
resp@gsnPanelFigureStringsBackgroundFillColor = -1 ;方框背景
resp@gsnPanelRight    = 0.66
; resp@pmLabelBarOrthogonalPosF = -0.092
; resp@pmLabelBarParallelPosF = 0.55
resp@pmLabelBarWidthF = 0.6;colorbar width
resp@pmLabelBarHeightF = 0.06;colorbar height
resp@lbLabelFontHeightF = 0.01
resp@gsnMaximize          = True
resp@lbBoxLineThicknessF = 5

resp2 = True
resp2@gsnFrame = False
resp2@gsnPanelLabelBar = False
resp2@gsnPanelFigureStrings = (/"c"/)
resp2@gsnPanelFigureStringsPerimOn = False ;去掉方框的边框
resp2@gsnPanelFigureStringsBackgroundFillColor = -1 ;方框背景
resp2@gsnPanelLeft    = 0.67
; resp2@gsnPanelDebug      = True
resp2@gsnMaximize          = True

resp2@vpYF = 0.1

gsn_panel(wks, plot1, (/1,2/), resp)
gsn_panel(wks, plot2, (/1,1/), resp2)

frame(wks)

end
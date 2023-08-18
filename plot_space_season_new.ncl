begin

f1 = addfile("/mnt/hgfs/iap_down0711/t2m_CRU.nc","r")
f2 = addfile("/mnt/hgfs/iap_down0711/t2m_new.nc","r")
; print(f1)
; print(f2)
t2m_C = f1->t2m(420:779,{50:90},:) ;420:779 后三十年
t2m_M1 = f2->t2m(420:779,{50:90},:)

t2m_M1(:,{85:90},:) = -9999.0
t2m_M = t2m_M1
t2m_M(:,:,0:359) = t2m_M1(:,:,360:719)
t2m_M(:,:,360:719) = t2m_M1(:,:,0:359)
copy_VarMeta(t2m_M1,t2m_M)

t2m_C_DJF = month_to_season(t2m_C, "DJF")
t2m_C_JJA = month_to_season(t2m_C, "JJA")
t2m_M_DJF = month_to_season(t2m_M, "DJF")
t2m_M_JJA = month_to_season(t2m_M, "JJA")
t2m_C_MAM = month_to_season(t2m_C, "MAM")
t2m_C_SON = month_to_season(t2m_C, "SON")
t2m_M_MAM = month_to_season(t2m_M, "MAM")
t2m_M_SON = month_to_season(t2m_M, "SON")

t2m_C_DJF_ave = dim_avg_n_Wrap(t2m_C_DJF,0)
t2m_C_JJA_ave = dim_avg_n_Wrap(t2m_C_JJA,0)
t2m_M_DJF_ave = dim_avg_n_Wrap(t2m_M_DJF,0)
t2m_M_JJA_ave = dim_avg_n_Wrap(t2m_M_JJA,0)
t2m_C_MAM_ave = dim_avg_n_Wrap(t2m_C_MAM,0)
t2m_C_SON_ave = dim_avg_n_Wrap(t2m_C_SON,0)
t2m_M_MAM_ave = dim_avg_n_Wrap(t2m_M_MAM,0)
t2m_M_SON_ave = dim_avg_n_Wrap(t2m_M_SON,0)


t2m_delta_DJF_ave = t2m_M_DJF_ave - t2m_C_DJF_ave
t2m_delta_MAM_ave = t2m_M_MAM_ave - t2m_C_MAM_ave
t2m_delta_JJA_ave = t2m_M_JJA_ave - t2m_C_JJA_ave
t2m_delta_SON_ave = t2m_M_SON_ave - t2m_C_SON_ave

copy_VarMeta(t2m_M_DJF_ave,t2m_delta_DJF_ave)
copy_VarMeta(t2m_M_MAM_ave,t2m_delta_MAM_ave)
copy_VarMeta(t2m_M_JJA_ave,t2m_delta_JJA_ave)
copy_VarMeta(t2m_M_SON_ave,t2m_delta_SON_ave)

wks_type = "png"
wks_type@wkWidth = 5000
wks_type@wkHeight = 5000

wks = gsn_open_wks(wks_type,"t2m_space_season_MAM")

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
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -32
res@cnMaxLevelValF = 32
res@cnLevelSpacingF = 4
res@lbBoxLineThicknessF = 5
res@cnLinesOn = False

res1 = res
res1@gsnLeftString = "OBS(MAM)"
res1@gsnRightString = "~S~o~N~C"
res1@gsnCenterString = "~S~"
res1@gsnCenterStringOrthogonalPosF = 0.122
; res1@vpXF = 0.01
; res1@vpYF = 0.4

res2 = res
res2@gsnLeftString = "CAS-ESM(MAM)"
res2@gsnRightString = "~S~o~N~C"
res2@gsnCenterString = "~S~"
res2@gsnCenterStringOrthogonalPosF = 0.122
; res2@vpXF = 0.34
; res2@vpYF = 0.4

res3 = res
res3@gsnLeftString = "Bias(b-a)"
res3@gsnRightString = "~S~o~N~C"
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

res3@cnLevelSelectionMode = "ManualLevels"
res3@cnMinLevelValF = -12
res3@cnMaxLevelValF = 12
res3@cnLevelSpacingF = 2

plot1(0) = gsn_csm_contour_map_polar(wks, t2m_C_MAM_ave, res1)


plot1(1) = gsn_csm_contour_map_polar(wks, t2m_M_MAM_ave, res2)


plot2 = gsn_csm_contour_map_polar(wks, t2m_delta_MAM_ave, res3)

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
begin

f1 = addfile("/mnt/hgfs/iap_down0711/t2m_CRU.nc","r")
f2 = addfile("/mnt/hgfs/iap_down0711/t2m_new.nc","r")
; print(f1)
; print(f2)

t2m_C = f1->t2m
t2m_M = f2->t2m
; printVarSummary(t2m_C(419:780,:,:))
t2m_C_ave = dim_avg_n_Wrap(t2m_C(420:779,:,:), 0)
t2m_M_ave1 = dim_avg_n_Wrap(t2m_M(420:779,:,:), 0)

t2m_M_ave = t2m_M_ave1
t2m_M_ave(:,0:359) = t2m_M_ave1(:,360:719)
t2m_M_ave(:,360:719) = t2m_M_ave1(:,0:359)
copy_VarCoords(t2m_M_ave1, t2m_M_ave)

t2m_M_ave({85:90},:) = -9999.0
t2m_M_ave!0 = "lat"
t2m_M_ave!1 = "lon"

; lon = f2->LON
; lat = f2->LAT
; t2m_C_ave&lat = lat
; t2m_C_ave&lon = lon

printVarSummary(t2m_C_ave)
printVarSummary(t2m_M_ave)

t2m_delta_ave = t2m_M_ave - t2m_C_ave
copy_VarMeta(t2m_M_ave,t2m_delta_ave)
printVarSummary(t2m_delta_ave)
; lon = f2->lon
; lat = f2->lat
; t2m_delta_ave!0 = "lat"
; t2m_delta_ave!1 = "lon"
; t2m_delta_ave&lat = lat
; t2m_delta_ave&lon = lon

wks_type = "png"
wks_type@wkWidth = 5000
wks_type@wkHeight = 5000

wks = gsn_open_wks(wks_type,"test")

plot1 =  new(2,graphic)

res = True
res@cnLinesOn = False
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

res@lbBoxLineThicknessF = 5
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -28
res@cnMaxLevelValF = 28
res@cnLevelSpacingF = 4

res1 = res
res1@gsnLeftString = "OBS"
res1@gsnRightString = "~S~o~N~C"
res1@gsnCenterString = "~S~"
res1@gsnCenterStringOrthogonalPosF = 0.122
; res1@vpXF = 0.01
; res1@vpYF = 0.4

res2 = res
res2@gsnLeftString = "CAS-ESM"
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

plot1(0) = gsn_csm_contour_map_polar(wks, t2m_C_ave, res1)


plot1(1) = gsn_csm_contour_map_polar(wks, t2m_M_ave, res2)


plot2 = gsn_csm_contour_map_polar(wks, t2m_delta_ave, res3)

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
; maximize_output(wks, resp)
;请注意，NCL 始终将其图形绘制到方形窗口。如果您设置 这些资源变成非正方形，生成的图像将被调整大小 以适合给定矩形的最大正方形。
;所以除非改成pdf，用wkPaperHeightF还可以自定义下大小，不然难搞

; 1.  我使用ncl Panel 做了2行3列的展板，（1，1）和（2，1）y轴显示，其余图片y轴标签不显示。
; 由于坐标轴的有无导致每个图大小不一样，因此要控制每个子图属性res@gsnMaximize= False  ， 展板属性resPanel@gsnMaximize = True   。

; 2. 展板中各子图空白间距过大。如何调节？
; 注意的是，        resPanel@gsnPanelYWhiteSpacePercent         = 0.5  
;                         resPanel@gsnPanelXWhiteSpacePercent         = 0.2  
; 这两个是用来增加空白间隙距离的，官网上的例子有写。

; 可以通过调整各子图的位置减小间隙：
; resPanel@gsnPanelXF                                         = (/0.1,0.38,0.66,0.1,0.38,0.66/)     ；各个子图的横向位置
; resPanel@gsnPanelYF                                         = (/-1,-1,-1,-1,-1,-1/)    ；-1是默认值


; 3. NCL子图写上标号：
; resPanel@gsnPanelFigureStrings                                      = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)
; resPanel@gsnPanelFigureStringsFontHeightF                      = 0.02
; resPanel@amJust                                                        = "TopLeft"    ；位置，默认是在右下方
; resPanel@gsnPanelFigureStringsPerimOn                    = False          ;去掉方框的边框
; resPanel@gsnPanelFigureStringsBackgroundFillColor    =-1                ;方框背景设置为透明

end
begin

f1 = addfile("/mnt/hgfs/iap_down0711/t2m_CRU.nc","r")
f2 = addfile("/mnt/hgfs/iap_down0711/t2m_new.nc","r")
; print(f1)
; print(f2)

t2m_C = f1->t2m(180:779,{50:90},:) ;420:779 后三十年
t2m_M1 = f2->t2m(180:779,{50:90},:)

t2m_M1(:,{85:90},:) = -9999.0
t2m_M = t2m_M1
t2m_M(:,:,0:359) = t2m_M1(:,:,360:719)
t2m_M(:,:,360:719) = t2m_M1(:,:,0:359)

lat = f1->lat({50:90})
t2m_C_ave = wgt_areaave_Wrap(t2m_C, lat, 1, 1)
t2m_M_ave = wgt_areaave_Wrap(t2m_M, lat, 1, 1)

t2m_C_annual = month_to_annual(t2m_C_ave, 1);0是求和
t2m_M_annual = month_to_annual(t2m_M_ave, 1);0是求和

t2m_C_month = wgt_areaave_Wrap(clmMonTLL(t2m_C),lat,1,1)
t2m_M_month = wgt_areaave_Wrap(clmMonTLL(t2m_M),lat,1,1)

T_annual = new((/2,50/), "float")
T_annual(0,:) = t2m_C_annual
T_annual(1,:) = t2m_M_annual

T_month = new((/2,12/), "float")
T_month(0,:) = t2m_C_month
T_month(1,:) = t2m_M_month

wks_type = "png"
wks_type@wkWidth = 5000
wks_type@wkHeight = 5000
;;;;;;;;;;;;;;;;;;;;annual;;;;;;;;;;;;;;;;;;;;;;

x = ispan(0, 49, 1)

wks = gsn_open_wks(wks_type,"t2m_annual")

res = True
res@gsnDraw = False
res@gsnFrame = False

res@vpWidthF=0.85
res@vpHeightF=0.3
res@vpXF = 0.1
res@vpYF = 0.65

res@xyLineThicknessF = 15
res@tmBorderThicknessF = 20
; res@xyDashPatterns = (/0,1/);line patterns

res@trXMaxF = 49

res@tmXBMode = "explicit"
res@tmXBLabelsOn = True
res@tmXBValues = (/5,15,25,35,45/)
res@tmXBLabels = (/1970,1980,1990,2000,2010/)
res@tmXBMajorLengthF = 0.015 ;major 长度
res@tmXBLabelDeltaF = -0.8
res@tmXBMajorThicknessF = 5

res@tmYLMode = "explicit"
res@tmYLLabelsOn = True
res@tmYLValues = (/-10,-9,-8,-7,-6,-5,-4,-3/)
res@tmYLLabels = (/-10,-9,-8,-7,-6,-5,-4,-3/)
res@tmYLMajorLengthF = 0.015 ;major 长度
res@tmYLLabelDeltaF = -0.8
res@tmYLMajorThicknessF = 5

res@tmXTOn = False
res@tmYROn = False

res@gsnLeftString = "Annual Mean 2m Temperature"
res@gsnRightString = "~S~o~N~C"
res@gsnLeftStringFontHeightF = 0.015
res@gsnRightStringFontHeightF = 0.015

;legend
res@pmLegendDisplayMode = True
;size
res@pmLegendWidthF = 0.16
res@pmLegendHeightF = 0.08
res@pmLegendSide = "right"
res@lgPerimOn = False
;position
res@pmLegendOrthogonalPosF = -0.98
res@pmLegendParallelPosF = 0.87
res@xyExplicitLegendLabels = (/"OBS","CAS-ESM"/)

plot = gsn_csm_xy(wks,x,T_annual,res)

; y = conform_dims((/2,50/),dim_avg_n(T_annual,1),0)
; print(y(0,0)-y(1,0))
rc_C = regline(x,T_annual(0,:))
y_C = rc_C*x + rc_C@yintercept
rc_M = regline(x,T_annual(1,:))
y_M = rc_M*x + rc_M@yintercept
y1 = array_append_record(y_C,y_M,0)
y = reshape(y1, (/2,50/))
printVarSummary(y)

res1 = True
res1@gsnDraw = False
res1@gsnFrame = False

res1@xyLineThicknessF = 5
; res1@xyDashPattern =
plot1 = gsn_csm_xy(wks,x,y,res1)

overlay(plot, plot1)
draw(plot)
frame(wks)

;;;;;;;;;;;;;;;;;;;;month;;;;;;;;;;;;;;;;;;;;;;
; x = ispan(0, 11, 1)

; wks = gsn_open_wks(wks_type,"t2m_month")

; res = True
; ; res@gsnFrame = False
; ; res@gsnDraw = False
; res@vpWidthF=0.85
; res@vpHeightF=0.3
; res@vpXF = 0.1
; res@vpYF = 0.65

; res@xyLineThicknessF = 15
; res@tmBorderThicknessF = 20
; ; res@xyDashPatterns = (/0,1/);line patterns

; res@trXMaxF = 11
; res@tiYAxisString = ""

; res@tmXBMode = "explicit"
; res@tmXBLabelsOn = True
; res@tmXBValues = (/0,1,2,3,4,5,6,7,8,9,10,11/)
; res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; res@tmXBMajorLengthF = 0.015 ;major 长度
; res@tmXBLabelDeltaF = -0.8
; res@tmXBMajorThicknessF = 5

; ; res@tmYLMode = "explicit"
; ; res@tmYLLabelsOn = True
; ; res@tmYLValues = (/-8,-7,-6,-5,-4,-3,-2,-1/)
; ; res@tmYLLabels = (/-8,-7,-6,-5,-4,-3,-2,-1/)
; res@tmYLMinorOn = False
; res@tmYLMajorLengthF = 0.015 ;major 长度
; res@tmYLLabelDeltaF = -0.8
; res@tmYLMajorThicknessF = 5

; res@tmXTOn = False
; res@tmYROn = False

; res@gsnLeftString = "Monthly Average 2m Temperature"
; res@gsnRightString = "~S~o~N~C"
; res@gsnLeftStringFontHeightF = 0.015
; res@gsnRightStringFontHeightF = 0.015

; ;legend
; res@pmLegendDisplayMode = True
; ;size
; res@pmLegendWidthF = 0.16
; res@pmLegendHeightF = 0.08
; res@pmLegendSide = "right"
; res@lgPerimOn = False
; ;position
; res@pmLegendOrthogonalPosF = -1
; res@pmLegendParallelPosF = 0.87
; res@xyExplicitLegendLabels = (/"OBS","CAS-ESM"/)

; plot = gsn_csm_xy(wks,x,T_month,res)

; ; y = conform_dims((/2,12/),dim_avg_n(T_month,1),0)
; ; print(y(0,0)-y(1,0))

; ; res1 = True
; ; res1@gsnDraw = False
; ; res1@gsnFrame = False

; ; res1@xyLineThicknessF = 5
; ; ; res1@xyDashPattern = 2
; ; plot1 = gsn_csm_xy(wks,x,y,res1)

; ; overlay(plot, plot1)
; ; draw(plot)
; ; frame(wks)






end
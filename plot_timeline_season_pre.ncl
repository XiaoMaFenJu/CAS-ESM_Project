begin

f1 = addfile("/mnt/hgfs/iap_down0711/pre_CRU.nc","r")
f2 = addfile("/mnt/hgfs/iap_down0711/pre_new.nc","r")
; print(f1)
; print(f2)

t2m_C = f1->pre(180:779,{50:90},:) ;420:779 后三十年
t2m_M1 = f2->pre(180:779,{50:90},:)

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

t2m_C_DJF_ave = wgt_areaave_Wrap(t2m_C_DJF, lat, 1, 1) 
t2m_M_DJF_ave = wgt_areaave_Wrap(t2m_M_DJF, lat, 1, 1) 
t2m_C_JJA_ave = wgt_areaave_Wrap(t2m_C_JJA, lat, 1, 1) 
t2m_M_JJA_ave = wgt_areaave_Wrap(t2m_M_JJA, lat, 1, 1) 

T_DJF = new((/2,50/), "float")
T_DJF(0,:) = t2m_C_DJF_ave
T_DJF(1,:) = t2m_M_DJF_ave

T_JJA = new((/2,50/), "float")
T_JJA(0,:) = t2m_C_JJA_ave
T_JJA(1,:) = t2m_M_JJA_ave

t2m_C_MAM = month_to_season(t2m_C, "MAM")
t2m_C_SON = month_to_season(t2m_C, "SON")
t2m_M_MAM = month_to_season(t2m_M, "MAM")
t2m_M_SON = month_to_season(t2m_M, "SON")

t2m_C_MAM_ave = wgt_areaave_Wrap(t2m_C_MAM, lat, 1, 1) 
t2m_M_MAM_ave = wgt_areaave_Wrap(t2m_M_MAM, lat, 1, 1) 
t2m_C_SON_ave = wgt_areaave_Wrap(t2m_C_SON, lat, 1, 1) 
t2m_M_SON_ave = wgt_areaave_Wrap(t2m_M_SON, lat, 1, 1) 

T_MAM = new((/2,50/), "float")
T_MAM(0,:) = t2m_C_MAM_ave
T_MAM(1,:) = t2m_M_MAM_ave

T_SON = new((/2,50/), "float")
T_SON(0,:) = t2m_C_SON_ave
T_SON(1,:) = t2m_M_SON_ave

wks_type = "png"
wks_type@wkWidth = 5000
wks_type@wkHeight = 5000

;;;;;;;;;;;;;;;;;;;;DJF;;;;;;;;;;;;;;;;;;;;;;

x = ispan(0, 49, 1)

wks = gsn_open_wks(wks_type,"pre_timeline_SON")

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
res@tiYAxisString = ""

res@trXMaxF = 49

res@tmXBMode = "explicit"
res@tmXBLabelsOn = True
res@tmXBValues = (/5,15,25,35,45/)
res@tmXBLabels = (/1970,1980,1990,2000,2010/)
res@tmXBMajorLengthF = 0.015 ;major 长度
res@tmXBLabelDeltaF = -0.8
res@tmXBMajorThicknessF = 5

; res@tmYLMode = "explicit"
; res@tmYLLabelsOn = True
; res@tmYLValues = (/-8,-7,-6,-5,-4,-3,-2,-1/)
; res@tmYLLabels = (/-8,-7,-6,-5,-4,-3,-2,-1/)
res@tmYLMajorLengthF = 0.015 ;major 长度
res@tmYLLabelDeltaF = -0.8
res@tmYLMajorThicknessF = 5
res@tmYLMinorOn = False

res@tmXTOn = False
res@tmYROn = False

res@gsnRightString = "mm/month"
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

res@gsnLeftString = "Average Precipitation (SON)"
plot = gsn_csm_xy(wks,x,T_SON,res)


; y = conform_dims((/2,50/),dim_avg_n(T_DJF,1),0)
; print(y(0,0)-y(1,0))
rc_C = regline(x,T_SON(0,:))
y_C = rc_C*x + rc_C@yintercept
rc_M = regline(x,T_SON(1,:))
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
end
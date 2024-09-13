load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin

syear    = 1870
eyear    = 2022
ymdStrt1 = syear*100+1
ymdLast1 = eyear*100+12
ymdbetween=2000*100+1
mm=(eyear-syear+1)*12
mm1=(2018-2000)*12;(2000-1980)*12;

opendir="/data0/fjiang/data/ocean/"
resultdir="/data0/fjiang/scripts/contrast/"


a= addfile(opendir+"HadISST_sst.nc","r");sst.mnmean.nc
time1      = a->time
ymd1      = cd_calendar(time1, -1)
iStrt3    = ind(ymd1.eq.ymdStrt1)
iLast3    = ind(ymd1.eq.ymdLast1)
ibe=ind(ymd1.eq.ymdbetween)
delete(time1)
delete(ymd1)
ha   = a->sst(iStrt3:iLast3,::-1,:)
ha:=lonPivot(ha,0.5)

hor= where(abs(ha).ge.50,ha@_FillValue,ha)
hor= where(ismissing(hor),0,hor)
hor=rmMonAnnCycTLL(hor)
copy_VarMeta(ha, hor)
hor:=dtrend_n(hor,False,0)
copy_VarMeta(ha, hor)
h=hor

sstannual=month_to_annual(hor,1)
printVarSummary(sstannual)


syear1=1980
syear2=1958
year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)


rc1    = regCoef_n(year1,sstannual(syear1-syear:eyear-syear,:,:),0,0)
rc2    = regCoef_n(year2,sstannual(syear2-syear:eyear-syear,:,:),0,0)

copy_VarMeta(sstannual(0,:,:),rc1)
copy_VarMeta(sstannual(0,:,:),rc2)
dims=dimsizes(sstannual)
confi1=student_t(rc1@tval,rc1@nptxy-2)
confi1:=onedtond(confi1,(/dims(1),dims(2)/))
rc1@long_name = "Trend"
rc1:= rc1*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(sstannual(0,:,:),rc1)
copy_VarMeta(rc1,confi1)

confi2=student_t(rc2@tval,rc2@nptxy-2)
confi2:=onedtond(confi2,(/dims(1),dims(2)/))
rc2@long_name = "Trend"
rc2:= rc2*10;(eyear-syear2+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(sstannual(0,:,:),rc2)
copy_VarMeta(rc2,confi2)



rc1_all=new((/111,180,360/),float)
rc2_all=new((/90,180,360/),float)
confi1_all=new((/111,180,360/),float)
confi2_all=new((/90,180,360/),float)

patterncor1=new((/153/),float)
patterncor1_a=new((/153/),float)
patterncor1_b=new((/153/),float)

patterncor2=new((/153/),float)
patterncor2_a=new((/153/),float)
patterncor2_b=new((/153/),float)


lats=-60
late=60
lons=120
lone=290
xx1=year1
printVarSummary(xx1)
do i=0,110
  rc1_temp=regCoef_n(year1,sstannual(0+i:42+i,:,:),0,0)
  copy_VarMeta(sstannual(0,:,:),rc1_temp)
  rc1_all(i,:,:)=rc1_temp
  confi1_temp=student_t(rc1_temp@tval,rc1_temp@nptxy-2)
  confi1_all(i,:,:)=onedtond(confi1_temp,(/dims(1),dims(2)/))
  copy_VarMeta(rc1,confi1_all(i,:,:))

  patterncor1(i+42)=pattern_cor(rc1({lats:late},{lons:lone}),rc1_all(i,{lats:late},{lons:lone}), 1.0, 0)
  patterncor1_a(i+42)=pattern_cor(rc1({lats+30:late-30},{lons:lone}),rc1_all(i,{lats+30:late-30},{lons:lone}), 1.0, 0)
  patterncor1_b(i+42)=pattern_cor(rc1({lats+50:late-50},{lons:lone}),rc1_all(i,{lats+50:late-50},{lons:lone}), 1.0, 0)

end do

xx2=year2
do i=0,88
  rc2_temp=regCoef_n(year2,sstannual(0+i:64+i,:,:),0,0)
  copy_VarMeta(sstannual(0,:,:),rc2_temp)
  rc2_all(i,:,:)=rc2_temp
  confi2_temp=student_t(rc2_temp@tval,rc2_temp@nptxy-2)
  confi2_all(i,:,:)=onedtond(confi2_temp,(/dims(1),dims(2)/))
  copy_VarMeta(rc2,confi2_all(i,:,:))
  patterncor2(i+64)=pattern_cor(rc2({lats:late},{lons:lone}),rc2_all(i,{lats:late},{lons:lone}), 1.0, 0)
  patterncor2_a(i+64)=pattern_cor(rc2({lats+30:late-30},{lons:lone}),rc2_all(i,{lats+30:late-30},{lons:lone}), 1.0, 0)
  patterncor2_b(i+64)=pattern_cor(rc2({lats+50:late-50},{lons:lone}),rc2_all(i,{lats+50:late-50},{lons:lone}), 1.0, 0)

end do

print(patterncor1)
print(patterncor2)


rc1_all=rc1_all*10
copy_VarMeta(rc2_all,rc1_all)
rc2_all=rc2_all*10
copy_VarMeta(rc1_all,rc2_all)

;hor=dtrend_msg_n(hor&time,hor,False,False,0)

copy_VarMeta(ha, hor)
h=hor

;sstannual=month_to_annual(h,1)
;printVarSummary(sstannual)
syear1=1870
eyear1=2022
year1:= ispan(syear1,eyear1,1)


ipo_obs  =wgt_areaave_Wrap(h(:,{-10:10},{170:270}), 1.0,1.0,0)-0.5*(wgt_areaave_Wrap(h(:,{25:45},{140:215}), 1.0,1.0,0)+wgt_areaave_Wrap(h(:,{-50:-15},{150:200}), 1.0,1.0,0))
fcbt =1/(12.*13)
cft1 =ezfftf_n(ipo_obs, 0)
fcft =fspan(0, 0.5,mm/2)
ifcft1=ind(fcft.gt.fcbt)  ;;decadaal
cft1(:,ifcft1)=0.0
ipo_obs:=ezfftb_n(cft1,cft1@xbar,1);
printVarSummary(ipo_obs)



    ca    = 99999.        ; band start (longer period)
    cb    = 12*13.0        ; band end

    fca   = 1.0/ca      ; 'left'  frequency
    fcb   = 1.0/cb      ; 'right' frequency

    dimss  = 0           ; 'time' dimension

    opt   = True        ; options to set
    opt@m  = 6 ; time series of filtered and envelope values
    opt@return_envelope  = False ; time series of filtered and envelope values

;ipo_obs:= bw_bandpass_filter (ipo_obs,fca,fcb,opt,dimss)
printVarSummary(ipo_obs)

rc1_ipo= regCoef_n(ipo_obs,h,0,0)
copy_VarMeta(h(0,:,:),rc1_ipo)

rc1_avg=wgt_areaave_Wrap(rc1_ipo({-5:5},{140:170}), 1.0,1.0,0)-wgt_areaave_Wrap(rc1_ipo({-5:5},{190:270}), 1.0,1.0,0)



ew_avg =wgt_areaave_Wrap(h(:,{-5:5},{140:170}), 1.0,1.0,0)-wgt_areaave_Wrap(h(:,{-5:5},{190:270}), 1.0,1.0,0)



dims=dimsizes(h)
confi1_ipo=student_t(rc1_ipo@tval,rc1_ipo@nptxy-2)
confi1_ipo:=onedtond(confi1_ipo,(/dims(1),dims(2)/))
rc1_ipo@long_name = "Trend"
;rc1=-rc1
copy_VarMeta(h(0,:,:),rc1_ipo)
copy_VarMeta(rc1_ipo,confi1_ipo)


fname1 = resultdir + "figuredata/IPO.nc";;eq.
system("/bin/rm -f " + fname1)
fout1 = addfile(fname1, "c")
fout1->rc1_ipo = rc1_ipo
fout1->ipo_obs =ipo_obs


ipo_ew=rc1_avg*month_to_annual(ipo_obs,1)
ew=month_to_annual(ew_avg,1)

rc1_all_ipo=new((/153/),float)
rc2_all_ipo=new((/153/),float)
rc1_raw_ipo=new((/153/),float)
rc2_raw_ipo=new((/153/),float)
syear1=1980
syear2=1958
year1:=ispan(syear1,eyear1,1)
year2:=ispan(syear2,eyear1,1)

do i=0,110
  rc1_all_ipo(i+42)=10*regCoef_n(year1,ipo_ew(0+i:42+i),0,0)
  rc1_raw_ipo(i+42)=10*regCoef_n(year1,ew(0+i:42+i),0,0)

end do

do i=0,88
  rc2_all_ipo(i+64)=10*regCoef_n(year2,ipo_ew(0+i:64+i),0,0)
  rc2_raw_ipo(i+64)=10*regCoef_n(year2,ew(0+i:64+i),0,0)

end do


d  =0
sm =3
pct=0.10
sdof =specx_anal(ipo_obs,d, sm, pct)
splt =specx_ci(sdof,0.05,0.95)


ipo_obs!0="time"
time=ispan(0,mm-1,1)
ipo_obs&time=time

year=time
N  = mm
  ;; 小波计算
  mother  = 0       ; 母小波类型，通常为0，即'Morlet'小波。其余两中被分别为1，'Paul'小波和2，'DOG' (derivative of Gaussian)小波
  dt      = 1.       ; 数组中数值之间的时间间隔，通常为1。本例中表示间隔1年。
  param   = 6     ; 母小波参数。 如果param < 0，则使用默认数值，即采用'Morlet'小波时为6；Paul'小波为4;'DOG'小波为2
  s0      = dt      ; 'Morlet'小波s0 = dt ; 'Paul'小波s0 = dt/4
  dj      = 0.25    ; 常用设定
  jtot    = 1+floattointeger(((log10(N*dt/s0))/dj)/log10(2.)) ; 常用设定
  npad    = N       ; 常用设定
  nadof   = 0       ; 常用设定
  noise   = 1       ; 常用设定,h红噪声检验
  siglvl  = .05  ; 置信度水平
  isigtest= 0       ; 采用chi-square 检验；若为1则是对全部波谱进行时间平均检验

w1=wavelet(ipo_obs,mother,dt,param,s0,dj,jtot,npad,noise,isigtest,siglvl,nadof)
power            = onedtond(w1@power,(/jtot,N/))  ; ¹¦ÂÊÆ×
power!0          = "period"                      ; Y axis
power&period     = w1@period

power!1          = "year"                        ; X axis
power&year      = year

power@long_name  = "Power Spectrum"
power@units      = "1/unit-freq"

;计算显著性 ( >= 1 则显著)
SIG              = power                         ; 复制元数据
SIG              = power/conform (power,w1@signif,0)
SIG@long_name    = "Significance"
SIG@units        = " "

  wks = gsn_open_wks("pdf",resultdir+"v6-s1");1979-2018")
   gsn_define_colormap(wks,"MPL_RdBu")
   ;gsn_define_colormap(wks,"MPL_RdBu");"CBR_drywet");
  gsn_reverse_colormap(wks);

   res = True
   res@gsnDraw   = False
   res@gsnFrame  = False
   res@gsnAddCyclic      = True
   res@mpMinLonF  = 0 ;lon2d(0,0)
   res@mpCenterLonF=180
   res@mpMaxLonF = 360 ;lon2d(0,ilon-1)
   res@mpFillOn=True
   res@mpLandFillColor="gray70"
   res@mpOceanFillColor=-1
   res@mpFillDrawOrder="PostDraw"
   res@tmXBMinorOn=False
   ;res@tmYLLabels =(/"30S","20S","10S","0","10N","20N","30N"/)
   res@tmYLMinorOn=False
   res@cnInfoLabelOn=False
   res@mpGeophysicalLineThicknessF=1
   res@mpGeophysicalLineColor="gray70"
   ;res@pmTickMarkDisplayMode = "Always"
   res@tmXBLabelFontHeightF  = 0.007
   res@tmYLLabelFontHeightF  = 0.007

   res@cnFillOn    = True
   res@cnLinesOn = False

   res@tmXBMajorLengthF=0.004
   res@tmYLMajorLengthF=0.004
   res@tmXBMajorOutwardLengthF=0.004
   res@tmYLMajorOutwardLengthF=0.004
   res@mpProjection="Robinson"
   res@mpPerimOn=False
   res@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
   ;res@cnFillColors=(/111,100,96,87,77,67,64,54,44,34,30,20/);
   ;res@cnFillColors=(/10,15,20,25,30,35,40,45,50,57,0,0,73,80,85,90,95,100,105,110,115,120/)
   res@cnFillColors=(/12,20,35,50,57,0,0,73,80,95,105,115/)

   res@gsnLeftStringFontHeightF=0.01
   res@gsnLeftStringOrthogonalPosF=0.02;0.005
   res@gsnRightStringFontHeightF=0.012
   res@gsnRightStringOrthogonalPosF=-0.003;5;0.005
   res@cnLineLabelsOn=False
   res@gsnRightString="";(~S~o~N~C/decade)"
   res@gsnLeftString=""
   res@tiMainString="(a) IPO Spatial Pattern ";+"~F34~4 "+"~F22~-1 "
   res@tiMainOffsetYF=-0.01
   res@tiMainFontHeightF=0.012
   plot1=new(8,graphic)
   ;dum1=new(4,graphic)
   res@vpXF=-0.14
   res@vpYF=0.97
   res@vpHeightF= 0.2            ; change aspect ratio of plot
   res@cnLevels             = 0.1*ispan(-10,10,2 )

   res@lbLabelBarOn=False
   res@lbLabelBarOn=True
   res@lbOrientation="Vertical"
   ;res@lbBoxEndCapStyle="TriangleBothEnds"
   res@pmLabelBarHeightF=0.16
   res@pmLabelBarWidthF=0.03
   res@lbLabelStride=2.
   res@lbBoxLinesOn=True
   res@pmLabelBarOrthogonalPosF=0.05
   res@lbLabelFontHeightF=0.01
   rc1_p=-rc1_ipo
   copy_VarMeta(rc1,rc1_p)
   plot1(0) = gsn_csm_contour_map(wks,rc1_p,res)

   ;res@pmLabelBarOrthogonalPosF=0.2
   ;res@pmLabelBarParallelPosF=-0.11
   ;res@cnLevels             = 0.02*ispan(-10,10,1 )

   res@gsnLeftStringOrthogonalPosF=-0.02;0.005


  sres1=True
  sres1@gsnDraw=False
  sres1@gsnFrame=False
  sres1@cnFillOn=True
  sres1@cnLinesOn=False
  sres1@cnMonoFillPattern=False
  sres1@gsnLeftString=""
  sres1@gsnRightString=""
  sres1@cnLevelSelectionMode="ExplicitLevels"
  sres1@cnLevels=(/0,0.05/)
  sres1@cnFillPatterns=(/-1,17,-1/)
  sres1@cnFillColors=(/"black","black","black"/)
  sres1@cnLineLabelsOn=False
  sres1@cnInfoLabelOn=False
  sres1@lbLabelBarOn=False
  plot2=new(6,graphic)
;
  plot2(0)=gsn_csm_contour(wks,confi1_ipo,sres1)



  overlay(plot1(0),plot2(0))


  draw(plot1(0))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
res1                 = True
res1@gsnDraw         = False
res1@gsnFrame        = False
res1@vpHeightF = 0.18
res1@vpWidthF  = 0.39
res1@trXMinF   =  0.;1870.
res1@trXMaxF   =  mm-1;2022.
res1@trYMaxF   = 1.
res1@trYMinF   = -1.
res1@tiYAxisString = " "
res1@tiXAxisString = " "
res1@tiMainString  = " "
res1@tmXBMajorThicknessF     = 1.5
res1@tmYLMajorThicknessF     = 1.5
res1@tmXBMinorThicknessF     = 1.0
res1@tmYLMinorThicknessF     = 1.0
res1@tmYLMinorOn=False
res1@tmYLValues    = 0.5*ispan(-5,5,1)
res1@tmYLLabels    = 0.5*ispan(-5,5,1);
res1@tmXBValues    = ispan((1870-1870),(2020-1870)*12,30*12)
res1@tmXBLabels    = ispan(1870,2020,30);(/"1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015","2020"/)
res1@tmXBLabelsOn=True
res1@tmXBMode="ExplicitLevels"
res1@gsnYRefLine=0.
res1@tmXBMajorLengthF=0.006
res1@tmYLMajorLengthF=0.004
res1@tmXBMajorOutwardLengthF=0.006
res1@tmYLMajorOutwardLengthF=0.004
res1@tmXBMinorLengthF        = 0.004
res1@tmYLMinorLengthF        = 0.004
res1@tmXBMinorOn=True
;res@tmXBMinorValues=ispan(1958,2022,1)

res1@gsnLeftString=" "
res1@gsnLeftStringFontHeightF= 0.012
res1@gsnLeftStringOrthogonalPosF=0.00;5
res1@tmXBOn = True
;res@tmYLOn = True
;res@tmYLLabelsOn = True
;res@tmYLMinorOn = True
;res@tmYLLabelFontColor = "black"
res1@tiXAxisString         = "";End year of the trend"
res1@tiXAxisFontHeightF=0.01
res1@tiYAxisFontColor      = "black"
res1@tmYLOn    = True
res1@tmYLMode          = "Explicit"         ; explicit labels
res1@tmXBMode          = "Explicit"         ; explicit labels
res1@tmXBOn    = True

res1@tmYLLabelsOn=True

res1@tmXBLabelFontHeightF  = 0.008
res1@tmYLLabelFontHeightF  = 0.008
res1@tiYAxisSide           = "Left"
res1@tiYAxisFontHeightF    = 0.01
res1@xyLineColor      = "orangered3"
res1@tmXBBorderOn=True
res1@tiYAxisFontHeightF=0.012

res1@tmXTOn=False
res1@tmYROn=False
res1@gsnLeftString=""
res1@tiMainString="(b) IPO Timeseries"
res1@tiMainOffsetYF=-0.01
res1@tiMainFontHeightF=0.012
res1@gsnRightString=""
res1@xyLineThicknessF = 2.
res1@vpXF=0.065
res1@vpYF=0.725
res1@tiMainFont=22
res1@vpXF=0.55
res1@vpYF=0.97
res1@gsnXYBarChart = True            ; Create bar plot
res1@gsnXYBarWidth=0.2
res1@gsnBelowYRefLineBarColors="gray"
res1@gsnAboveYRefLineBarColors="gray"
res1@xyLineColor=-1
res1@tiXAxisFontHeightF=0.01
res1@tiXAxisString="Year"
res1@tiYAxisFontHeightF=0.01
res1@tiYAxisString="[~S~o~N~C]"


plot1(1) = gsn_csm_xy(wks,ispan(0,mm-1,1),ipo_obs,res1)


draw(plot1(1))
lgres=True
amres=True
lgres@lgLineColors=(/"deepskyblue3","deepskyblue3","orangered3","orangered3"/)
lgres@lgLeftMarginF=0.1
lgres@lgRightMarginF=0.1

lgres@lgDashIndexes=(/2.,0.,2.,0./)
lgres@lgLineThicknesses=(/1.5,1.,1.5,1./)
lgres@lgMonoLineThickness=False
lgres@vpWidthF=0.162
lgres@vpHeightF=0.07
lgres@lgLabelFontHeightF=0.12
lgres@lgLineDashSegLenF=0.04
lgres@lgPerimThicknessF=2.
lgres@lgPerimOn=False
lgres@lgPerimFill="SolidFill"
lgres@lgPerimFillColor="white"
;lgres@lgLabelOffsetF=-0.1
lgres@lbBoxMajorExtentF=0.3
lbid1=gsn_create_legend(wks,4,(/"   IPO Contribution to 65-year trend","   65-year trend","   IPO Contribution to 43-year trend","   43-year trend"/),lgres)

amres@amParallelPosF=-0.28
amres@amOrthogonalPosF=-0.28
;annoid2=gsn_add_annotation(plot1(3),lbid2,amres)

;draw(plot1(3))


res@lbLabelBarOn=True
aa=-rc1_all(42-42,:,:)
copy_VarMeta(rc1_all(0,:,:),aa)
res@vpXF=-0.14
res@vpYF=0.7
res@gsnLeftString=""
res@tiMainString="(c) HadISST SST Trend "+"~F34~4 "+"~F22~-1 P1 (1870-1912)  "

res@tiMainOffsetYF=-0.01
res@tiMainFontHeightF=0.012
res@cnLevels:= 0.02*ispan(-12,12,2 )
res@cnFillColors:=(/12,20,30,40,50,57,0,0,73,80,90,96,105,115/)

plot1(2) = gsn_csm_contour_map(wks,aa,res)

res@tiMainString="(d) HadISST SST Trend P2 (1914-1956)  "
;res@cnLevels             = 0.02*ispan(-10,10,1 )
res@vpXF=0.34
res@vpYF=0.7
res@lbLabelBarOn=True

plot1(3) = gsn_csm_contour_map(wks,rc1_all(86-42,:,:),res)

res@tiMainString="(e)  HadISST SST Trend "+"~F34~4 "+"~F22~-1 P3 (1942-1984)  "
bb=-rc1_all(114-42,:,:)
copy_VarMeta(rc1_all(0,:,:),bb)
res@vpXF=-0.14
res@vpYF=0.43

plot1(4) = gsn_csm_contour_map(wks,bb,res)
res@vpXF=0.34
res@vpYF=0.43

res@lbLabelBarOn=True
;res@lbBoxEndCapStyle="TriangleBothEnds"
;res@pmLabelBarParallelPosF=1.2
res@tiMainString="(f)  HadISST SST Trend "+"~F34~4 "+"~F22~-1 P4 (1870-1934) "
cc=-rc2_all(64-64,:,:)
copy_VarMeta(rc1_all(0,:,:),cc)
;res@tiMainOffsetYF=-0.03


plot1(5) = gsn_csm_contour_map(wks,cc,res)



plot2(2)=gsn_csm_contour(wks,confi1_all(42-42,:,:), sres1)
plot2(3)=gsn_csm_contour(wks,confi1_all(86-42,:,:), sres1)
plot2(4)=gsn_csm_contour(wks,confi1_all(114-42,:,:),sres1)
plot2(5)=gsn_csm_contour(wks,confi2_all(64-64,:,:), sres1)

overlay(plot1(2),plot2(2))
overlay(plot1(3),plot2(3))
overlay(plot1(4),plot2(4))
overlay(plot1(5),plot2(5))


draw(plot1(2))
draw(plot1(3))
draw(plot1(4))
draw(plot1(5))


txres=True
txres@txFontHeightF=0.0105
gsn_text_ndc(wks, "[~S~o~N~C]",0.456,0.77, txres)
gsn_text_ndc(wks, "[~S~o~N~C/decade]",0.46,0.5, txres)
gsn_text_ndc(wks, "[~S~o~N~C/decade]",0.46,0.234, txres)
gsn_text_ndc(wks, "[~S~o~N~C/decade]",0.945,0.5, txres)
gsn_text_ndc(wks, "[~S~o~N~C/decade]",0.945,0.234, txres)

frame(wks)

  pres=True
  pres@gsnPanelCenter=False
  pres@gsnPanelYWhiteSpacePercent=7
  pres@gsnPanelXWhiteSpacePercent=3
  pres@gsnPanelRight=0.95

  ;pres@gsnPanelLabelBar =True
  ;pres@gsnPanelFigureStrings =(/"a-spring","b-spring","c-summer","d-summer","e-autumn","f-autumn","g-winter","h-winter"/);,"c","d"
  ;pres@gsnPanelFigureStrings =(/"a","b","c"/)
  ;pres@amJust="TopLeft"
  ;gsn_panel(wks,plot1,(/4,2/),pres)
  ;delete(wks)

   end

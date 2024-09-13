load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;;;ipo_obs
syear    = 1958
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

ipo_obs  =wgt_areaave_Wrap(h(:,{-10:10},{170:270}), 1.0,1.0,0)-0.5*(wgt_areaave_Wrap(h(:,{25:45},{140:215}), 1.0,1.0,0)+wgt_areaave_Wrap(h(:,{-50:-15},{150:200}), 1.0,1.0,0))
fcbt =1/(12.*13)
cft1 =ezfftf_n(ipo_obs, 0)
fcft =fspan(0, 0.5,mm/2)
ifcft1=ind(fcft.gt.fcbt)  ;;decadaal
cft1(:,ifcft1)=0.0
ipo_obs:=-ezfftb_n(cft1,cft1@xbar,1);
printVarSummary(ipo_obs)
ipo_obs:=dim_standardize_n_Wrap(ipo_obs,1,0)
rc_ipo=120*regCoef_n(ispan(0,43*12-1,1),ipo_obs(22*12::),0,0)



opendir="/data0/fjiang/data/atmos/NCEP/"
fname1 = opendir + "air.mon.mean.nc";;eq.
a = addfile(fname1, "r")
syear1=1980
syear2=1958

eyear    = 2022
ymdStrt1 = syear1*100+1
ymdStrt2 = syear2*100+1

ymdLast1 = eyear*100+12
mm1=(eyear-syear1+1)*12
mm2=(eyear-syear2+1)*12
year1    = ispan(syear2,eyear,1)


time1      = a->time
ymd1      = cd_calendar(time1, -1)
iStrt1    = ind(ymd1.eq.ymdStrt1)
iStrt2    = ind(ymd1.eq.ymdStrt2)
iLast1    = ind(ymd1.eq.ymdLast1)


tt=a->air(iStrt2:iLast1,:,::-1,:)
printVarSummary(tt)



temperature=dim_avg_n_Wrap(tt(:,{1000:50},{-5:5},:),2)
temperature:=dim_avg_n_Wrap(temperature(:,:,{120:270})   ,2)


level=temperature&level

temperature_annual=month_to_annual(temperature,1)
temperature_Detrend=dtrend_msg_n (ispan(0,779,1),temperature,True,False,0)
copy_VarMeta(temperature,temperature_Detrend)
temperature=temperature_Detrend
temp_trend  =regCoef_n(year1,temperature_annual,0,0)
temp_ipo  =regCoef_n(ipo_obs,temperature,0,0)
copy_VarMeta(temperature(0,:),temp_ipo)
copy_VarMeta(temperature(0,:),temp_trend)
printVarSummary(temp_ipo)

dims=dimsizes(temperature)
confi1=student_t(temp_ipo@tval,temp_ipo@nptxy-2)
;confi1:=onedtond(confi1,(/dims(1),dims(2)/))
;copy_VarMeta(temperature(0,:,:),confi1)

confi2=student_t(temp_trend@tval,temp_trend@nptxy-2)
;confi2:=onedtond(confi2,(/dims(1),dims(2)/))
;copy_VarMeta(temperature(0,:,:),confi2)

temp_ipo   =temp_ipo*rc_ipo
temp_trend =temp_trend*10
copy_VarMeta(temperature(0,:),temp_ipo)
copy_VarMeta(temperature(0,:),temp_trend)


temp_ipo_sig=where(confi1.lt.0.05,temp_ipo,temp_ipo@_FillValue)
temp_trend_sig=where(confi2.lt.0.05,temp_trend,temp_ipo@_FillValue)

print(temp_trend_sig)
print(temp_ipo_sig)

  wks = gsn_open_wks("eps",resultdir+"v7-s4-mean-avg");1979-2018")
   gsn_define_colormap(wks,"MPL_RdBu")
   ;gsn_define_colormap(wks,"MPL_RdBu");"CBR_drywet");
  gsn_reverse_colormap(wks);

  res                   = True                       ; plot mods desired
  res@vpWidthF=0.6
  res@vpHeightF=0.49

  res@gsnDraw=False
  res@gsnFrame  = False
  res@tiMainFont=22
  res@tiMainString      = "Equatorial Pacific Temperature Trend"             ; add title
  res@tiMainOffsetYF=-0.01
  res@tiMainFontHeightF=0.018
  res@gsnLeftStringFontHeightF=0.014
  res@trYReverse        = True                       ; reverse Y-axis
  res@xyDashPatterns    = 0                         ; choose dash patterns

  res@tmXBMajorLengthF=0.004
  res@tmYLMajorLengthF=0.004
  res@tmXBMajorOutwardLengthF=0.004
  res@tmYLMajorOutwardLengthF=0.004
  res@tmXBOn=True
  res@tmYLOn=True
  res@tmXTOn=False
  res@tmYROn=False

  res@tmXBMinorOn=False
  ;res@tmYLLabels =(/"30S","20S","10S","0","10N","20N","30N"/)
  res@tmYLMinorOn=False
  res@tmYLLabelsOn=True
  res@tmXBLabelsOn=True
  res@tmXBMode ="Explicit"
  res@tmYLMode ="Explicit"
  ;res@trYMinF=50;
  ;res@trYMaxF=1000;
  res@trYLog=True
  res@trXMaxF=.6;
  res@trXMinF=-.6;
  res@gsnXRefLine=0.
  res@tmYLValues            = level;(/10,20,50,80,100,150/)
  res@tmYLLabels            = level;(/10,20,50,80,100,150/)
  res@tmXBValues=.1*ispan(-6,6,3)
  res@tmXBLabels=.1*ispan(-6,6,3)
  res@tmXBLabelFontHeightF     = 0.011
  res@tmYLLabelFontHeightF     = 0.011
  res@tmYLValues            = level
  res@tmYLabels             = level
  res@gsnLeftStringFontHeightF=0.0112
  res@gsnLeftStringOrthogonalPosF=0.012;0.005
  res@gsnRightStringFontHeightF=0.012
  res@gsnRightStringOrthogonalPosF=0.005;0.005
  res@xyLineColor="orangered3"
  res@xyDashPattern=16
  res@xyLineThicknessF=3.
  res@tiYAxisString    = "level"             ; turn off axis label
  res@tiXAxisString    = "[~S~o~N~C/decade]"             ; turn off axis label
  res@tiXAxisFontHeightF=0.015

  res@tiYAxisFontHeightF=0.015
  res@xyMarkLineMode="MarkLines"
  res@xyMarkerColor="orangered3"

  res@xyMarker=4
  plot=new(4,graphic)
  res@vpXF=0.15
  res@vpYF=0.65
  print(level)
   plot(0)   = gsn_csm_xy (wks,temp_ipo,level,res) ; create plot
   res@xyLineColor="deepskyblue3"
   res@xyMarkerColor="deepskyblue3"

   plot(1)   = gsn_csm_xy (wks,temp_trend,level,res) ; create plot
overlay(plot(0),plot(1))
res@xyMarker=16
res@xyMarkerColor="orangered3"
res@xyMarkLineMode="Markers"

plot(2)   = gsn_csm_xy (wks,temp_ipo_sig,level,res) ; create plot
res@xyMarkerColor="deepskyblue3"

plot(3)   = gsn_csm_xy (wks,temp_trend_sig,level,res) ; create plot
overlay(plot(0),plot(2))
overlay(plot(0),plot(3))



lgres=True
amres=True
lgres@lgLineColors=(/"deepskyblue3","orangered3"/)
lgres@lgLeftMarginF=0.1
lgres@lgRightMarginF=0.1

lgres@lgDashIndexes=(/0,0/)
lgres@lgLineThicknesses=(/4,4/)
lgres@lgMonoLineThickness=False
lgres@vpWidthF=0.17
lgres@vpHeightF=0.063
lgres@lgLabelFontHeightF=0.13
lgres@lgLineDashSegLenF=0.04
lgres@lgPerimThicknessF=2.
lgres@lgPerimOn=False
lgres@lgPerimFill="SolidFill"
lgres@lgPerimFillColor="white"
;lgres@lgLabelOffsetF=-0.1
lgres@lbBoxMajorExtentF=0.3
lbid1=gsn_create_legend(wks,2,(/"  Emerging Trend","  IPO-related Trend"/),lgres)
amres@amParallelPosF=0.3
amres@amOrthogonalPosF=-0.41
annoid1=gsn_add_annotation(plot(0),lbid1,amres)



draw(plot(0))


  pres=True
  pres@gsnPanelCenter=False
  pres@gsnPanelYWhiteSpacePercent=10
  ;pres@gsnPanelLabelBar =True
  ;pres@gsnPanelFigureStrings =(/"a-spring","b-spring","c-summer","d-summer","e-autumn","f-autumn","g-winter","h-winter"/);,"c","d"
  ;pres@gsnPanelFigureStrings =(/"a","b","c"/)
  pres@amJust="TopLeft"
  ;gsn_panel(wks,plot,(/1,2/),pres)
  frame(wks)

   end

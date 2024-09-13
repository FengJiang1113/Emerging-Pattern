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

;hor=dtrend_msg_n(hor&time,hor,False,False,0)

copy_VarMeta(ha, hor)
h=hor

;sstannual=month_to_annual(h,1)
;printVarSummary(sstannual)
syear1=1870
eyear1=2022
year1    = ispan(syear1,eyear1,1)


ipo_obs  =wgt_areaave_Wrap(h(:,{-10:10},{170:270}), 1.0,1.0,0)-0.5*(wgt_areaave_Wrap(h(:,{25:45},{140:215}), 1.0,1.0,0)+wgt_areaave_Wrap(h(:,{-50:-15},{150:200}), 1.0,1.0,0))
fcbt =1/(12.*13)
cft1 =ezfftf_n(ipo_obs, 0)
fcft =fspan(0, 0.5,mm/2)
ifcft1=ind(fcft.gt.fcbt)  ;;decadaal
cft1(:,ifcft1)=0.0
ipo_obs:=ezfftb_n(cft1,cft1@xbar,1);
ipo_obs:=dim_standardize(ipo_obs,0)

ipo_obs:=-ipo_obs((1958-1870)*12::)


opendir="/data0/fjiang/data/ocean/oras5/postprocess/"
resultdir="/data0/fjiang/scripts/contrast/"
a = addfile("/data0/fjiang/data/ocean/oras5/oras5.Thermocline.all.int.nc", "r")
b = addfile(opendir+"zs.remap.nc", "r")
c = addfile(opendir+"ms.remap.nc", "r")
d = addfile(opendir+"w.remap.cal.nc", "r")

syear=1958
eyear    = 2022
ymdStrt1 = syear*100+1
ymdLast1 = eyear*100+12
mm=(eyear-syear+1)*12


mm1=eyear-syear+1
time1      = a->time_counter
ymd1      = cd_calendar(time1, -1)
iStrt1    = ind(ymd1.eq.ymdStrt1)
iLast1    = ind(ymd1.eq.ymdLast1)

delete(time1)
delete(ymd1)

lats=-3
latn=3
dd=27
temp      =a->D20(iStrt1:iLast1,{-21:21},{120:290})
tempannual =month_to_annual(temp,1)
uo        =b->sozotaux(iStrt1:iLast1,{-21:21},{120:290})
uoannual  =month_to_annual(uo,1)
vo        =c->sometauy(iStrt1:iLast1,{-21:21},{120:290})
voannual  =month_to_annual(vo,1)

lat=uo&lat
lon=uo&lon
ws=sqrt(uo*uo+vo*vo)
copy_VarMeta(uo,ws)



rc_temp=regCoef_n(ipo_obs,temp,0,0)
copy_VarMeta(temp(0,:,:),rc_temp)

rc_uo=regCoef_n(ipo_obs,uo,0,0)
copy_VarMeta(uo(0,:,:),rc_uo)

rc_vo=regCoef_n(ipo_obs,vo,0,0)
copy_VarMeta(vo(0,:,:),rc_vo)

rc_ws=regCoef_n(ipo_obs,ws,0,0)
copy_VarMeta(uo(0,:,:),rc_ws)


dimst=dimsizes(temp)
confi_temp=student_t(rc_temp@tval,rc_temp@nptxy-2)
confi_temp:=onedtond(confi_temp,(/dimst(1),dimst(2)/))
rc_temp@long_name = "Trend"
copy_VarMeta(temp(0,:,:),rc_temp)
copy_VarMeta(rc_temp,confi_temp)


dimsu=dimsizes(uo)
confi_uo=student_t(rc_uo@tval,rc_uo@nptxy-2)
confi_uo:=onedtond(confi_uo,(/dimsu(1),dimsu(2)/))
rc_uo@long_name = "Trend"
copy_VarMeta(uo(0,:,:),rc_uo)
copy_VarMeta(rc_uo,confi_uo)


confi_vo=student_t(rc_vo@tval,rc_vo@nptxy-2)
confi_vo:=onedtond(confi_vo,(/dimsu(1),dimsu(2)/))
rc_vo@long_name = "Trend"
copy_VarMeta(vo(0,:,:),rc_vo)
copy_VarMeta(rc_vo,confi_vo)

confi_ws=student_t(rc_ws@tval,rc_ws@nptxy-2)
confi_ws:=onedtond(confi_ws,(/dimsu(1),dimsu(2)/))
rc_ws@long_name = "Trend"
copy_VarMeta(vo(0,:,:),rc_ws)
copy_VarMeta(rc_ws,confi_ws)



  wks = gsn_open_wks("pdf",resultdir+"v7-s3");1979-2018")
   gsn_define_colormap(wks,"MPL_RdBu")
   ;gsn_define_colormap(wks,"MPL_RdBu");"CBR_drywet");
  gsn_reverse_colormap(wks);

   res = True
   res@gsnDraw   = False
   res@gsnFrame  = False
   res@gsnAddCyclic      = False
   res@lbLabelBarOn=True
   ;res@vpWidthF =0.9
   res@mpFillOn           = False
   res@mpMinLonF  = 120 ;lon2d(0,0)
   res@mpCenterLonF=180
   res@mpMaxLonF = 290 ;lon2d(0,ilon-1)
   ;res@lbBoxEndCapStyle="";TriangleBothEnds"
   res@pmLabelBarHeightF=0.07/2.
   res@pmLabelBarWidthF=0.55/2.
   res@pmLabelBarOrthogonalPosF=0.1
   ;res@pmLabelBarParallelPosF=0.2
   res@lbLabelStride=2.5
    res@lbBoxLinesOn=False
   res@mpMinLatF            = -20       ; zoom in on map
   res@mpMaxLatF            = 20
   res@vpHeightF = 0.15
   res@vpWidthF  = 0.34
   res@mpShapeMode="FreeAspect"
   res@pmTickMarkDisplayMode="Always"
   res@mpFillOn=True
   res@mpLandFillColor="gray60"
   res@mpOceanFillColor=-1
   res@mpFillDrawOrder="PostDraw"
   res@tmXBMinorOn=False
   ;res@tmYLLabels =(/"30S","20S","10S","0","10N","20N","30N"/)
   res@tmYLMinorOn=False
   res@cnInfoLabelOn=False
   res@mpGeophysicalLineThicknessF=1
   res@mpGeophysicalLineColor="gray60"

   res@pmTickMarkDisplayMode = "Always"
   res@cnFillOn    = True
   res@cnLinesOn = False

   res@tmXBMajorLengthF=0.004
   res@tmYLMajorLengthF=0.004
   res@tmXBMajorOutwardLengthF=0.004
   res@tmYLMajorOutwardLengthF=0.004
;
   res@cnFillOn             = True
   res@cnMaxLevelValF  =1.5
   res@cnMaxLevelValF  =-1.5

   res@cnLinesOn =False
   res@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
   ;res@cnFillColors=(/111,100,96,87,77,67,64,54,44,34,30,20/);
   res@cnFillColors=(/10,15,20,25,30,35,40,50,0,0,80,90,95,100,105,110,115,120/)

;   res@cnLineLabelsOn =True
;   res@cnLineLabelBackgroundColor="-1"
;   res@cnLineLabelInterval =0.5
;   res@cnLineThicknessF=1.5
;   res@cnInfoLabelOn=False
;
;   res@gsnContourNegLineDashPattern=2
;   res@gsnContourPosLineDashPattern=0
;   res@gsnContourZeroLineThicknessF=0
;   res@lbBoxEndCapStyle="TriangleBothEnds"
;   res@lbBoxLinesOn=False
;   res@pmLabelBarHeightF=0.1
   res@gsnLeftStringFontHeightF=0.01
   res@gsnLeftStringOrthogonalPosF=0.01;0.005
   res@gsnRightStringFontHeightF=0.01
   res@gsnRightStringOrthogonalPosF=0.005;0.005

   res@gsnRightString="(~S~o~N~C/65 yrs)"
   res@gsnLeftString="ORAs5 Temp trend 1958-2022  "

;


guess     = 1                ; use zonal means
is_cyclic = True             ; cyclic [global]
nscan     = 1500             ; usually much less than this
eps       = 1.e-2            ; variable dependent
relc      = 0.6              ; relaxation coefficient
opt       = 0                ; not used

poisson_grid_fill(rc_temp, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill(rc_ws, is_cyclic, guess, nscan, eps, relc, opt)
;
;poisson_grid_fill(temp_decadal1, is_cyclic, guess, nscan, eps, relc, opt)
;poisson_grid_fill(temp_decadal2, is_cyclic, guess, nscan, eps, relc, opt)

   plot1=new(6,graphic)
   ;dum1=new(4,graphic)
   res@vpXF=0.12
   res@vpYF=0.9


   res@cnLevels             = 0.1*ispan(-8,8,1 )

   res@tiMainFont=22
   res@tiMainOffsetYF=-0.014
   res@tiMainString="(a) Thermocline P1 minus P2 "
   res@tiMainFontHeightF=0.012
   res@gsnLeftString="";;;';'(a) Thermocline P1 minus P2 "
   res@cnLevels             = 2.*ispan(-8,8,1 )
   res@vpXF=0.12
   res@vpYF=0.9
   res@gsnRightString=""
   ;res@pmLabelBarHeightF=0.16
   ;res@pmLabelBarWidthF=0.03
   res@lbLabelStride=4.
   res@lbBoxLinesOn=True
   res@pmLabelBarOrthogonalPosF=0.15
   res@lbLabelFontHeightF=0.008

   res@tiMainString="(a) IPO-related Thermocline Depth"

   plot1(0) = gsn_csm_contour_map(wks,rc_temp,res)
   res@lbLabelBarOn=False

   ;plot1(1) = gsn_csm_contour_map(wks,temp_decadal2,res)

   res@cnLevels             = 0.001*ispan(-8,8,1 )
   res@vpXF=0.12
   res@vpYF=0.4
   ;res@gsnRightString="(N/m~S~2~N~)"
   res@gsnRightString=""

   res@tiMainString="(b) IPO-related Surface Wind Stress "
   res@lbLabelBarOn=True
   res@gsnRightString=""

   plot1(1) = gsn_csm_contour_map(wks,rc_ws,res)


   ;plot1(3) = gsn_csm_contour_map(wks,wspeed_decadal2,res)

   ;res@vpXF=0.12
   ;res@vpYF=0.44
   ;res@cnLevels             = .03*ispan(-8,8,1 )
   ;res@gsnLeftString="(e) Vertical Current P1 minus P2"
   ;res@gsnRightString="(m/day)"
;
   ;plot1(4) = gsn_csm_contour_map(wks,wo_decadal1,res)
   ;res@vpXF=0.5
   ;res@gsnLeftString="(f) Vertical Current P3 minus P2"
;
   ;plot1(5) = gsn_csm_contour_map(wks,wo_decadal2,res)

;
  res1=True
  res1@gsnDraw=False
  res1@gsnFrame=False
  res1@cnFillOn=True
  res1@cnLinesOn=False
  res1@cnMonoFillPattern=False
  res1@gsnLeftString=""
  res1@gsnRightString=""
  res1@cnLevelSelectionMode="ExplicitLevels"
  res1@cnLevels=(/0,0.1/)
  res1@cnFillPatterns=(/-1,17,-1/)
  res1@cnFillColors=(/"black","black","black"/)
  res1@cnLineLabelsOn=False
  res1@cnInfoLabelOn=False
  res1@lbLabelBarOn=False

  plot2=new(4,graphic)
;
  plot2(0)=gsn_csm_contour(wks,confi_temp,res1)
  plot2(1)=gsn_csm_contour(wks,confi_ws,res1)
  ;plot2(2)=gsn_csm_contour(wks,confi4,res1)
  ;;plot2(3)=gsn_csm_contour(wks,confi4,res1)
;
; ; plot2(1)=gsn_csm_contour(wks,confislp,res1)
;
  overlay(plot1(0),plot2(0))
  overlay(plot1(1),plot2(1))
  ;overlay(plot1(2),plot2(2))
  ;overlay(plot1(3),plot2(3))

  txres=True
  txres@txFontHeightF=0.011
  gsn_text_ndc(wks, "[m]",0.2668,0.343, txres)
  gsn_text_ndc(wks, "[N/m~S~2~N~]",0.768,0.343, txres)


    vcres                         = True             ; plot mods desired
    vcres@gsnDraw          = False           ; don't draw
    vcres@gsnFrame         = False           ; don't advance frame
    ;vcres@vcPositionMode="ArrowTail"
;   vcres@vcGlyphStyle     = "CurlyVector"   ; curly vectors
;   vcres@gsnRightString   = ""            ; turn off right string
;   vcres@gsnLeftString    = ""            ; turn off left string
;   vcres@tiXAxisString    = ""            ; turn off axis label

    vcres@vcLineArrowColor="black"
    vcres@vcRefAnnoOrthogonalPosF = -1.38           ; move ref vector up
    vcres@vcRefAnnoParallelPosF=0.997;-0.0001
    vcres@vcRefLengthF            = 0.045            ; define length of vec ref
    vcres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
    vcres@vcMinDistanceF          = 0.03

    vcres@vcLineArrowThicknessF=1.5
    ;  vcres@vcLineArrowThicknessF      = 1.5
    vcres@vcRefAnnoString2On =True
    vcres@vcRefAnnoString1On =False
    vcres@vcRefAnnoFontHeightF=0.01
    vcres@vcRefAnnoFontThicknessF=0.5
    vcres@vcRefAnnoPerimOn=False
    vcres@vcRefAnnoSide="Bottom"
    vcres@vcRefAnnoBackgroundColor="white"
    vcres@vcRefMagnitudeF  = 0.01    ; define vector ref mag
    vcres@vcRefAnnoString2 = "(10~S~-2~N~ "+"N/m~S~2~N~)"
    vcres@gsnRightString=""
    vcres@gsnLeftString=""

    vector=new(2,graphic)
    vector(0) = gsn_csm_vector(wks,rc_uo,rc_vo,vcres)
    ;vcres@vcRefMagnitudeF  = 0.15    ; define vector ref mag
    ;vcres@vcRefAnnoString2 = "0.15 m/s"
    ;vector(1) = gsn_csm_vector(wks,uo_decadal1,vo_decadal1,vcres)
    overlay(plot1(1),vector(0))
    ;overlay(plot1(3),vector(1))

  ;draw(plot1(0))
  ;draw(plot1(1))
  ;draw(plot1(2))
  ;draw(plot1(3))
  ;draw(plot1(4))
  ;draw(plot1(5))

  ;draw(plot1(3))





;frame(wks)

  pres=True
  pres@gsnPanelCenter=False
  pres@gsnPanelYWhiteSpacePercent=3
  ;pres@gsnPanelLabelBar =True
  ;pres@gsnPanelFigureStrings =(/"a-spring","b-spring","c-summer","d-summer","e-autumn","f-autumn","g-winter","h-winter"/);,"c","d"
  ;pres@gsnPanelFigureStrings =(/"a","b","c"/)
  ;pres@amJust="TopLeft"
  gsn_panel(wks,plot1,(/1,2/),pres)
  ;delete(wks)

   end

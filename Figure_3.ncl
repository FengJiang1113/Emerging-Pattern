load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
resultdir="/data0/fjiang/scripts/contrast/"
;;;;;ORAs5 d20 pattern
opendir1="/data0/fjiang/data/ocean/oras5/"
a = addfile(opendir1+"oras5.d20.all.nc", "r");d20.all.nc", "r");;oras5.d20.all.int.nc

syear1=1980
syear2=1958

eyear    = 2022
ymdStrt1 = syear1*100+1
ymdStrt2 = syear2*100+1

ymdLast1 = eyear*100+12
mm1=(eyear-syear1+1)*12
mm2=(eyear-syear2+1)*12


time1      = a->time_counter
ymd1      = cd_calendar(time1, -1)
iStrt1    = ind(ymd1.eq.ymdStrt1)
iStrt2    = ind(ymd1.eq.ymdStrt2)
iLast1    = ind(ymd1.eq.ymdLast1)

delete(time1)
delete(ymd1)

oras5_d20       =a->D20(:,:,:)

oras5_d20_annual1 =month_to_annual(oras5_d20(iStrt1:iLast1,:,:),1)
oras5_d20_annual2 =month_to_annual(oras5_d20(iStrt2:iLast1,:,:),1)


year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)

dims=dimsizes(oras5_d20_annual1)
rc1    = regCoef_n(year1,oras5_d20_annual1,0,0)
rc2    = regCoef_n(year2,oras5_d20_annual2,0,0)
confi1=student_t(rc1@tval,rc1@nptxy-2)
confi1:=onedtond(confi1,(/dims(1),dims(2)/))
rc1@long_name = "Trend"
rc1:= rc1*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(oras5_d20(0,:,:),rc1)
copy_VarMeta(rc1,confi1)

confi2=student_t(rc2@tval,rc2@nptxy-2)
confi2:=onedtond(confi2,(/dims(1),dims(2)/))
rc2@long_name = "Trend"
rc2:= rc2*10;(eyear-syear2+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(oras5_d20(0,:,:),rc2)
copy_VarMeta(rc2,confi2)


;;;;;
syear    = 1870
eyear    = 1979
ymdStrt1 = syear*100+1
ymdLast1 = eyear*100+12
mm=(eyear-syear+1)*12
opendir="/data0/fjiang/data/ocean/soda224/download/"


b= addfile(opendir+"soda224.d20.all.remap.nc","r");sst.mnmean.nc
time1      = b->time
ymd1      = cd_calendar(time1, -1)
iStrt3    = ind(ymd1.eq.ymdStrt1)
iLast3    = ind(ymd1.eq.ymdLast1)
delete(time1)
delete(ymd1)
soda224_d20   = b->D20;(iStrt3:iLast3,:,:)
;printVarSummary(soda224_d20);
soda224_d20_annual=month_to_annual(soda224_d20,1)
soda224_d20_annual:=soda224_d20_annual(0:108,:,:)
opendir="/data0/fjiang/data/ocean/soda332/"
c= addfile(opendir+"soda332.d20.all.nc","r");sst.mnmean.nc
time1      = c->time
ymd1      = cd_calendar(time1, -1)
iStrt3    = ind(ymd1.eq.ymdStrt1)
iLast3    = ind(ymd1.eq.ymdLast1)
delete(time1)
delete(ymd1)
soda332_d20   = c->D20;(iStrt3:iLast3,:,:)
soda332_d20_annual=month_to_annual(soda332_d20,1)

soda_d20_annual=new((/153,181,360/),float)
soda_d20_annual(1871-1870:1979-1870,:,:)=soda224_d20_annual
soda_d20_annual(1980-1870:2018-1870,:,:)=soda332_d20_annual
patterncor1=new((/153/),float)
patterncor2=new((/153/),float)
patterncor1_ora=new((/153/),float)

rc3_all=new((/106,181,360/),float)
rc4_all=new((/84,181,360/),float)
rc5_all=new((/23,181,360/),float)

lats=-30
late=30
lons=120
lone=270
xx1=year1
printVarSummary(xx1)
do i=0,105
  rc3_all(i,:,:)=regCoef_n(year1,soda_d20_annual(0+i:42+i,:,:),0,0)
  copy_VarMeta(soda_d20_annual(0,:,:),rc3_all(i,:,:))
  patterncor1(i+43)=pattern_cor(rc1({lats:late},{lons:lone}),rc3_all(i,{lats:late},{lons:lone}), 1.0, 0)
end do

do i=0,22
  rc5_all(i,:,:)=regCoef_n(year1,oras5_d20_annual2(0+i:42+i,:,:),0,0)
  copy_VarMeta(soda_d20_annual(0,:,:),rc5_all(i,:,:))

patterncor1_ora(i+43+87)=pattern_cor(rc1({lats:late},{lons:lone}),rc5_all(i,{lats:late},{lons:lone}), 1.0, 0)
end do
print(patterncor1_ora)


xx2=year2
do i=0,83
  rc4_all(i,:,:)=regCoef_n(year2,soda_d20_annual(0+i:64+i,:,:),0,0)
  copy_VarMeta(soda_d20_annual(0,:,:),rc4_all(i,:,:))
  patterncor2(i+65)=pattern_cor(rc2({lats:late},{lons:lone}),rc4_all(i,{lats:late},{lons:lone}), 1.0, 0)
end do

print(patterncor1)
print(patterncor2)


;;;ssh
opendir1="/data0/fjiang/data/ocean/oras5/postprocess/"
ss = addfile(opendir1+"ssh.remap.nc", "r")
syear1=1980
syear2=1958
eyear    = 2022
ymdStrt1 = syear1*100+1
ymdStrt2 = syear2*100+1

ymdLast1 = eyear*100+12
mm1=(eyear-syear1+1)*12
mm2=(eyear-syear2+1)*12


time1      = ss->time_counter
ymd1      = cd_calendar(time1, -1)
iStrt1    = ind(ymd1.eq.ymdStrt1)
iStrt2    = ind(ymd1.eq.ymdStrt2)
iLast1    = ind(ymd1.eq.ymdLast1)

delete(time1)
delete(ymd1)

ssh       =ss->sossheig(:,:,:)
ssh_annual1 =month_to_annual(ssh(iStrt1:iLast1,:,:),1)
ssh_annual2 =month_to_annual(ssh(iStrt2:iLast1,:,:),1)


year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)

rc1_ssh    = regCoef_n(year1,ssh_annual1,0,0)
rc2_ssh    = regCoef_n(year2,ssh_annual2,0,0)

dims_ssh=dimsizes(ssh)
confi1_ssh=student_t(rc1_ssh@tval,rc1_ssh@nptxy-2)
confi1_ssh:=onedtond(confi1_ssh,(/dims_ssh(1),dims_ssh(2)/))
copy_VarMeta(ssh(0,:,:),confi1_ssh)
confi2_ssh=student_t(rc2_ssh@tval,rc2_ssh@nptxy-2)
confi2_ssh:=onedtond(confi2_ssh,(/dims_ssh(1),dims_ssh(2)/))
copy_VarMeta(ssh(0,:,:),confi2_ssh)


rc1_ssh:= rc1_ssh*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(ssh(0,:,:),rc1_ssh)
rc2_ssh:= rc2_ssh*10;(eyear-syear2+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(ssh(0,:,:),rc2_ssh)





;;;copy
syear    = 1870
eyear    = 1979
ymdStrt1 = syear*100+1
ymdLast1 = eyear*100+12
mm=(eyear-syear+1)*12
opendir="/data0/fjiang/data/ocean/soda224/"


s_1= addfile(opendir+"ssh/soda224.ssh.remap.nc","r");sst.mnmean.nc

time1      = s_1->time
ymd1      = cd_calendar(time1, -1)
iStrt3    = ind(ymd1.eq.ymdStrt1)
iLast3    = ind(ymd1.eq.ymdLast1)
delete(time1)
delete(ymd1)

soda224_ssh=s_1->ssh
soda224_ssh_annual=month_to_annual(soda224_ssh,1)
soda224_ssh_annual:=soda224_ssh_annual(0:108,:,:)



opendir="/data0/fjiang/data/ocean/soda332/"
s_2= addfile(opendir+"soda332.ssh.remap.nc","r");sst.mnmean.nc

time1      = s_2->time
ymd1      = cd_calendar(time1, -1)
iStrt3    = ind(ymd1.eq.ymdStrt1)
iLast3    = ind(ymd1.eq.ymdLast1)
delete(time1)
delete(ymd1)

soda332_ssh=s_2->ssh
soda332_ssh_annual=month_to_annual(soda332_ssh,1)

soda_ssh_annual=new((/153,181,360/),float)
soda_ssh_annual(1871-1870:1979-1870,:,:)=soda224_ssh_annual
soda_ssh_annual(1980-1870:2018-1870,:,:)=soda332_ssh_annual
year3    = ispan(1954,2018,1)

rc1_test   = 10*regCoef_n(year3,soda_ssh_annual(1954-1870:2018-1870,:,:),0,0)
copy_VarMeta(soda332_ssh_annual(0,:,:),rc1_test)

patterncor1_ssh=new((/153/),float)
patterncor2_ssh=new((/153/),float)
patterncor1_ora_ssh=new((/153/),float)

rc3_all_ssh=new((/106,181,360/),float)
rc4_all_ssh=new((/84,181,360/),float)
rc5_all_ssh=new((/23,181,360/),float)


xx1=year1
printVarSummary(xx1)
do i=0,105
  rc3_all_ssh(i,:,:)=regCoef_n(year1,soda_ssh_annual(0+i:42+i,:,:),0,0)
  copy_VarMeta(soda_ssh_annual(0,:,:),rc3_all_ssh(i,:,:))
  patterncor1_ssh(i+43)=pattern_cor(rc1_ssh({lats:late},{lons:lone}),rc3_all_ssh(i,{lats:late},{lons:lone}), 1.0, 0)
end do

do i=0,22
  rc5_all_ssh(i,:,:)=regCoef_n(year1,ssh_annual2(0+i:42+i,:,:),0,0)
  copy_VarMeta(soda_ssh_annual(0,:,:),rc5_all_ssh(i,:,:))

patterncor1_ora_ssh(i+43+87)=pattern_cor(rc1_ssh({lats:late},{lons:lone}),rc5_all_ssh(i,{lats:late},{lons:lone}), 1.0, 0)
end do


xx2=year2
do i=0,83
  rc4_all_ssh(i,:,:)=regCoef_n(year2,soda_ssh_annual(0+i:64+i,:,:),0,0)
  copy_VarMeta(soda_ssh_annual(0,:,:),rc4_all_ssh(i,:,:))
  patterncor2_ssh(i+65)=pattern_cor(rc2_ssh({lats:late},{lons:lone}),rc4_all_ssh(i,{lats:late},{lons:lone}), 1.0, 0)
end do



;;;;



  wks = gsn_open_wks("eps",resultdir+"revise-01/r01-F3");1979-2018")
   gsn_define_colormap(wks,"MPL_RdBu")
   ;gsn_define_colormap(wks,"MPL_RdBu");"CBR_drywet");
  gsn_reverse_colormap(wks);


  guess     = 1                ; use zonal means
  is_cyclic = True             ; cyclic [global]
  nscan     = 1500             ; usually much less than this
  eps       = 1.e-2            ; variable dependent
  relc      = 0.6              ; relaxation coefficient
  opt       = 0                ; not used

  poisson_grid_fill(rc1, is_cyclic, guess, nscan, eps, relc, opt)
  poisson_grid_fill(rc2, is_cyclic, guess, nscan, eps, relc, opt)


   res = True
   res@gsnDraw   = False
   res@gsnFrame  = False
   res@gsnAddCyclic      = True
   res@mpShapeMode="FreeAspect"
   res@mpFillOn           = False
   res@mpMinLonF  = 120 ;lon2d(0,0)
   res@mpCenterLonF=180
   res@mpMaxLonF = 290 ;lon2d(0,ilon-1)
   res@lbBoxEndCapStyle="";TriangleBothEnds"
   res@pmLabelBarHeightF=0.07/2.
   res@pmLabelBarWidthF=0.7
   res@vpHeightF = 0.18
   res@vpWidthF  = 0.44
   res@lbLabelStride=2.
    res@lbBoxLinesOn=True
   res@mpMinLatF            = -30       ; zoom in on map
   res@mpMaxLatF            = 30
;
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
   res@tmXBLabelFontHeightF  = 0.008
   res@tmYLLabelFontHeightF  = 0.008

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
   ;res@cnFillColors=(/10,15,20,25,30,35,40,45,50,57,0,0,73,80,85,90,95,100,105,110,115,120/)
   res@cnFillColors=(/12,15,20,25,30,35,40,45,50,60,0,0,70,80,85,90,93,98,103,108,112,115/)

   res@gsnLeftStringFontHeightF=0.01
   res@gsnLeftStringOrthogonalPosF=0.02;0.005
   res@gsnRightStringFontHeightF=0.012
   res@gsnRightStringOrthogonalPosF=-0.003;5;0.005
   res@cnLineLabelsOn=False
   res@gsnRightString="";(~S~o~N~C/decade)"
   res@gsnLeftString=" "
   res@tiMainString="(a) ORAs5 Thermocline Trend 1980-2022  "
   res@tiMainFont=22
   res@tiMainOffsetYF=-0.014
   res@tiMainFontHeightF=0.0115
   plot1=new(8,graphic)
   ;dum1=new(4,graphic)
   res@vpXF=0.06
   res@vpYF=0.9

             ; change aspect ratio of plot

   res@cnLevels             = 1.*ispan(-10,10,1 )

   res@lbLabelBarOn=False

   plot1(0) = gsn_csm_contour_map(wks,rc1,res)
   res@lbLabelBarOn=True
   res@pmLabelBarOrthogonalPosF=0.2
   res@pmLabelBarParallelPosF=-0.11
   res@tiMainString="(b) ORAs5 Thermocline Trend 1958-2022  "
   res@cnLevels             = 1.*ispan(-10,10,1 )
   res@vpXF=0.5
   res@vpYF=0.9
   plot1(1) = gsn_csm_contour_map(wks,rc2,res)



  res1=True
  res1@gsnDraw=False
  res1@gsnFrame=False
  res1@cnFillOn=True
  res1@cnLinesOn=False
  res1@cnMonoFillPattern=False
  res1@gsnLeftString=""
  res1@gsnRightString=""
  res1@cnLevelSelectionMode="ExplicitLevels"
  res1@cnLevels=(/0,0.05/)
  res1@cnFillPatterns=(/-1,17,-1/)
  res1@cnFillColors=(/"black","black","black"/)
  res1@cnLineLabelsOn=False
  res1@cnInfoLabelOn=False
  res1@lbLabelBarOn=False
  plot2=new(4,graphic)
;
  plot2(0)=gsn_csm_contour(wks,confi1,res1)
  plot2(1)=gsn_csm_contour(wks,confi2,res1)

  overlay(plot1(0),plot2(0))
  overlay(plot1(1),plot2(1))
;;  overlay(plot1(2),plot2(2))
;;  overlay(plot1(3),plot2(3))

;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
delete(res)
res                 = True
res@gsnDraw         = False
res@gsnFrame        = False
res@vpHeightF = 0.18
res@vpWidthF  = 0.44
res@trXMinF   =  1870.
res@trXMaxF   =  2022.
res@trYMaxF   = 1.
res@trYMinF   = -1.
res@tiYAxisString = " "
res@tiXAxisString = "End year of the trend"
res@tiXAxisFontHeightF=0.01

res@tiMainString  = " "
res@tmXBMajorThicknessF     = 1.5
res@tmYLMajorThicknessF     = 1.5
res@tmXBMinorThicknessF     = 1.0
res@tmYLMinorThicknessF     = 1.0
res@tmYLMinorOn=False
res@tmYLValues    = 0.5*ispan(-5,5,1)
res@tmYLLabels    = 0.5*ispan(-5,5,1);
res@tmXBValues    = ispan(1870,2020,10)
res@tmXBLabels    = ispan(1870,2020,10); (/"1960","1965","1970","1975","1980","1985","1985","1990","1995","2000","2005","2010","2015","2020"/)
;
res@gsnYRefLine=0.
res@tmXBMajorLengthF=0.006
res@tmYLMajorLengthF=0.004
res@tmXBMajorOutwardLengthF=0.006
res@tmYLMajorOutwardLengthF=0.004
res@tmXBMinorLengthF        = 0.004
res@tmYLMinorLengthF        = 0.004
res@tmXBMinorOn=True
;res@tmXBMinorValues=ispan(1958,2022,1)

res@gsnLeftString=" "
res@gsnLeftStringFontHeightF= 0.012
res@gsnLeftStringOrthogonalPosF=0.00;5
res@tmXBOn = True
;res@tmYLOn = True
;res@tmYLLabelsOn = True
;res@tmYLMinorOn = True
;res@tmYLLabelFontColor = "black"
res@tiYAxisString         = ""
res@tiYAxisFontColor      = "black"
res@tmYLOn    = True
res@tmYLMode          = "Explicit"         ; explicit labels
res@tmYLLabelsOn=True

res@tmXBLabelFontHeightF  = 0.008
res@tmYLLabelFontHeightF  = 0.008
res@tiYAxisSide           = "Left"
res@tiYAxisFontHeightF    = 0.01
res@xyDashPattern    = 0
res@xyLineColor      = "orangered3"
res@tmXBBorderOn=True
res@tiYAxisFontHeightF=0.012
;res_low=res
;res_low@xyLineThicknessF=0.8
;res_low@xyDashPattern    = 0
;
;res@xyLineOpacityF=0.8
;
;res_trend=res_low
;res_trend@xyLineThicknessF=2.
res@tmXTOn=False
res@tmYROn=False
;
;
res@vpXF=0.068
res@vpYF=0.63
;print(ts_sst2)

;plot22 = new(4,graphic)
;plot33 = new(4,graphic)
res@tiMainFont=22
res@tiMainOffsetYF=-0.01
res@tiMainFontHeightF=0.0115
res@tiMainString="(e) Pattern Correlation with tropical Pacific Trend 1980-2022"
res@gsnRightString=""
res@xyLineThicknessF = 2.

plot1(4) = gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1,res)
res@tiMainString="(f) Pattern Correlation with tropical Pacific Trend 1958-2022"

plot1(5) = gsn_csm_xy(wks,ispan(1870,2022,1),patterncor2,res)
res@xyLineColor      = "orangered3"
res@gsnLeftString=""

res@xyLineThicknessF = 1.5
res@xyDashPattern    = 2

plota=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1_ora,res)
;res@xyLineColor      = "deepskyblue3"
;res@xyDashPattern    = 1
;
;plotb=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1_b,res)
;res@xyLineColor      = "gray40"
;res@xyDashPattern    = 2
;
;plotc=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor2_a,res)
;res@xyLineColor      = "deepskyblue3"
;res@xyDashPattern    = 1
;
;plotd=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor2_b,res)
overlay(plot1(4),plota)

delete(res)



pgres=True
pgres@gsLineColor      = "orangered3"
pgres@gsLineThicknessF = 2.0
pgres@tfPolyDrawOrder ="PostDraw"
pgres@xyDashPattern    = 0


txres=True
txres@txFontHeightF=0.011
gsn_text_ndc(wks, "[m/decade]",0.51,0.635, txres)


pgres=True
pgres@gsMarkerIndex      = 12

pgres@gsMarkerSizeF = 4.
pgres@tfPolyDrawOrder ="PostDraw"
pgres@gsMarkerColor="orangered3"
;dum_gon_1 = gsn_add_polymarker (wks, plot1(2),151.2+1870,patterncor1_oras5, pgres)
;pgres@gsMarkerColor="gray40"
;pgres@gsMarkerIndex      = 16
;dum_gon_2 = gsn_add_polymarker (wks, plot1(2),151.2+1870,patterncor1_a_oras5, pgres)
;pgres@gsMarkerColor="deepskyblue3"
;pgres@gsMarkerIndex      = 14
;dum_gon_3 = gsn_add_polymarker (wks, plot1(2),151.2+1870,patterncor1_b_oras5, pgres)
;pgres@gsMarkerColor="orangered3"
;pgres@gsMarkerIndex      = 12
;dum_gon_4 = gsn_add_polymarker (wks, plot1(3),151.2+1870,patterncor2_oras5, pgres)
;pgres@gsMarkerColor="gray40"
;pgres@gsMarkerIndex      = 16
;dum_gon_5 = gsn_add_polymarker (wks, plot1(3),151.2+1870,patterncor2_a_oras5, pgres)
;pgres@gsMarkerColor="deepskyblue3"
;pgres@gsMarkerIndex      = 14
;dum_gon_6 = gsn_add_polymarker (wks, plot1(3),151.2+1870,patterncor2_b_oras5, pgres)

pgres@gsMarkerColor="orangered3"
pgres@gsMarkerIndex      = 16
pgres@gsMarkerSizeF = 8.

;
print(patterncor1)
print(patterncor2)
txres@txFontHeightF=0.0105


;;;;;ws
resultdir="/data0/fjiang/scripts/contrast/"
;;;;;ORAs5 wspeed pattern
opendir1="/data0/fjiang/data/ocean/oras5/postprocess/"
a = addfile(opendir1+"zs.remap.nc", "r")
b = addfile(opendir1+"ms.remap.nc", "r")

syear1=1980
syear2=1958

eyear    = 2022
ymdStrt1 = syear1*100+1
ymdStrt2 = syear2*100+1

ymdLast1 = eyear*100+12
mm1=(eyear-syear1+1)*12
mm2=(eyear-syear2+1)*12


time1      = a->time_counter
ymd1      = cd_calendar(time1, -1)
iStrt1    = ind(ymd1.eq.ymdStrt1)
iStrt2    = ind(ymd1.eq.ymdStrt2)
iLast1    = ind(ymd1.eq.ymdLast1)

delete(time1)
delete(ymd1)

oras5_taux      =a->sozotaux
oras5_tauy      =b->sometauy
oras5_wspeed=sqrt(oras5_taux*oras5_taux+oras5_tauy*oras5_tauy)
copy_VarMeta(oras5_taux,oras5_wspeed)
oras5_wspeed_annual1 =month_to_annual(oras5_wspeed(iStrt1:iLast1,:,:),1)
oras5_wspeed_annual2 =month_to_annual(oras5_wspeed(iStrt2:iLast1,:,:),1)

oras5_taux_annual1 =month_to_annual(oras5_taux(iStrt1:iLast1,:,:),1)
oras5_taux_annual2 =month_to_annual(oras5_taux(iStrt2:iLast1,:,:),1)


oras5_tauy_annual1 =month_to_annual(oras5_tauy(iStrt1:iLast1,:,:),1)
oras5_tauy_annual2 =month_to_annual(oras5_tauy(iStrt2:iLast1,:,:),1)

year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)

dims=dimsizes(oras5_wspeed_annual1)
rc1    = regCoef_n(year1,oras5_wspeed_annual1,0,0)
rc2    = regCoef_n(year2,oras5_wspeed_annual2,0,0)

rc1_taux    = regCoef_n(year1,oras5_taux_annual1,0,0)
rc2_taux    = regCoef_n(year2,oras5_taux_annual2,0,0)

rc1_tauy    = regCoef_n(year1,oras5_tauy_annual1,0,0)
rc2_tauy    = regCoef_n(year2,oras5_tauy_annual2,0,0)

;rc1_ssh:=rc1_test

rc1_taux:= rc1_taux*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
rc2_taux:= rc2_taux*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
rc1_tauy:= rc1_tauy*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
rc2_tauy:= rc2_tauy*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(oras5_wspeed(0,:,:),rc1_taux)
copy_VarMeta(oras5_wspeed(0,:,:),rc2_taux)
copy_VarMeta(oras5_wspeed(0,:,:),rc1_tauy)
copy_VarMeta(oras5_wspeed(0,:,:),rc2_tauy)


delete(confi1)
delete(confi2)
confi1=student_t(rc1@tval,rc1@nptxy-2)
confi1:=onedtond(confi1,(/dims(1),dims(2)/))
rc1@long_name = "Trend"
rc1:= rc1*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(oras5_wspeed(0,:,:),rc1)
copy_VarMeta(rc1,confi1)

confi2=student_t(rc2@tval,rc2@nptxy-2)
confi2:=onedtond(confi2,(/dims(1),dims(2)/))
rc2@long_name = "Trend"
rc2:= rc2*10;(eyear-syear2+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(oras5_wspeed(0,:,:),rc2)
copy_VarMeta(rc2,confi2)


;;;;;
syear    = 1870
eyear    = 1979
ymdStrt1 = syear*100+1
ymdLast1 = eyear*100+12
mm=(eyear-syear+1)*12
opendir="/data0/fjiang/data/ocean/soda224/"


b_x= addfile(opendir+"taux/soda224.taux.remap.nc","r");sst.mnmean.nc
b_y= addfile(opendir+"tauy/soda224.tauy.remap.nc","r");sst.mnmean.nc

time1      = b_x->time
ymd1      = cd_calendar(time1, -1)
iStrt3    = ind(ymd1.eq.ymdStrt1)
iLast3    = ind(ymd1.eq.ymdLast1)
delete(time1)
delete(ymd1)

soda224_taux=b_x->taux
soda224_tauy=b_y->tauy
soda224_wspeed=sqrt(soda224_taux*soda224_taux+soda224_tauy*soda224_tauy)
copy_VarMeta(soda224_taux,soda224_wspeed)
printVarSummary(soda224_wspeed)

soda224_wspeed_annual=month_to_annual(soda224_wspeed,1)
soda224_wspeed_annual:=soda224_wspeed_annual(0:108,:,:)
opendir="/data0/fjiang/data/ocean/soda332/"
c_x= addfile(opendir+"soda332.taux.remap.nc","r");sst.mnmean.nc
c_y= addfile(opendir+"soda332.tauy.remap.nc","r");sst.mnmean.nc

time1      = c_x->time
ymd1      = cd_calendar(time1, -1)
iStrt3    = ind(ymd1.eq.ymdStrt1)
iLast3    = ind(ymd1.eq.ymdLast1)
delete(time1)
delete(ymd1)

soda332_taux=c_x->taux
soda332_tauy=c_y->tauy
soda332_wspeed=sqrt(soda332_taux*soda332_taux+soda332_tauy*soda332_tauy)
copy_VarMeta(soda332_taux,soda332_wspeed)
soda332_wspeed_annual=month_to_annual(soda332_wspeed,1)


soda_wspeed_annual=new((/153,181,360/),float)
soda_wspeed_annual(1871-1870:1979-1870,:,:)=soda224_wspeed_annual
soda_wspeed_annual(1980-1870:2018-1870,:,:)=soda332_wspeed_annual
patterncor1=new((/153/),float)
patterncor2=new((/153/),float)
patterncor1_ora=new((/153/),float)

rc3_all=new((/106,181,360/),float)
rc4_all=new((/84,181,360/),float)
rc5_all=new((/23,181,360/),float)


xx1=year1
printVarSummary(xx1)
do i=0,105
  rc3_all(i,:,:)=regCoef_n(year1,soda_wspeed_annual(0+i:42+i,:,:),0,0)
  copy_VarMeta(soda_wspeed_annual(0,:,:),rc3_all(i,:,:))
  patterncor1(i+43)=pattern_cor(rc1({lats:late},{lons:lone}),rc3_all(i,{lats:late},{lons:lone}), 1.0, 0)
end do

do i=0,22
  rc5_all(i,:,:)=regCoef_n(year1,oras5_wspeed_annual2(0+i:42+i,:,:),0,0)
  copy_VarMeta(soda_wspeed_annual(0,:,:),rc5_all(i,:,:))

patterncor1_ora(i+43+87)=pattern_cor(rc1({lats:late},{lons:lone}),rc5_all(i,{lats:late},{lons:lone}), 1.0, 0)
end do



xx2=year2
do i=0,83
  rc4_all(i,:,:)=regCoef_n(year2,soda_wspeed_annual(0+i:64+i,:,:),0,0)
  copy_VarMeta(soda_wspeed_annual(0,:,:),rc4_all(i,:,:))
  patterncor2(i+65)=pattern_cor(rc2({lats:late},{lons:lone}),rc4_all(i,{lats:late},{lons:lone}), 1.0, 0)
end do


   res = True
   res@gsnDraw   = False
   res@gsnFrame  = False
   res@gsnAddCyclic      = True
   res@mpShapeMode="FreeAspect"
   ;res@vpWidthF =0.9
   res@mpFillOn           = False
   res@mpMinLonF  = 120 ;lon2d(0,0)
   res@mpCenterLonF=180
   res@mpMaxLonF = 290 ;lon2d(0,ilon-1)
   res@lbBoxEndCapStyle="";TriangleBothEnds"
   res@pmLabelBarHeightF=0.07/2.
   res@pmLabelBarWidthF=0.7
   res@vpHeightF = 0.18
   res@vpWidthF  = 0.44
   res@lbLabelStride=2.
    res@lbBoxLinesOn=True
   res@mpMinLatF            = -30       ; zoom in on map
   res@mpMaxLatF            = 30
;
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
   res@tmXBLabelFontHeightF  = 0.008
   res@tmYLLabelFontHeightF  = 0.008

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

   res@cnFillColors=(/12,15,20,25,30,35,40,45,50,60,0,0,70,80,85,90,93,98,103,108,112,115/)

   res@gsnLeftStringFontHeightF=0.01
   res@gsnLeftStringOrthogonalPosF=0.005;0.005
   res@gsnRightStringFontHeightF=0.012
   res@gsnRightStringOrthogonalPosF=-0.003;5;0.005
   res@cnLineLabelsOn=False
   res@gsnRightString="";(~S~o~N~C/decade)"
   res@gsnLeftString="";
   res@tiMainString="(c) ORAs5 Wind Stress & SSH Trend 1980-2022  "
   res@tiMainFont=22
   res@tiMainOffsetYF=-0.014
   res@tiMainFontHeightF=0.0115

   ;dum1=new(4,graphic)
   res@vpXF=0.06
   res@vpYF=0.9

   res@cnLevels             = 0.004*ispan(-10,10,1 )

   res@lbLabelBarOn=False


     poisson_grid_fill(rc1_ssh, is_cyclic, guess, nscan, eps, relc, opt)
     poisson_grid_fill(rc2_ssh, is_cyclic, guess, nscan, eps, relc, opt)


   plot1(2) = gsn_csm_contour_map(wks,rc1_ssh,res)
   res@lbLabelBarOn=True
   res@pmLabelBarOrthogonalPosF=0.2
   res@pmLabelBarParallelPosF=-0.11
   res@tiMainString="(d) ORAs5 Wind Stress Trend & SSH 1958-2022  "
   res@cnLevels             = 0.004*ispan(-10,10,1 )
   res@vpXF=0.5
   res@vpYF=0.9
   plot1(3) = gsn_csm_contour_map(wks,rc2_ssh,res)


  res1=True
  res1@gsnDraw=False
  res1@gsnFrame=False
  res1@cnFillOn=True
  res1@cnLinesOn=False
  res1@cnMonoFillPattern=False
  res1@gsnLeftString=""
  res1@gsnRightString=""
  res1@cnLevelSelectionMode="ExplicitLevels"
  res1@cnLevels=(/0,0.05/)
  res1@cnFillPatterns=(/-1,17,-1/)
  res1@cnFillColors=(/"black","black","black"/)
  res1@cnLineLabelsOn=False
  res1@cnInfoLabelOn=False
  res1@lbLabelBarOn=False
;  plot2=new(2,graphic)
;
;  plot2(0)=gsn_csm_contour(wks,confi1,res1)
;  plot2(1)=gsn_csm_contour(wks,confi2,res1)
  plot2(2)=gsn_csm_contour(wks,confi1_ssh,res1)
  plot2(3)=gsn_csm_contour(wks,confi2_ssh,res1)
;;
;;;  plot2(1)=gsn_csm_contour(wks,confislp,res1)

  overlay(plot1(2),plot2(2))
  overlay(plot1(3),plot2(3))

;
;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
delete(res)
res                 = True
res@gsnDraw         = False
res@gsnFrame        = False
res@trXMinF   =  1870.
res@trXMaxF   =  2022.
res@trYMaxF   = 1.
res@trYMinF   = -1.
res@tiYAxisString = " "
res@tiXAxisString = "End year of the trend"
res@tiXAxisFontHeightF=0.01
res@tiMainString  = " "
res@tmXBMajorThicknessF     = 1.5
res@tmYLMajorThicknessF     = 1.5
res@tmXBMinorThicknessF     = 1.0
res@tmYLMinorThicknessF     = 1.0
res@tmYLMinorOn=False
res@tmYLValues    = 0.5*ispan(-5,5,1)
res@tmYLLabels    = 0.5*ispan(-5,5,1);
res@tmXBValues    = ispan(1870,2020,10)
res@tmXBLabels    = ispan(1870,2020,10); (/"1960","1965","1970","1975","1980","1985","1985","1990","1995","2000","2005","2010","2015","2020"/)
;
res@gsnYRefLine=0.
res@tmXBMajorLengthF=0.006
res@tmYLMajorLengthF=0.004
res@tmXBMajorOutwardLengthF=0.006
res@tmYLMajorOutwardLengthF=0.004
res@tmXBMinorLengthF        = 0.004
res@tmYLMinorLengthF        = 0.004
res@tmXBMinorOn=True
;res@tmXBMinorValues=ispan(1958,2022,1)

res@gsnLeftString=" "
res@gsnLeftStringFontHeightF= 0.01
res@gsnLeftStringOrthogonalPosF=0.00;5
res@tmXBOn = True
;res@tmYLOn = True
;res@tmYLLabelsOn = True
;res@tmYLMinorOn = True
;res@tmYLLabelFontColor = "black"
res@tiYAxisString         = ""
res@tiYAxisFontColor      = "black"
res@tmYLOn    = True
res@tmYLMode          = "Explicit"         ; explicit labels
res@tmYLLabelsOn=True

res@tmXBLabelFontHeightF  = 0.008
res@tmYLLabelFontHeightF  = 0.008
res@tiYAxisSide           = "Left"
res@tiYAxisFontHeightF    = 0.01
res@xyDashPattern    = 0
res@tmXBBorderOn=True
res@tiYAxisFontHeightF=0.012
res@xyLineColor      = "darkseagreen"

res@tmXTOn=False
res@tmYROn=False
;
;
res@vpXF=0.068
res@vpYF=0.63
;print(ts_sst2)

;plot22 = new(4,graphic)
;plot33 = new(4,graphic)
res@gsnLeftString=""
res@tiMainString=""
res@gsnRightString=""
res@xyLineThicknessF = 2.

res@tiMainFont=22
res@tiMainOffsetYF=-0.01
res@tiMainFontHeightF=0.0115

plotb = gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1,res)

plotc = gsn_csm_xy(wks,ispan(1870,2022,1),patterncor2,res)

res@xyLineColor      = "deepskyblue3"

plotd = gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1_ssh,res)
plote = gsn_csm_xy(wks,ispan(1870,2022,1),patterncor2_ssh,res)

res@gsnLeftString=""

res@xyLineThicknessF = 1.5
res@xyDashPattern    = 2
res@xyLineColor      = "darkseagreen"

plotf=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1_ora,res)
res@xyLineColor      = "deepskyblue3"

plotg=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1_ora_ssh,res)


overlay(plot1(4),plotb)
overlay(plot1(5),plotc)
overlay(plot1(4),plotd)
overlay(plot1(5),plote)
overlay(plot1(4),plotf)
overlay(plot1(4),plotg)


;dum_box3 = gsn_add_polyline   (wks, plot1(4),(/lons,lone,lone,lons,lons/),(/lats,lats,late,late,lats/), pgres)
;dum_box4 = gsn_add_polyline   (wks, plot1(5),(/lons,lone,lone,lons,lons/),(/lats,lats,late,late,lats/), pgres)

txres=True
txres@txFontHeightF=0.011
gsn_text_ndc(wks, "[m/decade]",0.51,0.332, txres)


pgres=True
pgres@gsMarkerIndex      = 12

pgres@gsMarkerSizeF = 4.
pgres@tfPolyDrawOrder ="PostDraw"
pgres@gsMarkerColor="orangered3"

  txres@txFontHeightF=0.0105


  vcres                         = True             ; plot mods desired
  vcres@gsnDraw          = False           ; don't draw
  vcres@gsnFrame         = False           ; don't advance frame
  ;vcres@vcPositionMode="ArrowTail"
  ;vcres@vcLevelSpacingF  = 3
  vcres@gsnRightString   = " "             ; turn off right string
  vcres@gsnLeftString    = " "             ; turn off left string
  vcres@tiXAxisString    = " "             ; turn off axis label
  vcres@vcLineArrowColor="black"
  vcres@vcRefLengthF            = 0.045            ; define length of vec ref
  ;vcres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  vcres@vcMinDistanceF          = 0.035
  vcres@vcLineArrowThicknessF=1.3
  ;  vcres@vcLineArrowThicknessF      = 1.5
  vcres@vcRefAnnoString2On =True
  vcres@vcRefAnnoString1On =False
  vcres@vcRefAnnoFontHeightF=0.009
  vcres@vcRefAnnoFontThicknessF=0.5
  vcres@vcRefAnnoPerimOn=False
  vcres@vcRefAnnoSide="Bottom"
  vcres@vcRefAnnoBackgroundColor="white"
  vcres@vcRefMagnitudeF  = 0.005    ; define vector ref mag
  vcres@vcRefAnnoString2 =    "(0.005 "+"N/m~S~2~N~/decade)"
  vcres@vcGlyphStyle         = "LineArrow"


  vector=new(2,graphic)

  vcres@vcRefAnnoOrthogonalPosF = -1.34+0.19+1.2      ; move ref vector up
  vcres@vcRefAnnoParallelPosF=0.15;-0.0001

  vector(0) = gsn_csm_vector(wks,rc1_taux,rc1_tauy,vcres)
  vcres@vcRefAnnoParallelPosF=0.997;-0.0001

  vcres@vcRefAnnoOrthogonalPosF = -1.64+0.2+1.2      ; move ref vector up
  vector(1) = gsn_csm_vector(wks,rc2_taux,rc2_tauy,vcres)

  overlay(plot1(2),vector(0))
  overlay(plot1(3),vector(1))


   lgres=True
   amres=True
   lgres@lgLineColors=(/"darkseagreen","darkseagreen","deepskyblue3","deepskyblue3","orangered3","orangered3"/);"deepskyblue3",
   lgres@lgLeftMarginF=0.1
   lgres@lgRightMarginF=0.1

   lgres@lgDashIndexes=(/2,0,2,0,2,0/);1,
   lgres@lgLineThicknesses=(/1.,2.5,1.,2.5,1.,2.5/);1.5,
   lgres@lgMonoLineThickness=False
   lgres@vpWidthF=0.14
   lgres@vpHeightF=0.085;55
   lgres@lgLabelFontHeightF=0.15
   lgres@lgLineDashSegLenF=0.04
   lgres@lgPerimThicknessF=2.
   lgres@lgPerimOn=False
   lgres@lgPerimFill="SolidFill"
   lgres@lgPerimFillColor="white"
   ;lgres@lgLabelOffsetF=-0.1
   lgres@lbBoxMajorExtentF=0.3


   lbid1=gsn_create_legend(wks,6,(/"  ORAs5 Zonal Wind Stress","  SODA  Zonal Wind Stress","  ORAs5 SSH","  SODA SSH","  ORAs5 Thermocline Depth","  SODA Thermocline Depth"/),lgres);"  EN4",
   lgres@vpHeightF=0.085/2;55
   lgres@lgDashIndexes:=(/0,0,0/);1,
   lgres@lgLineThicknesses:=(/2.5,2.5,2.5/);1.5,
   lgres@lgLineColors:=(/"darkseagreen","deepskyblue3","orangered3"/);"deepskyblue3",
   lbid2=gsn_create_legend(wks,3,(/"  SODA  Zonal Wind Stress","  SODA SSH","  SODA Thermocline Depth"/),lgres);"  EN4",

   amres@amParallelPosF=-0.33
   amres@amOrthogonalPosF=-0.26
   annoid1=gsn_add_annotation(plot1(4),lbid1,amres)
   amres@amOrthogonalPosF=-0.33
   annoid2=gsn_add_annotation(plot1(5),lbid2,amres)


  pres=True
  pres@gsnPanelCenter=False
  pres@gsnPanelYWhiteSpacePercent=20
  pres@gsnPanelXWhiteSpacePercent=3.5

  ;pres@gsnPanelLabelBar =True
  ;pres@gsnPanelFigureStrings =(/"a-spring","b-spring","c-summer","d-summer","e-autumn","f-autumn","g-winter","h-winter"/);,"c","d"
  ;pres@gsnPanelFigureStrings =(/"a","b","c"/)
  ;pres@amJust="TopLeft"
  gsn_panel(wks,plot1,(/3,2/),pres)
  ;delete(wks)

   end

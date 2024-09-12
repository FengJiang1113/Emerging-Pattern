load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

function arrow(wks,plt,xpts,ypts,res)
local x0, x1, y0, y1, xp1, xp2, yp1, yp2, xp1n, xp2n, yp1n, yp2n, pi, \
rad2deg, theta, asp, fhead, ahead, phi, resp, lambda1, lambda2, arr_str, \
x1an, ya1n, xa1, ya1, xpts1, ypts1, xa2n, ya2n, xa2, ya2, xpts2, ypts2
begin

 res2 = res    ; make a copy so we can modify it if needed.

 getvalues plt
  "trXMinF": x0
  "trXMaxF": x1
  "trYMinF": y0
  "trYMaxF": y1
 end getvalues

 arr_str = unique_string("arr")   ; need to make sure this is unique
 plt@$arr_str$ = gsn_add_polyline(wks,plt,xpts,ypts,res2)

; determines how broad should the arrow be
 asp = 0.4

 xp1 = xpts(0)
 xp2 = xpts(1)
 yp1 = ypts(0)
 yp2 = ypts(1)

 xp1n = (xp1 - x0)/(x1-x0)
 xp2n = (xp2 - x0)/(x1-x0)
 yp1n = (yp1 - y0)/(y1-y0)
 yp2n = (yp2 - y0)/(y1-y0)

 pi = get_pi("float")
 rad2deg = get_r2d("float")

 theta = atan2((yp2n-yp1n),(xp2n-xp1n))

 if(theta.lt.0 .and. yp2.gt.yp1)then
  theta = pi + theta
 else
  if(theta.gt.0 .and. yp2.lt.yp1)then
   theta = pi + theta
  end if
 end if

; length of the line of arrowhead as a fraction of the length of the given line
 fhead = 0.2
 ahead = sqrt((xp2n-xp1n)*(xp2n-xp1n) + (yp2n-yp1n)*(yp2n-yp1n))*fhead
 phi = atan(asp)

 lambda1 = theta + phi
 xa1n = xp2n - ahead*cos(lambda1)
 ya1n = yp2n - ahead*sin(lambda1)
 xa1 = x0 + xa1n*(x1-x0)
 ya1 = y0 + ya1n*(y1-y0)

;---Don't add line labels to arrowhead.
 if(isatt(res2,"gsLineLabelString")) then
    delete(res2@gsLineLabelString)
 end if

 xpts1 = (/xa1,xp2/)
 ypts1 = (/ya1,yp2/)
 arr_str = unique_string("arr")   ; need to make sure this is unique
 plt@$arr_str$ = gsn_add_polyline(wks,plt,xpts1,ypts1,res2)

 lambda2 = theta - phi
 xa2n = xp2n - ahead*cos(lambda2)
 ya2n = yp2n - ahead*sin(lambda2)
 xa2 = x0 + xa2n*(x1-x0)
 ya2 = y0 + ya2n*(y1-y0)
 xpts2 = (/xa2,xp2/)
 ypts2 = (/ya2,yp2/)
 arr_str = unique_string("arr")   ; need to make sure this is unique
 plt@$arr_str$ = gsn_add_polyline(wks,plt,xpts2,ypts2,res2)

 return (plt)
end



begin

syear    = 1870
eyear    = 2022
ymdStrt1 = syear*100+1
ymdLast1 = eyear*100+12
ymdbetween=2000*100+1
mm=(eyear-syear+1)*12
mm1=(2018-2000)*12;(2000-1980)*12;

opendir="/data0/fjiang/data/ocean/"
resultdir="/data0/fjiang/scripts/contrast/revise-01/"


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
;ha:=lonPivot(ha,180.5)

hor= where(abs(ha).ge.50,ha@_FillValue,ha)
hor= where(ismissing(hor),0,hor)
copy_VarMeta(ha, hor)
h=hor

sstannual=month_to_annual(h,1)
printVarSummary(sstannual)


sstannual_rm=sstannual
sstannual_30=dim_avg_n_Wrap(sstannual(0:49,:,:),0)
do i=0,eyear-syear
sstannual_rm(i,:,:)=sstannual(i,:,:)-sstannual_30;dim_rmvmean_n_Wrap(sstannual,0)
end do
copy_VarMeta(sstannual,sstannual_rm)
sstannual=sstannual_rm

latS=-5
latN=5

sstannual_ep=wgt_areaave_Wrap(sstannual(:,{latS:latN},{190:270}),1.0,1.0,0)
sstannual_cp=wgt_areaave_Wrap(sstannual(:,{latS:latN},{140:170}),1.0,1.0,0)
sstannual_gr=wgt_areaave_Wrap(sstannual(:,{latS:latN},{140:170}),1.0,1.0,0)

sstannual_gr=sstannual_ep-sstannual_cp


sstannual_ep_10=runave_n_Wrap(sstannual_ep,15,0,0)
sstannual_cp_10=runave_n_Wrap(sstannual_cp,15,0,0)
sstannual_gr_10=runave_n_Wrap(sstannual_gr,15,0,0)


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

fname1 = "/data0/fjiang/scripts/contrast/figuredata/IPO.nc";;eq.
fout1 = addfile(fname1, "r")
rc1_ipo=fout1->rc1_ipo
patterncor1=new((/153/),float)
patterncor1_ipo=new((/153/),float)

;patterncor1_a=new((/153/),float)
;patterncor1_b=new((/153/),float)
patterncor2=new((/153/),float)

patterncor2_ipo=new((/153/),float)
;patterncor2_a=new((/153/),float)
;patterncor2_b=new((/153/),float)


lats=-30
late=30
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
;;;;
  patterncor1_ipo(i+42)=-1*pattern_cor(rc1_ipo({lats:late},{lons:lone}),rc1_all(i,{lats:late},{lons:lone}), 1.0, 0)

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
  ;;;;
  patterncor2_ipo(i+64)=-1*pattern_cor(rc1_ipo({lats:late},{lons:lone}),rc2_all(i,{lats:late},{lons:lone}), 1.0, 0)

end do

print(patterncor1)
print(patterncor2)


rc1_all=rc1_all*10
copy_VarMeta(rc2_all,rc1_all)
rc2_all=rc2_all*10
copy_VarMeta(rc1_all,rc2_all)


syear    = 1870
eyear    = 2022
ymdStrt1 = syear*100+1
ymdLast1 = eyear*100+12
ymdbetween=2000*100+1
mm=(eyear-syear+1)*12
mm1=(2018-2000)*12;(2000-1980)*12;

opendir="/data0/fjiang/data/ocean/"
resultdir="/data0/fjiang/scripts/contrast/"


a= addfile(opendir+"sst.mnmean.ersstv5.remap.nc","r")
time1      = a->time
ymd1      = cd_calendar(time1, -1)
iStrt3    = ind(ymd1.eq.ymdStrt1)
iLast3    = ind(ymd1.eq.ymdLast1)
ibe=ind(ymd1.eq.ymdbetween)
delete(time1)
delete(ymd1)
delete(ha)
delete(sstannual)
ha   = a->sst(iStrt3:iLast3,::-1,:)


sstannual=month_to_annual(ha,1)
printVarSummary(sstannual)


syear1=1980
syear2=1958
year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)


rc1_c    = regCoef_n(year1,sstannual(syear1-syear:eyear-syear,:,:),0,0)
rc2_c    = regCoef_n(year2,sstannual(syear2-syear:eyear-syear,:,:),0,0)



  dims=dimsizes(sstannual)
  confi1_c=student_t(rc1_c@tval,rc1_c@nptxy-2)
  confi1_c:=onedtond(confi1_c,(/dims(1),dims(2)/))
  rc1_c@long_name = "Trend"
  copy_VarMeta(sstannual(0,:,:),rc1_c)
  copy_VarMeta(rc1_c,confi1_c)

  confi2_c=student_t(rc2_c@tval,rc2_c@nptxy-2)
  confi2_c:=onedtond(confi2_c,(/dims(1),dims(2)/))
  rc2_c@long_name = "Trend"
  copy_VarMeta(sstannual(0,:,:),rc2_c)
  copy_VarMeta(rc2_c,confi2_c)

  rc1_c_c:= rc1_c*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
  rc2_c_c:= rc2_c*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)

printVarSummary(sstannual)
copy_VarMeta(sstannual(0,:,:),rc1_c)
copy_VarMeta(sstannual(0,:,:),rc2_c)

rc1_all_c=new((/111,181,360/),float);;;89,180
rc2_all_c=new((/90,181,360/),float)
patterncor1_c=new((/153/),float)
patterncor1_a=new((/153/),float)
patterncor1_b=new((/153/),float)

patterncor2_c=new((/153/),float)
patterncor2_a=new((/153/),float)
patterncor2_b=new((/153/),float)


lats=-30
late=30
lons=120
lone=290
xx1=year1
do i=0,110
  rc1_all_c(i,:,:)=regCoef_n(year1,sstannual(0+i:42+i,:,:),0,0)
  copy_VarMeta(sstannual(0,:,:),rc1_all_c(i,:,:))
  patterncor1_c(i+42)=pattern_cor(rc1_c({lats:late},{lons:lone}),rc1_all_c(i,{lats:late},{lons:lone}), 1.0, 0)
  ;patterncor1_a(i+42)=pattern_cor(rc1({lats+30:late-30},{lons:lone}),rc1_all(i,{lats+30:late-30},{lons:lone}), 1.0, 0)
  ;patterncor1_b(i+42)=pattern_cor(rc1({lats+50:late-50},{lons:lone}),rc1_all(i,{lats+50:late-50},{lons:lone}), 1.0, 0)

end do

xx2=year2
do i=0,88
  rc2_all_c(i,:,:)=regCoef_n(year2,sstannual(0+i:64+i,:,:),0,0)
  copy_VarMeta(sstannual(0,:,:),rc2_all_c(i,:,:))
  patterncor2_c(i+64)=pattern_cor(rc2_c({lats:late},{lons:lone}),rc2_all_c(i,{lats:late},{lons:lone}), 1.0, 0)
  ;patterncor2_a(i+64)=pattern_cor(rc2({lats+30:late-30},{lons:lone}),rc2_all(i,{lats+30:late-30},{lons:lone}), 1.0, 0)
  ;patterncor2_b(i+64)=pattern_cor(rc2({lats+50:late-50},{lons:lone}),rc2_all(i,{lats+50:late-50},{lons:lone}), 1.0, 0)
end do

print(patterncor1_c)
print(patterncor2_c)



a = addfile(opendir+"sst.mon.anom.kaplan.remap.nc", "r")

time1      = a->time
ymd1      = cd_calendar(time1, -1)
iStrt3    = ind(ymd1.eq.ymdStrt1)
iLast3    = ind(ymd1.eq.ymdLast1)
delete(time1)
delete(ymd1)
ha_a   = a->sst(iStrt3:iLast3,::-1,:)

sstannual_a=month_to_annual(ha_a,1)
printVarSummary(sstannual_a)


syear1=1980
syear2=1958
year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)



rc1_a    = regCoef_n(year1,sstannual_a(syear1-syear:eyear-syear,:,:),0,0)
rc2_a    = regCoef_n(year2,sstannual_a(syear2-syear:eyear-syear,:,:),0,0)

confi1_a=student_t(rc1_a@tval,rc1_a@nptxy-2)
confi1_a:=onedtond(confi1_a,(/dims(1),dims(2)/))
rc1_a@long_name = "Trend"
copy_VarMeta(sstannual(0,:,:),confi1_a)

confi2_a=student_t(rc2_a@tval,rc2_a@nptxy-2)
confi2_a:=onedtond(confi2_a,(/dims(1),dims(2)/))
rc2@long_name = "Trend"
copy_VarMeta(sstannual(0,:,:),confi2_a)


  rc1_a:= rc1_a*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
  rc2_a:= rc2_a*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)

copy_VarMeta(sstannual_a(0,:,:),rc1_a)
copy_VarMeta(sstannual_a(0,:,:),rc2_a)

rc1_all_a=new((/111,181,360/),float);36,72
rc2_all_a =new((/90,181,360/),float)



lons=120
lone=290

do i=0,110
  rc1_all_a(i,:,:)=regCoef_n(year1,sstannual_a(0+i:42+i,:,:),0,0)
  copy_VarMeta(sstannual_a(0,:,:),rc1_all_a(i,:,:))
  patterncor1_a(i+42)=pattern_cor(rc1_a({lats:late},{lons:lone}),rc1_all_a(i,{lats:late},{lons:lone}), 1.0, 0)
end do

xx2=year2
do i=0,88
  rc2_all_a(i,:,:)=regCoef_n(year2,sstannual_a(0+i:64+i,:,:),0,0)
  copy_VarMeta(sstannual_a(0,:,:),rc2_all_a(i,:,:))
  patterncor2_a(i+64)=pattern_cor(rc2_a({lats:late},{lons:lone}),rc2_all_a(i,{lats:late},{lons:lone}), 1.0, 0)
end do

print(patterncor1_a)
print(patterncor2_a)



a= addfile(opendir+"sst.mon.mean.cobe.remap.nc","r")
syear    = 1891
eyear    = 2022
ymdStrt1 = syear*100+1
ymdLast1 = eyear*100+12
time1      = a->time
ymd1      = cd_calendar(time1, -1)
iStrt3    = ind(ymd1.eq.ymdStrt1)
iLast3    = ind(ymd1.eq.ymdLast1)
delete(time1)
delete(ymd1)

ha_b   = a->sst(iStrt3:iLast3,::-1,:)

delete(sstannual)
sstannual=month_to_annual(ha_b,1)
printVarSummary(sstannual)

printVarSummary(sstannual)
sstannual_b=new((/153,181,360/),float)
sstannual_b(21::,:,:)=sstannual

syear1=1980
syear2=1958
year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)

syear    = 1870

rc1_b    = regCoef_n(year1,sstannual_b(syear1-syear:eyear-syear,:,:),0,0)
rc2_b    = regCoef_n(year2,sstannual_b(syear2-syear:eyear-syear,:,:),0,0)


  confi1_b=student_t(rc1_b@tval,rc1_b@nptxy-2)
  confi1_b:=onedtond(confi1_b,(/dims(1),dims(2)/))
  rc1_b@long_name = "Trend"
  copy_VarMeta(sstannual(0,:,:),rc1_b)
  copy_VarMeta(rc1_b,confi1_b)

  confi2_b=student_t(rc2_b@tval,rc2_b@nptxy-2)
  confi2_b:=onedtond(confi2_b,(/dims(1),dims(2)/))
  rc2_b@long_name = "Trend"
  copy_VarMeta(sstannual(0,:,:),rc2_b)
  copy_VarMeta(rc2_b,confi2_b)



    rc1_b:= rc1_b*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
    rc2_b:= rc2_b*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)



printVarSummary(sstannual)
copy_VarMeta(sstannual_b(0,:,:),rc1_b)
copy_VarMeta(sstannual_b(0,:,:),rc2_b)

rc1_all_b=new((/111,181,360/),float)
rc2_all_b=new((/90,181,360/),float)


lats=-30
late=30
lons=120
lone=290
xx1=year1
do i=21,110
  rc1_all_b(i,:,:)=regCoef_n(year1,sstannual_b(0+i:42+i,:,:),0,0)
  copy_VarMeta(sstannual_b(0,:,:),rc1_all_b(i,:,:))
  patterncor1_b(i+42)=pattern_cor(rc1_b({lats:late},{lons:lone}),rc1_all_b(i,{lats:late},{lons:lone}), 1.0, 0)
  ;patterncor1_a(i+42)=pattern_cor(rc1({lats+30:late-30},{lons:lone}),rc1_all(i,{lats+30:late-30},{lons:lone}), 1.0, 0)
  ;patterncor1_b(i+42)=pattern_cor(rc1({lats+50:late-50},{lons:lone}),rc1_all(i,{lats+50:late-50},{lons:lone}), 1.0, 0)

end do

xx2=year2
do i=21,88
  rc2_all_b(i,:,:)=regCoef_n(year2,sstannual_b(0+i:64+i,:,:),0,0)
  copy_VarMeta(sstannual_b(0,:,:),rc2_all_b(i,:,:))
  patterncor2_b(i+64)=pattern_cor(rc2_b({lats:late},{lons:lone}),rc2_all_b(i,{lats:late},{lons:lone}), 1.0, 0)
  ;patterncor2_a(i+64)=pattern_cor(rc2({lats+30:late-30},{lons:lone}),rc2_all(i,{lats+30:late-30},{lons:lone}), 1.0, 0)
  ;patterncor2_b(i+64)=pattern_cor(rc2({lats+50:late-50},{lons:lone}),rc2_all(i,{lats+50:late-50},{lons:lone}), 1.0, 0)
end do





resultdir="/data0/fjiang/scripts/contrast/revise-01/"


  wks = gsn_open_wks("eps",resultdir+"r01-F1");1979-2018")
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
   res@cnFillColors=(/12,20,30,40,50,57,0,0,73,80,90,96,105,115/)

   res@gsnLeftStringFontHeightF=0.01
   res@gsnLeftStringOrthogonalPosF=0.02;0.005
   res@gsnRightStringFontHeightF=0.012
   res@gsnRightStringOrthogonalPosF=-0.003;5;0.005
   res@cnLineLabelsOn=False
   res@gsnRightString="";(~S~o~N~C/decade)"
   res@gsnLeftString=""
   res@tiMainString="(a) HadISST SST Trend (1980-2022)  "
   res@tiMainOffsetYF=-0.01
   res@tiMainFontHeightF=0.011
   plot1=new(8,graphic)
   ;dum1=new(4,graphic)
   res@vpXF=-0.119
   res@vpYF=0.97
   res@vpHeightF= 0.165           ; change aspect ratio of plot
   res@cnLevels             = 0.02*ispan(-12,12,2 )

   res@lbLabelBarOn=False

   plot1(0) = gsn_csm_contour_map(wks,rc1,res)
   res@lbLabelBarOn=True
   res@lbOrientation="Vertical"
   ;res@lbBoxEndCapStyle="TriangleBothEnds"
   res@pmLabelBarHeightF=0.16
   res@pmLabelBarWidthF=0.03
   res@lbLabelStride=4.
   res@lbBoxLinesOn=True
   res@pmLabelBarOrthogonalPosF=0.05
   res@lbLabelFontHeightF=0.01
   ;res@pmLabelBarOrthogonalPosF=0.2
   ;res@pmLabelBarParallelPosF=-0.11
   res@tiMainString="(b) HadISST SST Trend (1958-2022) "
   ;res@cnLevels             = 0.02*ispan(-10,10,1 )
   res@vpXF=0.397-0.17
   res@vpYF=0.97
   res@pmLabelBarHeightF=0.04
   res@pmLabelBarWidthF=0.4
   res@lbLabelStride=4.
   res@pmLabelBarParallelPosF=-0.03
   res@pmLabelBarOrthogonalPosF=0.02
   res@lbOrientation="Horizontal"

   plot1(1) = gsn_csm_contour_map(wks,rc2,res)
   res@gsnLeftStringOrthogonalPosF=0.02;0.005


  delete(res)
  res                 = True
  res@gsnDraw         = False
  res@gsnFrame        = False
  res@vpHeightF = 0.07
  res@vpWidthF  = 0.2
  res@trXMinF   =  1870.
  res@trXMaxF   =  2022.
  res@trYMaxF   = 2.;27.5
  res@trYMinF   = -2.;25.5
  res@tiYAxisString = " "
  res@tiXAxisString = " "
  res@tiMainString  = " "
  res@tmXBMajorThicknessF     = 1.5
  res@tmYLMajorThicknessF     = 1.5
  res@tmXBMinorThicknessF     = 1.0
  res@tmYLMinorThicknessF     = 1.0
  res@tmYLMinorOn=False
  res@tmXTOn=False

  res@tmXBValues    = ispan(1870,2020,30)
  res@tmXBLabels    = ispan(1870,2020,30);(/"1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015","2020"/)
  res@tmXBLabelsOn=True
  res@tmXBMode="ExplicitLevels"
  res@tmXBMajorLengthF=0.006
  res@tmYLMajorLengthF=0.004
  res@tmXBMajorOutwardLengthF=0.006
  res@tmYLMajorOutwardLengthF=0.004
  res@tmXBMinorLengthF        = 0.004
  res@tmYLMinorLengthF        = 0.004
  res@tmXBMinorOn=True
  ;res@tmXBMinorValues=ispan(1958,2018,1)

  res@gsnLeftString=" "
  res@gsnLeftStringFontHeightF= 0.014
  res@gsnLeftStringOrthogonalPosF=0.00;5
  ;res@tmXBOn = False
  ;res@tmYLOn = True
  ;res@tmYLLabelsOn = True
  ;res@tmYLMinorOn = True
  ;res@tmYLLabelFontColor = "black"
  res@tiYAxisString         = "";Tropical Pacific SST (~S~o~N~C)"
  res@tiYAxisFontColor      = "black"
  res@tmYROn    = False;True
  res@tmYRLabelsOn = False
  res@tmYLOn    = True
  res@tmYLMode          = "Explicit"         ; explicit labels
  res@tmYLLabelsOn=True
  res@tmYLValues        = 1.*ispan(-5,5,1);ispan(51,55,1)
  res@tmYLLabels        = 1.*ispan(-5,5,1);;(/"25.5","26.0","26.5","27.0","27.5"/)
  res@tmXBLabelFontHeightF  = 0.0075
  res@tmYLLabelFontHeightF  = 0.0075
  res@tiYAxisSide           = "Left"
  res@tiYAxisFontHeightF    = 0.0075
  res@xyDashPattern    = 2
  res@xyLineThicknessF = 1.
  res@xyLineColor      = "black"
  ;res@tmXBBorderOn=False
  res@tiYAxisFontHeightF=0.01
  res@gsnYRefLine=0.

  res_low=res
  res_low@xyLineThicknessF=0.8
  res_low@xyDashPattern    = 0

  res@xyLineOpacityF=0.8

  res_trend=res_low
  res_trend@xyLineThicknessF=2.
  res@tmXTBorderOn=False
  res@tmYRBorderOn=False
  res@tmXBBorderOn=False


  plot1a = new(3,graphic)
  plot2a = new(3,graphic)
  res@gsnLeftString=""
  res@vpXF=0.39+0.375
  res@vpYF=0.962
  res@tiMainString="(c) Equatorial Pacific SST Evolution"
  res@tiMainFont=22
  res@tiMainOffsetYF=-0.005;2
  res@tiMainFontHeightF=0.010
  res@tmXBLabelsOn=False;(/"1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015","2020"/)
  res@tmXBOn=False
  res@tiYAxisOffsetXF=-0.008

  res@tiYAxisString="Cold Tongue"
  plot1a(0) = gsn_csm_xy(wks,ispan(1870,2022,1),sstannual_ep,res)
  res@tiYAxisOffsetXF=0.0004;5

  res@tiMainString="";Warm Pool SST Anomalies "
  res@vpYF=0.874
  res@trYMaxF   = 1.;27.5
  res@trYMinF   = -1.;25.5
  res@tmXBLabelsOn=True;(/"1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015","2020"/)
  res@tmYLLabels        = 0.5*ispan(-5,5,1);;(/"25.5","26.0","26.5","27.0","27.5"/)
  res@tmYLValues        = .5*ispan(-5,5,1);ispan(51,55,1)
  res@tmXBBorderOn=True
  res@tmXBOn=True
  res@tiXAxisString         = "Year"
  res@tiYAxisString="Warm Pool"
  plot1a(1) = gsn_csm_xy(wks,ispan(1870,2022,1),sstannual_cp,res)
  res@gsnLeftString="(c) Cold tongue minus warm pool SST anomalies"

  ;plot1a(2) = gsn_csm_xy(wks,ispan(1870,2022,1),sstannual_gr,res)
  res@tiXAxisString         = ""

  res@gsnLeftString=""
  res@xyDashPattern    = 0
  res@xyLineThicknessF = 3.
  res@xyLineColor      = "orangered3"

  plot2a(0) = gsn_csm_xy(wks,ispan(1870,2022,1),sstannual_ep_10,res)
  plot2a(1) = gsn_csm_xy(wks,ispan(1870,2022,1),sstannual_cp_10,res)
  ;plot2a(2) = gsn_csm_xy(wks,ispan(1870,2022,1),sstannual_gr_10,res)

  do i=0,1
   overlay(plot1a(i),plot2a(i))
  end do







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
  plot2(0)=gsn_csm_contour(wks,confi1,sres1)
  plot2(1)=gsn_csm_contour(wks,confi2,sres1)





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
res1                 = True
res1@gsnDraw         = False
res1@gsnFrame        = False
res1@vpHeightF = 0.18
res1@vpWidthF  = 0.44
res1@trXMinF   =  1870.
res1@trXMaxF   =  2022.
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
res1@tmXBValues    = ispan(1870,2020,30)
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
res1@tiXAxisString         = "End year of the trend"
res1@tiXAxisFontHeightF=0.01
res1@tiYAxisFontColor      = "black"
res1@tmYLOn    = True
res1@tmYLMode          = "Explicit"         ; explicit labels
res1@tmYLLabelsOn=True

res1@tmXBLabelFontHeightF  = 0.008
res1@tmYLLabelFontHeightF  = 0.008
res1@tiYAxisSide           = "Left"
res1@tiYAxisFontHeightF    = 0.01
res1@xyDashPattern    = 2
res1@tmXBBorderOn=True
res1@tiYAxisFontHeightF=0.012

res1@tmXTOn=False
res1@tmYROn=False
res1@gsnLeftString=""
res1@tiMainString="(d) Pattern Correlation with Tropical Pacific Trend 1980-2022"
res1@tiMainOffsetYF=-0.01
res1@tiMainFontHeightF=0.0115
res1@gsnRightString=""
res1@xyLineThicknessF = 2.
res1@vpXF=0.04
res1@vpYF=0.71
res1@tiMainFont=22
res1@xyLineColor      = "gray30"
res1@xyDashPattern    = 0

plot1(2) = gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1,res1)
res1@tiMainString="(e) Pattern Correlation with Tropical Pacific Trend 1958-2022"
res1@vpXF=0.54

plot1(3) = gsn_csm_xy(wks,ispan(1870,2022,1),patterncor2,res1)
res1@xyDashPattern    = 1
res1@xyLineColor      = "orangered3"
res1@gsnLeftString=""
res1@xyLineThicknessF = 1.5

plota=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1_a,res1)
res1@xyLineColor      = "deepskyblue3"
res1@xyDashPattern    = 2
plotb=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1_b,res1)
res1@xyDashPattern    = 1
res1@xyLineColor      = "orangered3"
plotc=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor2_a,res1)
res1@xyLineColor      = "deepskyblue3"
res1@xyDashPattern    = 2
plotd=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor2_b,res1)
res1@xyDashPattern    = 1
res1@xyLineThicknessF = 1.

res1@xyLineColor      = "darkseagreen"
plote=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1_c,res1)
res1@xyLineColor      = "darkseagreen"
res1@xyDashPattern    = 1
plotf=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor2_c,res1)
overlay(plot1(2),plota)
overlay(plot1(2),plotb)
overlay(plot1(3),plotc)
overlay(plot1(3),plotd)
overlay(plot1(2),plote)
overlay(plot1(3),plotf)


res1@gsnAboveYRefLineColor = "gray80"             ; above ref line fill red
res1@gsnBelowYRefLineColor = "gray80"
res1@xyLineColor      = "gray80"

plotg=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor1_ipo,res1)
ploth=gsn_csm_xy(wks,ispan(1870,2022,1),patterncor2_ipo,res1)
overlay(plot1(2),plotg)
;overlay(plot1(3),ploth)

pgres=True
pgres@gsLineColor      = "black"
pgres@gsLineThicknessF = 2.0
pgres@tfPolyDrawOrder ="PostDraw"
pgres@xyDashPattern    = 0


print(patterncor1_a)
lat1s=-5
lat1n=5
lon1e=140
lon1w=170
lon2e=190
lon2w=270

  overlay(plot1(0),plot2(0))
  overlay(plot1(1),plot2(1))
dum_box1 = gsn_add_polyline   (wks, plot2(1),(/lon1e,lon1w,lon1w,lon1e,lon1e/),(/lat1s,lat1s,lat1n,lat1n,lat1s/), pgres)
dum_box2 = gsn_add_polyline   (wks, plot2(1),(/lon2e,lon2w,lon2w,lon2e,lon2e/),(/lat1s,lat1s,lat1n,lat1n,lat1s/), pgres)


  draw(plot1(0))
  draw(plot1(1))



;pgres@gsLineColor      = "gray40"
;pgres@gsLineThicknessF = 1.5
;pgres@gsLineDashPattern    = 2
;
;dum_box1_a = gsn_add_polyline (wks, plot2(0),(/lons,lone,lone,lons,lons/),(/lats+30,lats+30,late-30,late-30,lats+30/), pgres)
;pgres@gsLineColor      = "deepskyblue3"
;pgres@gsLineDashPattern    = 1
;
;dum_box1_b = gsn_add_polyline (wks, plot2(0),(/lons,lone,lone,lons,lons/),(/lats+50,lats+50,late-50,late-50,lats+50/), pgres)
;pgres@gsLineColor      = "gray40"
;pgres@gsLineDashPattern    = 2
;
;dum_box2_a = gsn_add_polyline (wks, plot2(1),(/lons,lone,lone,lons,lons/),(/lats+30,lats+30,late-30,late-30,lats+30/), pgres)
;pgres@gsLineColor      = "deepskyblue3"
;pgres@gsLineDashPattern    = 1
;
;dum_box2_b = gsn_add_polyline (wks, plot2(1),(/lons,lone,lone,lons,lons/),(/lats+50,lats+50,late-50,late-50,lats+50/), pgres)


txres=True
txres@txFontHeightF=0.0105
;gsn_text_ndc(wks, "[~S~o~N~C/decade]",0.955,0.79, txres)
gsn_text_ndc(wks, "[~S~o~N~C/decade]",0.35,0.755, txres)

;gsn_text_ndc(wks, "[~S~o~N~C/decade]",0.935,0.271, txres)
;gsn_text_ndc(wks, "[~S~o~N~C/decade]",0.935,0.031, txres)





pgres=True
pgres@gsMarkerIndex      = 12

pgres@gsMarkerSizeF = 4.
pgres@tfPolyDrawOrder ="PostDraw"
pgres@gsMarkerColor="black"
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
;
pgres@gsMarkerColor="gray30"
pgres@gsMarkerIndex      = 16
pgres@gsMarkerSizeF = 8.

dum_gon_7 = gsn_add_polymarker (wks,  plot1(2),1870+42, patterncor1(42), pgres)
dum_gon_8 = gsn_add_polymarker (wks,  plot1(2),1870+86, patterncor1(86), pgres)
dum_gon_9 = gsn_add_polymarker (wks,  plot1(2),1870+114,patterncor1(114), pgres)
dum_gon_10 = gsn_add_polymarker (wks, plot1(3),1870+64, patterncor2(64), pgres)

print(patterncor1)
print(patterncor2)
txres@txFontHeightF=0.0105

tx_a = gsn_add_text (wks, plot1(2),"P1",1870+42-20,patterncor1(42)+0.12,   txres)
tx_b = gsn_add_text (wks, plot1(2),"P2",1870+86-20,  patterncor1(86)+0.12,   txres)
tx_c = gsn_add_text (wks, plot1(2),"P3",1870+114-20, patterncor1(114)+0.12, txres)
tx_d = gsn_add_text (wks, plot1(3),"P4",1870+64-32,patterncor2(64)+0.12,   txres)



   lgres=True
   amres=True
   lgres@lgLineColors=(/"deepskyblue3","orangered3","darkseagreen","gray30"/)
   lgres@lgLeftMarginF=0.1
   lgres@lgRightMarginF=0.1

   lgres@lgDashIndexes=(/1,2,1,0/)
   lgres@lgLineThicknesses=(/1.5,1.5,1.5,2.5/)
   lgres@lgMonoLineThickness=False
   lgres@vpWidthF=0.097
   lgres@vpHeightF=0.07
   lgres@lgLabelFontHeightF=0.12
   lgres@lgLineDashSegLenF=0.04
   lgres@lgPerimThicknessF=2.
   lgres@lgPerimOn=False
   lgres@lgPerimFill="SolidFill"
   lgres@lgPerimFillColor="white"
   ;lgres@lgLabelOffsetF=-0.1
   lgres@lbBoxMajorExtentF=0.3
   lbid1=gsn_create_legend(wks,4,(/"  COBE","  Kaplan","  ERSSTv5","  HadISST"/),lgres)
   lbid2=gsn_create_legend(wks,4,(/"  COBE","  Kaplan","  ERSSTv5","  HadISST"/),lgres)

   amres@amParallelPosF=-0.366
   amres@amOrthogonalPosF=-0.28
   annoid1=gsn_add_annotation(plot1(2),lbid1,amres)
   annoid2=gsn_add_annotation(plot1(3),lbid2,amres)



   lgres@lgLineColors:=(/"black"/)
   lgres@lgDashIndexes:=(/2/)
   lgres@lgLineThicknesses:=(/1.5/)
   lgres@vpWidthF=0.05
   lgres@vpHeightF=0.03

   lbid3=gsn_create_legend(wks,1,(/"  Raw"/),lgres)
   lgres@lgLineColors:=(/"orangered3"/)
   lgres@lgDashIndexes:=(/0/)
   lgres@lgLineThicknesses:=(/2.5/)
   lgres@vpWidthF=0.11
   lbid4=gsn_create_legend(wks,1,(/"  15-yr running mean"/),lgres)
   amres@amParallelPosF=-0.31
   amres@amOrthogonalPosF=0.3
   annoid3=gsn_add_annotation(plot1a(1),lbid3,amres)
   amres@amParallelPosF=0.129

   annoid4=gsn_add_annotation(plot1a(1),lbid4,amres)



   draw(plot1a(0))
   draw(plot1a(1))


   ;---Values for three arrows.
    xarrow1a = (/1891,1912/)
    yarrow1a = (/patterncor1(42),patterncor1(42)/)
    xarrow1b = (/1891,1870/)
    yarrow1b = (/patterncor1(42),patterncor1(42)/)

    xarrow2a = (/1935,1956/)
    yarrow2a = (/patterncor1(86),patterncor1(86)/)
    xarrow2b = (/1935,1914/)
    yarrow2b = (/patterncor1(86),patterncor1(86)/)

    xarrow3a = (/1963,1984/)
    yarrow3a = (/patterncor1(114),patterncor1(114)/)
    xarrow3b = (/1963,1942/)
    yarrow3b = (/patterncor1(114),patterncor1(114)/)


    xarrow4a = (/1902,1934/)
    yarrow4a = (/patterncor2(64),patterncor2(64)/)
    xarrow4b = (/1902,1870/)
    yarrow4b = (/patterncor2(64),patterncor2(64)/)



    resp = True   ; Resource list for customizing arrows. Start with default arrows.
    resp@gsLineColor       = "black"
    resp@gsLineDashPattern = 0    ; solid line
    resp@gsLineThicknessF  = 2.0  ; 3x thicker
    plot1(2) = arrow(wks,plot1(2),xarrow1a,yarrow1a,resp)
    plot1(2) = arrow(wks,plot1(2),xarrow1b,yarrow1b,resp)
    plot1(2) = arrow(wks,plot1(2),xarrow2a,yarrow2a,resp)
    plot1(2) = arrow(wks,plot1(2),xarrow2b,yarrow2b,resp)
    plot1(2) = arrow(wks,plot1(2),xarrow3a,yarrow3a,resp)
    plot1(2) = arrow(wks,plot1(2),xarrow3b,yarrow3b,resp)

    plot1(3) = arrow(wks,plot1(3),xarrow4a,yarrow4a,resp)
    plot1(3) = arrow(wks,plot1(3),xarrow4b,yarrow4b,resp)


   draw(plot1(2))
   draw(plot1(3))



      res@lbLabelBarOn=False
      aa=-rc1_all(42-42,:,:)
      copy_VarMeta(rc1_all(0,:,:),aa)
      res@vpXF=-0.04
      res@vpYF=0.47
      res@gsnLeftString=""
      res@tiMainString="(e) HadISST SST Trend "+"~F34~4 "+"~F22~-1 P1 (1870-1912)  "

      res@tiMainOffsetYF=-0.01
      res@tiMainFontHeightF=0.0115
      plot1(4) = gsn_csm_contour_map(wks,aa,res)

      res@tiMainString="(f) HadISST SST Trend P2 (1904-1956)  "
      ;res@cnLevels             = 0.02*ispan(-10,10,1 )
      res@vpXF=0.33
      res@vpYF=0.47
      res@lbLabelBarOn=True

      plot1(5) = gsn_csm_contour_map(wks,rc1_all(86-42,:,:),res)

      res@tiMainString="(g)  HadISST SST Trend "+"~F34~4 "+"~F22~-1 P3 (1942-1984)  "
      bb=-rc1_all(117-42,:,:)
      copy_VarMeta(rc1_all(0,:,:),bb)
      res@vpXF=-0.04
      res@vpYF=0.23
      res@lbLabelBarOn=False

      plot1(6) = gsn_csm_contour_map(wks,bb,res)
      res@vpXF=0.33
      res@vpYF=0.23

      res@lbLabelBarOn=True
      ;res@lbBoxEndCapStyle="TriangleBothEnds"
      ;res@pmLabelBarParallelPosF=1.2
      res@tiMainString="(h)  HadISST SST Trend "+"~F34~4 "+"~F22~-1 P4 (1870-1934) "
      cc=-rc2_all(64-64,:,:)
      copy_VarMeta(rc1_all(0,:,:),cc)
      ;res@tiMainOffsetYF=-0.03


      plot1(7) = gsn_csm_contour_map(wks,cc,res)



      plot2(2)=gsn_csm_contour(wks,confi1_all(42-42,:,:), sres1)
      plot2(3)=gsn_csm_contour(wks,confi1_all(86-42,:,:), sres1)
      plot2(4)=gsn_csm_contour(wks,confi1_all(117-42,:,:),sres1)
      plot2(5)=gsn_csm_contour(wks,confi2_all(64-64,:,:), sres1)



      overlay(plot1(4),plot2(2))
      overlay(plot1(5),plot2(3))
      overlay(plot1(6),plot2(4))
      overlay(plot1(7),plot2(5))


;draw(plot1(4))
;draw(plot1(5))
;draw(plot1(6))
;draw(plot1(7))


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
